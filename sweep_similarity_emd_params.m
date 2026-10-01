% SWEEP_SIMILARITY_EMD_PARAMS
%
% With the new compute_EMD (R eliminated, constraint exact via correct_warp,
% lambdaR entering only as proj_linf(lambdaR)), does similarity_WH_EMD still
% need LambdaR=1e3 and Tol=1e-6, or can both be relaxed for speed?
%
% The old soft-constraint formulation amplified solver residual by lambdaR, so
% Tol had to track LambdaR. That coupling is gone. What remains:
%
%   * LambdaR still sets the tradeoff between transporting and discarding mass.
%     With unit-mass profiles of width L, LambdaR of order L is the natural
%     scale (pay for displacing a unit of mass up to ~L bins before preferring
%     to drop it).
%   * Tol still controls how accurately M itself is recovered.
%
% This script measures, on small synthetic pairs that mirror what the matcher
% sees after unit-mass normalisation:
%
%   1. relative error in ||M||_1 against a tight reference solve
%   2. wall time and TFOCS iterations
%   3. separation of a true (shifted) pair from a wrong (disjoint-support) pair
%      via the matcher's reference-cost ratio

clear
close all

root = fileparts(pwd);
addpath(fullfile(root, 'TFOCS'))
addpath(fullfile(root, 'Utils'))
addpath(genpath(fullfile(root, 'FlexMF')))

%% Synthetic unit-mass profiles, same scale as generate_data motifs
N = 15;
L = 35;
rng(0);

% Motif A: three neurons, staggered peaks (disjoint support from B)
WA = zeros(N, L);
WA(1, 5:7) = 1;
WA(2, 10:12) = 1;
WA(3, 15:17) = 1;
WA = WA / sum(WA(:));

% Motif B: different neurons, different lags
WB = zeros(N, L);
WB(8, 8:10) = 1;
WB(9, 16:18) = 1;
WB(10, 24:26) = 1;
WB = WB / sum(WB(:));

% Exact copy shifted by 3 bins (transport cost should be ~3)
WAshift = helper.shift_profiles(WA, 3, L);
WAshift = WAshift / max(sum(WAshift(:)), eps);  % shift drops edge mass; renormalise

% Mild amplitude perturbation on the shifted copy
WApert = WAshift;
WApert(WApert > 0) = WApert(WApert > 0) .* (1 + 0.1 * randn(nnz(WApert), 1));
WApert = max(WApert, 0);
WApert = WApert / sum(WApert(:));

cases = {
    'identical',     WA,      WA
    'shift3',        WA,      WAshift
    'shift3+pert',   WA,      WApert
    'wrong',         WA,      WB
    'empty-ref',     zeros(N,L), WA
    };

lambdaRs = [10, 35, 70, 100, 200, 1000];
tols     = [1e-3, 1e-4, 1e-5, 1e-6];

%% Reference solves: same lambdaR, very tight tol
fprintf('Building reference solves (tol=1e-8)...\n');
ref = struct();
for r = 1:numel(lambdaRs)
    opts = make_opts(1e-8);
    cont = make_cont();
    for c = 1:size(cases, 1)
        [d, M, R, out] = compute_EMD(cases{c, 2}, cases{c, 3}, opts, ...
            'continuationOptions', cont, 'lambdaR', lambdaRs(r));
        ref(r, c).d = d;
        ref(r, c).transport = norm(M(:), 1);
        ref(r, c).residual = lambdaRs(r) * norm(R(:), 1);
        ref(r, c).niter = niter_of(out);
    end
end

%% Sweep
fprintf('\n%8s %8s %14s %10s %10s %8s %8s %8s\n', ...
    'lambdaR', 'tol', 'case', 'rel_err_M', 'rel_err_d', 'niter', 'sec', 'cost/ref');
fprintf('%s\n', repmat('-', 1, 90));

rows = [];
for r = 1:numel(lambdaRs)
    for t = 1:numel(tols)
        opts = make_opts(tols(t));
        cont = make_cont();

        % empty-ref cost for this (lambdaR, tol)
        [dRef, ~, ~, ~] = compute_EMD(zeros(N, L), WA, opts, ...
            'continuationOptions', cont, 'lambdaR', lambdaRs(r));

        for c = 1:size(cases, 1)
            tic
            [d, M, R, out] = compute_EMD(cases{c, 2}, cases{c, 3}, opts, ...
                'continuationOptions', cont, 'lambdaR', lambdaRs(r));
            sec = toc;
            transport = norm(M(:), 1);
            niter = niter_of(out);

            denomM = max(ref(r, c).transport, 1e-12);
            denomD = max(ref(r, c).d, 1e-12);
            relM = abs(transport - ref(r, c).transport) / denomM;
            relD = abs(d - ref(r, c).d) / denomD;
            % identical short-circuits to 0; treat as exact
            if strcmp(cases{c, 1}, 'identical')
                relM = 0;
                relD = 0;
            end
            ratio = d / max(dRef, eps);

            fprintf('%8g %8.0e %14s %10.3g %10.3g %8d %8.3f %8.3f\n', ...
                lambdaRs(r), tols(t), cases{c, 1}, relM, relD, niter, sec, ratio);

            rows = [rows; struct( ...
                'lambdaR', lambdaRs(r), 'tol', tols(t), ...
                'case', string(cases{c, 1}), ...
                'relM', relM, 'relD', relD, ...
                'transport', transport, 'd', d, ...
                'refTransport', ref(r, c).transport, 'refD', ref(r, c).d, ...
                'niter', niter, 'sec', sec, 'ratio', ratio, 'dRef', dRef)]; %#ok<AGROW>
        end
    end
end

T = struct2table(rows);
save('sweep_similarity_emd_params.mat', 'T', 'lambdaRs', 'tols', 'cases', 'ref');

%% Summaries that answer the question
fprintf('\n=== worst relative ||M||_1 error over non-identical cases ===\n');
fprintf('%8s %8s %12s %10s\n', 'lambdaR', 'tol', 'worst_relM', 'median_sec');
for r = 1:numel(lambdaRs)
    for t = 1:numel(tols)
        mask = T.lambdaR == lambdaRs(r) & T.tol == tols(t) & T.case ~= "identical";
        fprintf('%8g %8.0e %12.3g %10.3f\n', ...
            lambdaRs(r), tols(t), max(T.relM(mask)), median(T.sec(mask)));
    end
end

fprintf('\n=== true/wrong separation via cost/ref (needs true << 0.7 << wrong) ===\n');
fprintf('%8s %8s %12s %12s %12s %s\n', ...
    'lambdaR', 'tol', 'shift3', 'shift+pert', 'wrong', 'gap?');
for r = 1:numel(lambdaRs)
    for t = 1:numel(tols)
        rShift = T.ratio(T.lambdaR == lambdaRs(r) & T.tol == tols(t) & T.case == "shift3");
        rPert  = T.ratio(T.lambdaR == lambdaRs(r) & T.tol == tols(t) & T.case == "shift3+pert");
        rWrong = T.ratio(T.lambdaR == lambdaRs(r) & T.tol == tols(t) & T.case == "wrong");
        ok = max(rShift, rPert) < 0.7 && rWrong > 0.7;
        fprintf('%8g %8.0e %12.3f %12.3f %12.3f %s\n', ...
            lambdaRs(r), tols(t), rShift, rPert, rWrong, ternary(ok, 'YES', 'NO'));
    end
end

fprintf('\n=== reference transport for shift3 (should be ~3) ===\n');
for r = 1:numel(lambdaRs)
    c = find(strcmp(cases(:, 1), 'shift3'));
    fprintf('  lambdaR=%g  ||M||_1 = %.4f  (ref tol=1e-8)\n', ...
        lambdaRs(r), ref(r, c).transport);
end

%% Local helpers
function opts = make_opts(tol)
opts = tfocs_SCD;
opts.continuation = 1;
opts.tol = tol;
opts.stopCrit = 4;
opts.maxIts = 100;
opts.printEvery = 0;
opts.alg = 'N83';
end

function cont = make_cont()
cont = continuation();
cont.verbose = 0;
end

function n = niter_of(out)
if isempty(out)
    n = 0;
elseif isfield(out, 'niter')
    n = out.niter;
elseif isfield(out, 'iterations')
    n = out.iterations;
else
    n = -1;
end
end

function s = ternary(cond, a, b)
if cond, s = a; else, s = b; end
end
