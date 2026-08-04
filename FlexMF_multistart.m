function [W, H, cost, errors, loadings, power, M, R, info] = FlexMF_multistart(X, varargin)
% FlexMF_multistart  Multi-restart FlexMF with optional SeqNMF warm-start
%
% USAGE:
%   [W, H, cost, errors, loadings, power, M, R, info] = FlexMF_multistart(X, ...
%       'K', K, 'L', L, 'EMD', 1, ...
%       'nRestarts', 10, 'warmStart', true, ...
%       'lambda', 0.1, 'lambda_M', 0.1, 'lambda_R', 1, ...);
%
% Multistart-specific options (all other name/value pairs are passed to FlexMF):
%   'nRestarts'       10       Number of independent FlexMF fits
%   'warmStart'       false    If true, each restart warm-starts FlexMF from a
%                              fresh SeqNMF fit (same K, L, lambda)
%   'seqNMF_maxiter'  50       Max iterations for each SeqNMF warm-start
%   'seqNMF_lambda'   []       SeqNMF lambda (default: FlexMF 'lambda')
%   'includeConstraint' false  If true and EMD=1, add ||Xcorr-R-W*H||_1 to the
%                              selection objective. Off by default: in EMD mode
%                              the reconstruction is a hard constraint solved
%                              only approximately by TFOCS, so its residual is
%                              solver slack rather than model quality. Adding it
%                              at unit weight can dominate the penalty terms and
%                              select the run that happened to converge tightest.
%                              Use info.constraints_rel to screen runs instead.
%   'verbal'          1        Print per-restart progress (FlexMF itself stays quiet)
%
% Best-run selection objective (EMD=1), matching FlexMF_test_init.m:
%   obj = lambda * reg_cross + lambda_M * ||M||_1 + lambda_R * ||R||_1
%        [+ ||Xcorr - R - reconstruct(W,H)||_1  if includeConstraint]
% For EMD=0:
%   obj = recon_err + lambda * reg_cross
%
% OUTPUTS:
%   W,H,cost,errors,loadings,power,M,R  — from the best FlexMF run
%   info  struct with fields:
%     .best_idx, .objs, .constraints, .constraints_rel, .warmStart
%     .nIters, .maxiter, .stoppedEarly   iterations run by each restart, and
%                                        whether it stopped before maxiter
%     .W_all, .H_all, .M_all, .R_all, .cost_all, .errors_all
%     .loadings_all, .power_all, .seqNMF_W, .seqNMF_H
%
% Before trusting the selected run, check that info.constraints_rel is small and
% comparable across restarts; a run with a much larger residual did not really
% solve the constrained problem and its objective is not comparable.
%
% See also: FlexMF, FlexMF_test_init, seqNMF

%% Parse multistart options; remaining args go to FlexMF
[ms, flexArgs] = parse_multistart_params(varargin);

nRestarts = ms.nRestarts;
warmStart = ms.warmStart;
seqNMF_maxiter = ms.seqNMF_maxiter;
includeConstraint = ms.includeConstraint;
verbal = ms.verbal;

% Pull FlexMF hyperparameters used in the selection objective / SeqNMF
[K, L, lambda, lambda_M, lambda_R, EMD, maxiter] = get_flex_defaults(flexArgs);
if isempty(ms.seqNMF_lambda)
    seqNMF_lambda = lambda;
else
    seqNMF_lambda = ms.seqNMF_lambda;
end

% FlexMF prints its own per-iteration diagnostics; keep it quiet and let the
% multistart 'verbal' flag control progress reporting
flexArgs = set_default_nv(flexArgs, 'showPlot', 0);
flexArgs = set_default_nv(flexArgs, 'verbal', 0);

% Drop any caller-provided W_init/H_init — restarts manage initialization
flexArgs = drop_nv(flexArgs, {'W_init', 'H_init'});

%% Allocate
objs = inf(nRestarts, 1);
constraints = nan(nRestarts, 1);
constraints_rel = nan(nRestarts, 1);
nIters = nan(nRestarts, 1);
normX1 = norm(X(:), 1);
W_all = cell(nRestarts, 1);
H_all = cell(nRestarts, 1);
M_all = cell(nRestarts, 1);
R_all = cell(nRestarts, 1);
cost_all = cell(nRestarts, 1);
errors_all = cell(nRestarts, 1);
loadings_all = cell(nRestarts, 1);
power_all = nan(nRestarts, 1);
seqNMF_W = cell(nRestarts, 1);
seqNMF_H = cell(nRestarts, 1);

%% Multi-restart loop
for n = 1:nRestarts
    if verbal
        fprintf('FlexMF_multistart: restart %d/%d', n, nRestarts);
        if warmStart
            fprintf(' (SeqNMF warm-start)');
        end
        fprintf('\n');
    end

    runArgs = flexArgs;
    if warmStart
        [W0, H0] = seqNMF(X, 'K', K, 'L', L, ...
            'lambda', seqNMF_lambda, 'maxiter', seqNMF_maxiter, ...
            'showPlot', 0);
        seqNMF_W{n} = W0;
        seqNMF_H{n} = H0;
        runArgs = [flexArgs, {'W_init', W0, 'H_init', H0}];
    end

    if EMD
        [Wn, Hn, costn, errorsn, loadingsn, powern, Mn, Rn] = FlexMF(X, runArgs{:});
    else
        [Wn, Hn, costn, errorsn, loadingsn, powern] = FlexMF(X, runArgs{:});
        Mn = [];
        Rn = [];
    end

    W_all{n} = Wn;
    H_all{n} = Hn;
    M_all{n} = Mn;
    R_all{n} = Rn;
    cost_all{n} = costn;
    errors_all{n} = errorsn;
    loadings_all{n} = loadingsn;
    power_all(n) = powern;

    % FlexMF trims cost to iter+1 entries when it stops, so the number of
    % iterations actually run is one less than the length of the cost trace
    nIters(n) = numel(costn) - 1;

    if EMD
        % Match FlexMF_test_init.m objective; optionally add constraint residual
        reg_cross = errorsn(end, 2);
        objs(n) = lambda * reg_cross + lambda_M * norm(Mn(:), 1) + lambda_R * norm(Rn(:), 1);
        Xcorr = helper.correct_warp(X, Mn);
        constraint = Xcorr - Rn - helper.reconstruct(Wn, Hn);
        constraints(n) = norm(constraint(:), 1);
        constraints_rel(n) = constraints(n) / normX1;
        if includeConstraint
            objs(n) = objs(n) + constraints(n);
        end
    else
        % Same form as the EMD branch: reconstruction + lambda * cross-orthogonality
        objs(n) = errorsn(end, 1) + lambda * errorsn(end, 2);
    end

    if verbal
        fprintf('  obj=%.6g', objs(n));
        if EMD
            fprintf(', constraint=%.6g (%.3g x ||X||_1)', constraints(n), constraints_rel(n));
        end
        fprintf(', iters=%d/%d', nIters(n), maxiter);
        if nIters(n) < maxiter
            fprintf(' (stopped early)');
        end
        fprintf('\n');
    end
end

%% Select best run
[~, best_idx] = min(objs);
W = W_all{best_idx};
H = H_all{best_idx};
cost = cost_all{best_idx};
errors = errors_all{best_idx};
loadings = loadings_all{best_idx};
power = power_all(best_idx);
M = M_all{best_idx};
R = R_all{best_idx};

if verbal
    fprintf('FlexMF_multistart: best restart = %d (obj=%.6g)\n', best_idx, objs(best_idx));
end

info = struct();
info.best_idx = best_idx;
info.objs = objs;
info.constraints = constraints;
info.constraints_rel = constraints_rel;
info.nIters = nIters;
info.maxiter = maxiter;
info.stoppedEarly = nIters < maxiter;
info.warmStart = warmStart;
info.nRestarts = nRestarts;
info.includeConstraint = includeConstraint;
info.W_all = W_all;
info.H_all = H_all;
info.M_all = M_all;
info.R_all = R_all;
info.cost_all = cost_all;
info.errors_all = errors_all;
info.loadings_all = loadings_all;
info.power_all = power_all;
info.seqNMF_W = seqNMF_W;
info.seqNMF_H = seqNMF_H;

end

%% ------------------------------------------------------------------------
function [ms, flexArgs] = parse_multistart_params(inputs)
    ms = struct();
    ms.nRestarts = 10;
    ms.warmStart = false;
    ms.seqNMF_maxiter = 50;
    ms.seqNMF_lambda = [];
    ms.includeConstraint = false;
    ms.verbal = 1;

    flexArgs = {};
    i = 1;
    while i <= numel(inputs)
        name = inputs{i};
        if i == numel(inputs)
            error('FlexMF_multistart:MissingValue', ...
                'Name/value pair for ''%s'' is missing a value.', name);
        end
        val = inputs{i+1};
        switch lower(char(name))
            case 'nrestarts'
                ms.nRestarts = val;
            case 'warmstart'
                ms.warmStart = logical(val);
            case 'seqnmf_maxiter'
                ms.seqNMF_maxiter = val;
            case 'seqnmf_lambda'
                ms.seqNMF_lambda = val;
            case 'includeconstraint'
                ms.includeConstraint = logical(val);
            case 'verbal'
                % Multistart verbal; also forward a quiet default to FlexMF separately
                ms.verbal = val;
            otherwise
                flexArgs(end+1:end+2) = {name, val};
        end
        i = i + 2;
    end

    if ~(isscalar(ms.nRestarts) && ms.nRestarts >= 1 && ms.nRestarts == floor(ms.nRestarts))
        error('FlexMF_multistart:BadNRestarts', 'nRestarts must be a positive integer.');
    end
end

function [K, L, lambda, lambda_M, lambda_R, EMD, maxiter] = get_flex_defaults(flexArgs)
    K = 10; L = 100; lambda = 1e-3; lambda_M = 1e-4; lambda_R = 1; EMD = 0;
    maxiter = 100;
    for i = 1:2:numel(flexArgs)
        switch lower(char(flexArgs{i}))
            case 'k', K = flexArgs{i+1};
            case 'l', L = flexArgs{i+1};
            case 'lambda', lambda = flexArgs{i+1};
            case 'lambda_m', lambda_M = flexArgs{i+1};
            case 'lambda_r', lambda_R = flexArgs{i+1};
            case 'emd', EMD = flexArgs{i+1};
            case 'maxiter', maxiter = flexArgs{i+1};
        end
    end
end

function args = set_default_nv(args, name, value)
    for i = 1:2:numel(args)
        if strcmpi(char(args{i}), name)
            return
        end
    end
    args = [args, {name, value}];
end

function args = drop_nv(args, names)
    keep = true(1, numel(args));
    for i = 1:2:numel(args)
        if any(strcmpi(char(args{i}), names))
            keep(i:i+1) = false;
        end
    end
    args = args(keep);
end
