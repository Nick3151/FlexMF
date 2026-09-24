function z0 = dual_warm_start(dual, names, scales, affineF)
% Build the TFOCS SCD dual initial point from stored unscaled duals.
% dual:    struct of unscaled duals keyed by block name (may be empty/missing fields)
% names:   block names in the order of affineF rows
% scales:  proxScale of each block (block operator is scaled by 1/scale)
% affineF: TFOCS affine cell array, used for block output sizes
% Returns [] when nothing is stored, so TFOCS falls back to zeros.

nBlocks = numel(names);
if isempty(dual) || ~any(isfield(dual, names))
    z0 = [];
    return
end

z0 = cell(1, nBlocks);
for i = 1:nBlocks
    sz = affineF{i,1}([], 0);
    sz = sz{2};
    if isfield(dual, names{i}) && isequal(size(dual.(names{i})), sz)
        z0{i} = dual.(names{i}) * scales(i);
    else
        z0{i} = zeros(sz);
    end
end

if nBlocks == 1
    z0 = z0{1};
end
end
