function dual = dual_unscale(dual, z, names, scales)
% Store TFOCS SCD duals z (unscaled by each block's proxScale) into struct dual.
% Fields for blocks not in names are kept, so duals shared across H/W updates
% ('fit', 'M') are overwritten while block-specific ones are preserved.

if isempty(dual)
    dual = struct();
end
if isa(z, 'tfocs_tuple')
    z = cell(z);
elseif ~iscell(z)
    z = {z};
end

for i = 1:numel(names)
    dual.(names{i}) = z{i} / scales(i);
end
end
