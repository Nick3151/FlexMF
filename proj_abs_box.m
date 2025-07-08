function op = proj_abs_box( W )

%PROJ_LINF   Projection onto the scaled box around center.
%    OP = PROJ_LINF( W ) returns an operator implementing the 
%    indicator function for { X | abs(Xi) <= Wi, any i }. 
% Dual: prox_l1(Wx)

if any(~isnumeric( W(:) )) || any(~isreal( W(:) )) || any(~W(:))
	error( 'Argument must be positive.' );
end
op = @(varargin)proj_abs_box_W( W, varargin{:} );

function [ v, x ] = proj_abs_box_W( W, x, t )
assert(isequal(size(W),size(x)), 'Weights should have the same size as input.')
v = 0;
switch nargin
	case 2
		if nargout == 2
			error( 'This function is not differentiable.' );
		elseif any( abs(x(:))>W(:)) 
			v = Inf;
		end
	case 3	
        x = x ./ max( 1, abs( x ./ W ) );
    otherwise
		error( 'Not enough arguments.' );
end
