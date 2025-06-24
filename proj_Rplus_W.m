function op = proj_Rplus_W(KL)
% Projection of W_flat to the nonnegative orthant
% x = [W_flat M R]
% W_flat: N*KL

narginchk(1,1)
op = @(varargin)proj_Rplus_W_impl(KL,varargin{:});

function [ v, x ] = proj_Rplus_W_impl( KL, x, t )
v = 0;
switch nargin
	case 2
		if nargout == 2
			error( 'This function is not differentiable.' );
        elseif any( x(:,1:KL) < 0, 'all' )
            v = Inf;
        end
	case 3
        W_flat = x(:,1:KL);
		x = [max( W_flat, 0 ) x(:,KL+1:end)];
    otherwise
        error( 'Wrong number of arguments' );
end
		