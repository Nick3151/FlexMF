function op = proj_Rplus_H(K)
% Projection of H to the nonnegative orthant
% x = [H; M; R]
% H: K*T

narginchk(1,1)
op = @(varargin)proj_Rplus_H_impl(K,varargin{:});

function [ v, x ] = proj_Rplus_H_impl( K, x, t )
v = 0;
switch nargin
	case 2
		if nargout == 2
			error( 'This function is not differentiable.' );
        elseif any( x(1:K,:) < 0, 'all')
            v = Inf;
        end
	case 3
        H = x(1:K,:);
		x = [max( H, 0 ); x(K+1:end,:)];
    otherwise
		error( 'Wrong number of arguments' );
end