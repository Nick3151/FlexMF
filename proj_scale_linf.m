function op = proj_scale_linf( Q )

%PROJ_SCALE_LINF   Projection onto the scaled infinity norm ball.
%    OP = PROJ_LINF( Q ) returns an operator implementing the 
%    indicator function for the infinity norm ball of size Q,
%    { X | |X| <= Q }. Q is optional; if omitted,
%    Q=1 is assumed. But if Q is supplied, each entry must be non-negative

% Dual: sum(Q.*|X|)
% See also: prox_l1, prox_linf, proj_l1

if nargin == 0,
	Q = 1;
elseif ~isnumeric( Q ) || ~isreal( Q ) || any(Q(:) < 0)
	error( 'Argument must be non-negative.' );
end
op = @(varargin)proj_scale_linf_Q( Q, varargin{:} );

function [ v, X ] = proj_scale_linf_Q( Q, X, t )
v = 0;
switch nargin,
	case 2,
		if nargout == 2,
			error( 'This function is not differentiable.' );
		elseif any(abs(X(:))-Q(:)>0),
			v = Inf;
		end
	case 3,			
        X = min(abs(X),Q).*sign(X);
	otherwise,
		error( 'Not enough arguments.' );
end

% TFOCS v1.3 by Stephen Becker, Emmanuel Candes, and Michael Grant.
% Copyright 2013 California Institute of Technology and CVX Research.
% See the file LICENSE for full license information.