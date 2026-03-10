function Xcorr = correct_warp(X, M)
T = size(X,2);
assert(size(M,2)==T, 'Size of data and M does not match!')

% Correct the temporal warpped/jittered part of X
D = eye(T) - diag(ones(T-1,1),-1);
D(T,T) = 0;
Xcorr = M*D'+X;

end