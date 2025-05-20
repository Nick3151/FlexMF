function X = generate_sequence(T,N,Dt,varargin)

p  = inputParser;
addOptional(p, 'noise', 0, @isscalar)
addOptional(p, 'shift', 0, @isscalar)
addOptional(p, 'warp', 0, @isscalar)
addOptional(p, 'len_burst', 1, @isscalar)
addOptional(p, 'dynamic', 0, @isscalar)
parse(p, varargin{:})

noise = p.Results.noise;
shift = p.Results.shift;
warp = p.Results.warp;
len_burst = p.Results.len_burst;
dynamic = p.Results.dynamic;

lseq = (Dt+warp).*N+len_burst;
assert(lseq<T, 'Not enough total length T!')
X_hat = zeros(N,T);

% Apply warping and shift to sequence
for n=1:N
    X_hat(n,(Dt+warp)*(n-1)+(1:len_burst)+shift) = 1;
end

% Add noise
X_noise = (rand(size(X_hat))<noise);
X = X_hat + (~X_hat).*X_noise;

% Apply calcium dynamics
if dynamic
    filtbio = exp(-(0:8)/4);
    X = conv2(X, filtbio);
end
