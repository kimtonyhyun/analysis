function [ Xnrm ] = normalize_tensor( X, meta )
%NORMALIZE_TENSOR Regress away main effect of day
%   X = NORMALIZE_TENSOR(X)

% dimensions
[N,T,K] = size(X);
Xnrm = X;

% get indices for day
days = unique(meta.day);
ndays = length(days);
d = zeros(K,1);
tn = zeros(K,1); % trial within day
for d_ = 1:ndays
    idx = meta.day == days(d_);
    d(idx) = d_;
    tn(idx) = 1:sum(idx);
end

% subtract off average
for d_ = 1:ndays
    idx = (d==d_);
    Xd = X(:,:,idx); % data on day d_
    avg_d = mean(Xd(:,:),2); % average for each cell on day d_
    Xnrm(:,:,idx) = Xd - repmat(avg_d, 1, T, sum(idx));
end

