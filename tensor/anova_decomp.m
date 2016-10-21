function [ Xdpca ] = anova_decomp( X, meta )
%ANOVA_DECOMP Analysis-of-variance decomposition of data tensor
%   Xphi = ANOVA_DECOMP(X, meta)

% data dimensions
[N,T,K] = size(X);

% get indices for starts
s = zeros(K,1);
s(strcmp(meta.start,'east')) = 1;
s(strcmp(meta.start,'west')) = 2;

% get indices for ends
e = zeros(K,1);
e(strcmp(meta.end,'north')) = 1;
e(strcmp(meta.end,'south')) = 2;

% get indices for day
days = unique(meta.day);
ndays = length(days);
d = zeros(K,1);
for d_ = 1:ndays
    d(meta.day == days(d_)) = d_;
end

% find averages
Xdpca = nan(N,T,2,2,ndays);
for s_ = 1:2
    for e_ = 1:2
        for d_ = 1:ndays
            idx = (s==s_) & (e==e_) & (d==d_);
            Xdpca(:,:,s_,e_,d_) = mean(X(:,:,idx),3);
        end
    end
end
