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

% subtract off averages
for d_ = 1:ndays
    idx = (d==d_);
    Xs = X(:,:,idx);
    Xnrm(:,:,idx) = Xnrm(:,:,idx) - repmat(mean(Xs(:,:),2),1,T,sum(idx));
    % Xnrm(:,:,idx) = Xnrm(:,:,idx) ./ repmat(std(Xs(:,:),[],2),1,T,sum(idx));
end

% % subtract off averages
% for tn_ = 1:max(tn)
%     idx = (tn==tn_);
%     Xnrm(:,:,idx) = Xnrm(:,:,idx) - repmat(mean(X(:,:,idx),3),1,1,sum(idx));
% end

