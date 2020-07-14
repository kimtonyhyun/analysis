function corrlist = corr_to_corrlist(corr)
% Convert correlation matrix 'corr' into a 'corrlist', where each row k:
%   corrlist(k,:) = [i, j, corr(i,j)]

% TODO: Upper diagonal option
[n_x, n_y] = size(corr);

X = repmat((1:n_x)', 1, n_y);
Y = repmat(1:n_y, n_x, 1);

corrlist = [X(:) Y(:) corr(:)];