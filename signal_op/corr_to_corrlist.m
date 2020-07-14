function corrlist = corr_to_corrlist(corr, varargin)
% Convert correlation matrix 'corr' into a 'corrlist', where each row k:
%   corrlist(k,:) = [i, j, corr(i,j)]

upper_triangular = false;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'upper' % Return only upper triangular portion of the correlation matrix
                upper_triangular = true;
        end
    end
end

% TODO: Upper diagonal option
[n_x, n_y] = size(corr);

X = repmat((1:n_x)', 1, n_y);
Y = repmat(1:n_y, n_x, 1);

corrlist = [X(:) Y(:) corr(:)];

if upper_triangular
    upper_mask = corrlist(:,2) > corrlist(:,1);
    corrlist = corrlist(upper_mask, :);
end