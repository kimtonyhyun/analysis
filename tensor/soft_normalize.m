function [Y] = soft_normalize(X)
% SOFT_NORMALIZE, apply soft-normalization to each neuron in data tensor X and
% return a new, normalized tensor Y.
%
%     [Y] = SOFT_NORMALIZE(X)

Y = zeros(size(X));
if ndims(X) == 3
	for c = 1:size(X,1)
		x = X(c,:,:);
		Y(c,:,:) = x ./ (1 + max(x(:)));
	end
else
	error('input has inappropriate dimensions')
end
