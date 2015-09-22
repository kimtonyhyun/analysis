function [U,S,V] = svds_custom(M, k)
% Functionally same as svds(), but does eigs() on covariance matrix instead.

[num_pixels,num_frames] = size(M);

fprintf('%s: Computing covariance matrix...\n',datestr(now));

if num_pixels>num_frames
    C = double(M'*M);
    fprintf('%s: Computing right (temporal) singular vectors...\n', datestr(now));
    options.issym = 'true';
    [V, S2] = eigs(C, k, 'LM', options);
    S = diag(diag(S2).^(1/2));
    fprintf('%s: Computing corresponding left (spatial) singular vectors...\n', datestr(now));
    U = (M * V) / S;    
    
else
    C = double(M*M');
    fprintf('%s: Computing left (spatial) singular vectors...\n', datestr(now));
    options.issym = 'true';
    [U, S2] = eigs(C, k, 'LM', options);
    S = diag(diag(S2).^(1/2));
    fprintf('%s: Computing corresponding right (temporal) singular vectors...\n', datestr(now));
    VT = S \ (U'*M);
    V = VT';
end

U =single(U);
V = single(V);

end