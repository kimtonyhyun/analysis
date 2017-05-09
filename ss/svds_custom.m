function [U,S,V] = svds_custom(M, k,verbose)
% Functionally same as svds(), but does eigs() on covariance matrix instead.

if ~exist('verbose','var'), verbose = 0; end

[num_pixels,num_frames] = size(M);

dispfun(sprintf('%s: Computing covariance matrix...\n',datestr(now)),verbose);

if num_pixels>num_frames
    C = double(M'*M);
    dispfun(sprintf('%s: Computing right (temporal) singular vectors...\n', datestr(now)),verbose);
    options.issym = 'true';
    [V, S2] = eigs(C, k, 'LM', options);
    S = diag(diag(S2).^(1/2));
    dispfun(sprintf('%s: Computing corresponding left (spatial) singular vectors...\n', datestr(now)),verbose);
    U = (M * V) / S;    
    
else
    C = double(M*M');
    dispfun(sprintf('%s: Computing left (spatial) singular vectors...\n', datestr(now)),verbose);
    options.issym = 'true';
    [U, S2] = eigs(C, k, 'LM', options);
    S = diag(diag(S2).^(1/2));
    dispfun(sprintf('%s: Computing corresponding right (temporal) singular vectors...\n', datestr(now)),verbose);
    VT = S \ (U'*M);
    V = VT';
end

U =single(U);
V = single(V);

function dispfun(str,state)
    if state == 1,
        fprintf(str);
    end
end

end