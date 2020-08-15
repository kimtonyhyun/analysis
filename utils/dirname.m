function name = dirname(recursion_idx, path)
% Returns the name of the current working directory. Can retrieve the name
% of upstream directories by using 'recursion_idx'. For example, to return 
% the name of the direct parent, use: dirname(1)
%   

if ~exist('recursion_idx', 'var')
    recursion_idx = 0;
end
if ~exist('path', 'var')
    path = pwd;
end

[path, top] = fileparts(path);
if (recursion_idx == 0)
    name = top; 
elseif (recursion_idx > 0)
    name = dirname(recursion_idx-1, path);
else
    error('Recursion index must be a positive integer');
end