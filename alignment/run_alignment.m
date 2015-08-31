function [match_1to2, match_2to1, info] = run_alignment(ds1, ds2, varargin)
% Align two sets of IC filters.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%
% Optional inputs:
%   'one-to-one': Each cell of Dataset1 will match at most one cell of
%       Dataset 2, and visa versa.
%
% Outputs:
%   match_XtoY: Cell that contains mapping information from source X to
%       source Y.
%
%       For example, match_1to2{k} is a [Nx2] matrix where N is the number
%       of ICs from source 2 that match IC k of source 1. The first column
%       of the matrix is the index of the matching IC in source 2; the
%       second column is the overlap score between the IC pairs.
%
% Example usage:
%   [m1to2, m2to1] = run_alignment('c9m7d07_ica001', 'c9m7d08_ica001');
%

bijective_matching = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case {'one-to-one', 'bijective'}
                bijective_matching = 1;
        end
    end
end

% Control point-based registration of two sets of ICs
%------------------------------------------------------------
fprintf('run_alignment: Beginning alignment...\n');    
[affine_info, masks1, masks2_tform] = compute_affine_transform(ds1, ds2);

% Compute matrix of mask overlaps
%   Note that the overlap matrix is non-symmetric!
%------------------------------------------------------------
input('run_alignment: Press any key to continue with mask matching >> ');
M = zeros(affine_info.num_cells1, affine_info.num_cells2);
for i = 1:affine_info.num_cells1
    if (mod(i,20)==0)
        fprintf('  %s: Computing overlaps (%.1f%%)...\n',...
            datestr(now), 100*i/affine_info.num_cells1);
    end
    for j = 1:affine_info.num_cells2
        M(i,j) = compute_mask_overlap(masks1{i}, masks2_tform{j});
    end
end
fprintf('  %s: Overlap matrix completed!\n', datestr(now));

% Find the nonzero elements of the mask overlap matrix
%------------------------------------------------------------
overlap_threshold = 1/3;
match_1to2 = cell(size(M,1),1);
for i = 1:size(M,1)
    m = M(i,:);
    matching_js = find(m>overlap_threshold);
    if isempty(matching_js)
        match_1to2{i} = [];
    else
        matching_overlaps = m(matching_js);
        match_1to2{i} = [matching_js' matching_overlaps']; % [j overlap]
        match_1to2{i} = sortrows(match_1to2{i},2); % Sort by overlap
        match_1to2{i} = flipud(match_1to2{i}); % Descending
    end
end

match_2to1 = cell(size(M,2),1); % Same in the opposite direction
for j = 1:size(M,2)
    m = M(:,j)';
    matching_is = find(m>overlap_threshold);
    if isempty(matching_is)
        match_2to1{j} = [];
    else
        matching_overlaps = m(matching_is);
        match_2to1{j} = [matching_is' matching_overlaps'];
        match_2to1{j} = sortrows(match_2to1{j},2);
        match_2to1{j} = flipud(match_2to1{j});
    end
end

% Set up auxiliary output
%------------------------------------------------------------
info.affine = affine_info;
info.masks1 = masks1;
info.masks2 = masks2_tform;
info.overlap_threshold = overlap_threshold;
info.overlap_matrix = M;

% Optional bijective filtering
%------------------------------------------------------------
if bijective_matching
    fprintf('run_alignment: Applying bijective filter...\n');
    [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1);
end