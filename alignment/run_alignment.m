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

input('run_alignment: Press any key to continue with mask matching >> ');
[match_1to2, match_2to1, M] = match_masks(masks1, masks2_tform);

% Set up auxiliary output
%------------------------------------------------------------
info.affine = affine_info;
info.masks1 = masks1;
info.masks2 = masks2_tform;
info.overlap_matrix = M;

% Optional bijective filtering
%------------------------------------------------------------
if bijective_matching
    fprintf('run_alignment: Applying bijective filter...\n');
    [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1);
end