function [match_1to2, match_2to1, affine_info] = run_alignment(ds1, ds2, varargin)
% Align two sets of cell filters.
%
% Inputs:
%   ds1/2: DaySummary object containing cell maps to be aligned
%
% Optional inputs:
%   'notrans': No transformation applied between ds1 and ds2
%   'keepall': A cell from ds1 may match to more than one cell in ds2, and
%              visa versa.
%   Additional optional arugments are passed into 'match_masks.m'
%
% Outputs:
%   match_XtoY: Cell that contains mapping information from source X to
%       source Y.
%
%       For example, match_1to2{k} is a [Nx2] matrix where N is the number
%       of ICs from source 2 that match cell k of source 1. The first column
%       of the matrix is the index of the matching cell in source 2; the
%       second column is the overlap score between the cell pairs.
%
%   affine_info: Struct with additional information (e.g. selected points,
%       etc.) regarding the affine transformation applied to match the
%       datasets.
%
% Example usage:
%   [m1to2, m2to1] = run_alignment('c9m7d07_ica001', 'c9m7d08_ica001');
%

% Default alignment options
%------------------------------------------------------------
use_transform = 1;
bijective_matching = 1;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch varargin{k}
            case 'notrans'
                % No affine transform needed (e.g. for matching multiple
                % extraction runs on the same movie). Requires the image
                % dimensions to be identical!
                assert(all(size(ds1.cells(1).mask)==size(ds2.cells(1).mask)),...
                    'notrans option requires cell image dimensions to be identical!');
                use_transform = 0;
            case 'keepall' % Keep all matches
                bijective_matching = 0;                
        end
    end
end

% Control point-based registration of two sets of ICs
%------------------------------------------------------------
if use_transform
    fprintf('run_alignment: Beginning alignment...\n');
    [affine_info, masks1, masks2] = compute_affine_transform(ds1, ds2);
else
    masks1 = {ds1.cells.mask};
    masks2 = {ds2.cells.mask};
    affine_info = [];
    
    figure;
    plot_boundaries_with_transform(ds1, 'b', 2, [], []);
    plot_boundaries_with_transform(ds2, 'r', 1, [], []);
    title('Dataset1 (blue) vs. Dataset2 (red)');
end
input('run_alignment: Press any key to continue with mask matching >> ');

[match_1to2, match_2to1] = match_masks(masks1, masks2, ds1, ds2, varargin{:});

% Optional bijective filtering
%------------------------------------------------------------
if bijective_matching
    fprintf('run_alignment: Applying bijective filter...\n');
    [match_1to2, match_2to1] = bijective_filter(match_1to2, match_2to1);
end