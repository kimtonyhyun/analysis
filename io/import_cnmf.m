function import_cnmf(cnmf_struct, varargin)
% Converts the CNMF output into a format that can be directly read by
% `classify_cells`, and saves it into a `rec_*.mat` file.
%
% Example usage:
% output = 
% 
%                          success: 1
%                           params: [1x1 struct]
%     datasetComponentProperties_P: [1x1 struct]
%                        movieList: ''
%                  extractedImages: [512x512x349 double]
%                 extractedSignals: [349x300 double]
%                   extractedPeaks: [349x300 double]
%                      cnmfOptions: [1x1 struct]
%                             time: [1x1 struct]
% 
% import_cnmf(output);

info.type = 'cnmf';

filters = cnmf_struct.extractedImages;
traces = cnmf_struct.extractedSignals'; %#ok<*NASGU>

info.num_pairs = size(filters, 3); %#ok<*STRNU>

% Save to file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');
