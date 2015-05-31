function import_cellmax(cellmax_struct)
% Converts the CellMax output into a format that can be directly read by
% `classify_cells`, and saves it into a `rec_*.mat` file.
%
% Example usage:
%     m7d07_output = 
% 
%             dsMovieFilename: 'F:\c9m7d07\c9m7d07_cr_mc_cr_norm_dff_ti2.hdf5'
%               movieFilename: []
%                  cellImages: [499x500x507 double]
%                dsCellTraces: [507x25462 double]
%                   centroids: [507x2 double]
%         dsScaledProbability: [507x25462 double]
%                   EMoptions: [1x1 struct]
%                dsFiltTraces: [507x25462 double]
%                dsEventTimes: {507x1 cell}
%                eventOptions: [1x1 struct]
% 
%     import_cellmax(m7d07_output);

info.type = 'cellmax';
info.cellmax_source = inputname(1); % Workspace variable name of `cellmax_struct`

filters = cellmax_struct.cellImages; %#ok<*NASGU>
traces = cellmax_struct.dsCellTraces';

% Save the CellMax traces & filters to a mat file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

save(rec_savename, 'info', 'filters', 'traces');