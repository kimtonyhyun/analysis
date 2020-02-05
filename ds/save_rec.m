function rec_savename = save_rec(info, filters, traces) %#ok<*INUSD>
% Data format should be as follows:
%   - info.type: Name of the cell extraction method (e.g. 'ica', 'cellmax')
%   - info.num_pairs = Number of filter-trace pairs in the file
%   - filters: [height x width x num_pairs]
%   - traces: [num_frames x num_pairs]
%
% TODO: Sanity checks on rec file structure!

% Save to file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);
save(rec_savename, 'info', 'filters', 'traces', '-v7.3');
