function [rec_savename, class_savename] = generate_cascade_rec(path_to_dff)
% Generate a "rec" file using CASCADE's spike rates as "traces"

input_rec_filename = get_most_recent_file(path_to_dff, 'rec_*.mat');
cascade_filename = get_most_recent_file(path_to_dff, 'cascade_*.mat');

rec = load(input_rec_filename);
cascade = load(cascade_filename);

info.type = 'cascade';
info.num_pairs = rec.info.num_pairs;

info.cascade.input_file = cascade.rec_file;
info.cascade.model_name = cascade.model_name;
info.cascade.dff_traces = rec.traces;

[rec_savename, timestamp] = save_rec(info, rec.filters, cascade.spike_probs);
class_savename = generate_class_file(info.num_pairs, 'timestamp', timestamp);