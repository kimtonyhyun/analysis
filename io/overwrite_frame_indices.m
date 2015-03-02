function overwrite_frame_indices(plusmaze_source, new_frame_indices)

% Read in information from the original PlusMaze source
[~, location_info, time] = parse_plusmaze(plusmaze_source);
num_trials = length(time);

% Generate a new PlusMaze output
[path_to, name, ~] = fileparts(plusmaze_source);
output_name = sprintf('%s_new.txt', name);
output_name = fullfile(path_to, output_name);

fid = fopen(output_name, 'w');
for i = 1:num_trials
    fprintf(fid, '%s %s %s %.3f %d %d %d %d\n',...
        location_info{i,1}, location_info{i,2}, location_info{i,3},...
        time(i),...
        new_frame_indices(i,1), new_frame_indices(i,2),...
        new_frame_indices(i,3), new_frame_indices(i,4));
end
fclose(fid);