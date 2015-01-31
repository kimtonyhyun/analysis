function frame_indices = get_trial_frame_indices(source)
% Compatible with the plus maze output as of Cohort 9

fid = fopen(source);
maze_data = textscan(fid, '%s %s %s %f %d %d %d %d');
fclose(fid);

frame_indices = [maze_data{5} maze_data{6} maze_data{7} maze_data{8}];