function [frame_indices,location_info,time] = getTrialIndicies(source)

fid = fopen(source);
maze_data = textscan(fid, '%s %s %s %f %d %d %d %d');
fclose(fid);

frame_indices = [maze_data{5} maze_data{6} maze_data{7} maze_data{8}];
location_info = [maze_data{1} maze_data{2} maze_data{3}];
time = maze_data{4};