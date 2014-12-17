function [ percent_correct ] = calc_behavior_score( file )
%calc_behavior_score Reads behavior text file and gives behavior score as % correct
%not for autotraining results!
%2014-12-17 Fori Wang

    fid = fopen(file);
    maze_data = textscan(fid,'%s %s %s %f %d %d');
    fclose(fid);

    end_arm = maze_data{2};
    actual_end_arm = maze_data{3};
    num_trials = length(actual_end_arm);
    
    num_correct_instances = length(find(strcmp(end_arm,actual_end_arm)));
    percent_correct = num_correct_instances*100/num_trials;

end