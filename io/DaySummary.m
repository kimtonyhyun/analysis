classdef DaySummary
% Summary of PlusMaze data for a single day.
%
% Inputs:
%   plusmaze_txt: Name of the PlusMaze output text file. Frame indices of
%       the text file needs to be consistent with the ICA traces!
%   ica_dir: Directory containing 
%       - ICA results in a "ica_*.mat" file
%       - Classification results in a "class_*.txt" file
%
% Output:
%   DaySummary object
%
% Example usage:
%   ds = DaySummary('c9m7d06_ti2.txt', 'ica001');
%
    properties
        cells
        trials
        
        num_cells
        num_trials
        
        trial_indices
    end
    
    methods
        function obj = DaySummary(plusmaze_txt, ica_dir, varargin)
            % Handle optional input
            exclude_probe_trials = 0;
            for k = 1:length(varargin)
                switch lower(varargin{k})
                    case 'excludeprobe'
                        exclude_probe_trials = 1;
                end
            end
            
            % Load data
            %------------------------------------------------------------
            [trial_indices, loc_info, trial_durations] =...
                parse_plusmaze(plusmaze_txt);
            
            data = load(get_most_recent_file(ica_dir, 'ica_*.mat'));
            
            % Parse trial data
            %   TODO: Bring in centroids corresponding to mouse position
            %------------------------------------------------------------
            if (exclude_probe_trials)
                is_probe = strcmp(loc_info(:,1), 'north') | ...
                           strcmp(loc_info(:,1), 'south');
                       
                trial_indices = trial_indices(~is_probe,:);
                loc_info = loc_info(~is_probe,:);
                trial_durations = trial_durations(~is_probe);
            end
            
            num_trials = size(trial_indices, 1);
            traces = cell(num_trials, 1);
            for k = 1:num_trials
                trial_frames = trial_indices(k,1):...
                               trial_indices(k,end);
                traces{k} = data.ica_traces(trial_frames, :)';
            end
            
            obj.num_trials = num_trials;
            obj.trial_indices = trial_indices;
            obj.trials = struct(...
                'start', loc_info(:,1),...
                'goal',  loc_info(:,2),...
                'end',   loc_info(:,3),...
                'time',  num2cell(trial_durations),...
                'traces', traces);
            
            % Parse cell data
            %------------------------------------------------------------
            obj.num_cells = data.ica_info.num_ICs;
            
            class = load_classification(get_most_recent_file(ica_dir, 'class_*.txt'));
            obj.cells = struct(...
                'im', squeeze(num2cell(data.ica_filters, [1 2])),...
                'label', class);
        end
        
        function [trace, frame_indices] = get_cell_trial(obj, cell_idx, trial_idx)
            trace = obj.trials(trial_idx).traces(cell_idx,:);
            frame_indices = obj.trial_indices(trial_idx,1):...
                            obj.trial_indices(trial_idx,end);
        end
    end
end