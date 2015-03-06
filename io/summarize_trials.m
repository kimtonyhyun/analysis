function summarize_trials(sources,traces,varargin)
% Outputs a .mat file with [# of trials] x 4 cell having the following 4 columns:
%   trial start arm
%   trial finish arm
%   trial correctness
%   cell traces belonging to the trial
%
% inputs:
%   sources : Struct containing maze output, trim and binning parameters
%   needed here.
%   traces : [# of frames] x [#of ICs] array of cell traces.
%
% variable input arguments:
%   includeProbes: Set to 0 to retain all the trials, set to 1 to get rid
%   of probe trials.
%   classFile: Path to the classification file. Specifying this will
%   generate output for only the ICs that are classified as cells.
%
% Hakan Inan (Mar 15)
%
include_probes = 1; % Default
if ~isempty(varargin)
    len = length(varargin);
    for k = 1:len
        switch varargin{k}
            case 'includeProbes'
                include_probes = varargin{k+1};
                if include_probes~=0 && include_probes~=1
                    error('includeProbes can be either 0 or 1.');
                end
            case 'classFile'
                class_file = varargin{k+1};
                if ~ischar(class_file)
                    error('classFile must be a string containing path to the classification file.');
                end
                classes = load_classification(class_file);
                na = strcmp(classes,'not a cell');
                traces(:,na) = [];
        end
    end
end

% Get trial info from maze output
[trial_frame_indices,location_info,~] = parse_plusmaze(sources.maze);
trial_start_arms  = location_info(:,1);
trial_finish_arms = location_info(:,2);
trial_correctness = strcmp(location_info(:,2),location_info(:,3));
compressed_indices = compress_frame_indices(trial_frame_indices,sources.trim);
compressed_indices = bin_frame_indices(compressed_indices,sources.time_bin);
num_trials = size(trial_frame_indices,1);

summary_session = cell(num_trials,4);
for trial = 1:num_trials
    frames = compressed_indices(trial,[1,4]);
    summary_session{trial,1} = trial_start_arms{trial};
    summary_session{trial,2} = trial_finish_arms{trial};
    summary_session{trial,3} = trial_correctness(trial);
    summary_session{trial,4} = traces(frames(1):frames(2),:)';
end

if include_probes==0
    probe_north = find(strcmp(trial_start_arms,'north'));
    probe_south = find(strcmp(trial_start_arms,'south'));
    summary_session([probe_north,probe_south],:) = [];
end

[~, folder, ~] = fileparts(pwd);
save_name = strcat('summary_',folder);
save(save_name,'summary_session');