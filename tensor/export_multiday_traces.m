function [matched_traces, meta, neuron_map, trial_map] = export_multiday_traces(md, varargin)
% Exports cross-day aligned cell traces (via MultiDay) into a lightweight
% format (description below) for tensor analysis
%
% Output format:
%   matched_traces: Column cell vector [M x 1] where M is the total number 
%      of trials across all days of the requested 'trial_type' (e.g. 'en')
%       
%      matched_traces{j} is a matrix of cell traces, in the format
%      [num_matched_cells x num_frames_in_trial_j].
%
%   neuron_map: Matrix [num_matched_cells x num_days] that maps each 
%      matched neuron to its per-day neuron index.
%
%   trial_map: Matrix [M x 2] where M is the total number of trials of 
%      requested 'trial_type'. For the i-th trial, K(i,1) is the day index 
%      of that trial, and K(i,2) is the trial index in the original day.
%   
max_num_trials = 0;
for day_idx = md.valid_days
    max_num_trials = max_num_trials + md.day(day_idx).num_trials;
end

neuron_map = md.matched_indices;

% get trial map (filtering out those specified)
trial_map = filter_trials(md, varargin{:});
K = size(trial_map,1);

% copy selected trials into lightweight cell array
matched_traces = cell(K,1);
meta.start = cell(K,1);
meta.end = cell(K,1);
meta.correct = cell(K,1);
meta.day = cell(K,1);
for k = 1:K
    d = trial_map(k,1); % day for this trial
    ni = neuron_map(:,md.valid_days == d);
    trial = md.day(d).trials(trial_map(k,2));
    matched_traces{k} = trial.traces(ni,:);
    meta.start{k} = trial.start;
    meta.end{k} = trial.end;
    meta.correct{k} = num2str(trial.correct);
    meta.day{k} = num2str(d);
end
