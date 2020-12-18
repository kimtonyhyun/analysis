function trace = filter_trace(ds, cell_idx, fps, cutoff_freq)
% Compute a low-pass filtered version of the fluorescence trace. Used in
% event computation.

% Default parameters comes from Wagner et al. Cell 2019, where we used
%   - 4 Hz cutoff frequency
%   - 30 Hz sampling frequency
if ~exist('cutoff_freq', 'var')
    cutoff_freq = 4;
end

% Don't apply LPF across trial boundaries
filt_traces = cell(1,ds.num_trials);
for k = 1:ds.num_trials
    tr_k = ds.get_trace(cell_idx, k);
    filt_traces{k} = filter_segment(tr_k, cutoff_freq, fps);
end
trace = cell2mat(filt_traces);

end

function trf = filter_segment(tr, cutoff_freq, sampling_freq)
% Zero-phase filtering by 2nd-order Butterworth lowpass filter

[b,a] = butter(2, cutoff_freq/(sampling_freq/2));
trf = filtfilt(b,a,double(tr));

end
