function plot_events_psth(ds, varargin)
% Plots the PSTH (peristimulus time histogram) for all cells (default).
%
% The PSTH will be computed over the frames that are common to all trials.

% Default parameters
cell_indices = find(ds.is_cell);
bin_width = 1; % Number of frames per bin
align_to = 3;

for j = 1:length(varargin)
    vararg = varargin{j};
    if ischar(vararg)
        switch lower(vararg)
            case {'cell', 'cells'} % Select cells to count
                cell_indices = varargin{j+1};

            case {'bin', 'bin_width'} % Number of bins to use over trial
                bin_width = varargin{j+1};
                
            case {'align', 'align_to'}
                align_to = varargin{j+1};
        end
    end
end

% Determine the number of frames that is common to all trials
[all_frames, common_frames] = compute_frame_bounds(ds.trial_indices, align_to);
num_frames = diff(all_frames) + 1;
num_bins = ceil(num_frames / bin_width);

f2b = @(f) frame2bin(f, all_frames, bin_width, num_bins);
bin0 = f2b(0); % Bin corresponding to the alignment frame

% Tally PSTH
x = (1:num_bins) - bin0;
y = zeros(size(x));
num_cells = length(cell_indices);
for k = 1:num_cells
    cell_idx = cell_indices(k);
    
    for trial_idx = 1:ds.num_trials
        eventdata = ds.get_events(cell_idx, trial_idx, 'align_to', align_to);
        if ~isempty(eventdata)
            event_times = eventdata(:,2); % Note: Using peak times!
            bin_locations = f2b(event_times);
            y(bin_locations) = y(bin_locations) + 1;
        end
    end
end
y_range = [0 max(y)];

% Display results
plot(x,y,'.-');
xlim(x([1 end]));
ylim(y_range);
xlabel(sprintf('Bin (each bin is %d frames; 0 contains the alignment frame)', bin_width));
ylabel(sprintf('Event tallies (over %d cells)', num_cells));
grid on;
hold on;
common_bin1 = f2b(common_frames(1)) - bin0;
common_bin2 = f2b(common_frames(2)) - bin0;
plot(common_bin1*[1 1], y_range, 'r--');
plot(common_bin2*[1 1], y_range, 'r--');
hold off;
end % plot_events_psth

function [all_frames, common_frames] = compute_frame_bounds(frame_indices, align_to)
	start_frames = frame_indices(:,1) - frame_indices(:,align_to);
    end_frames = frame_indices(:,4) - frame_indices(:,align_to);
    
    common_frames = [max(start_frames) min(end_frames)];
    all_frames = [min(start_frames) max(end_frames)];
end % compute_common_frames

function bin = frame2bin(frame, common_frames, frames_per_bin, num_bins)
    bin = floor((frame - common_frames(1))/frames_per_bin) + 1;
    bin = max(1, bin);
    bin = min(bin, num_bins);
end % frame2bin