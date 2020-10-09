function show_corr(ds1, idx1, ds2, idx2, corr_val, varargin)

% Default settings
display_mode = 'standard';

trace_norm_method = 'norm';
ds_labels = {'ds1', 'ds2'};
frames = [];
tform = [];

color1 = [0 0.4470 0.7410];
color2 = [0.85 0.325 0.098];

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'overlay'
                % Overlay the cellmaps of ds1 and ds2 after applying the
                % spatial transform 'tform' to the ds2 cellmap
                display_mode = 'overlay';
                tform = varargin{k+1};
            case 'zsc'
                trace_norm_method = 'zsc';
            case {'name', 'names', 'ds_name', 'ds_names'}
                if iscell(varargin{k+1})
                    ds_labels = varargin{k+1};
                elseif ischar(varargin{k+1})
                    ds_labels = varargin(k+1);
                end
            case 'frames' % Indicate frames with vertical bar
                frames = varargin{k+1};
        end
    end
end

% Display data
tr1 = ds1.get_trace(idx1, trace_norm_method)'; % Column vector
tr2 = ds2.get_trace(idx2, trace_norm_method)';

if isempty(corr_val)
    corr_val = corr(tr1, tr2);
end

% Show traces
%------------------------------------------------------------
subplot(311);
plot(tr1, 'Color', color1);
hold on;
plot(tr2, 'Color', color2);
for m = 2:ds1.num_trials % Trial boundaries
    xline(ds1.trial_indices(m,1), 'k:');
end
for m = 1:length(frames) % Extra vertical markers
    xline(frames(m), 'b:');
end
hold off;
legend(ds_labels, 'Location', 'NorthWest');
xlim([1 length(tr1)]);
xlabel('Frames');
set(gca, 'TickLength', [0 0]);

title(sprintf('%s cell=%d\n%s cell=%d\ncorr=%.4f',...
      ds_labels{1}, idx1, ds_labels{2}, idx2, corr_val));

switch display_mode
    case 'standard'
        subplot(3,3,[4 7]); % Cellmap for ds1
        plot_boundaries_with_transform(ds1, color1, 1, idx1, []);
        title(ds_labels{1});
        
        subplot(3,3,[5 8]); % Cellmap for ds2
        plot_boundaries_with_transform(ds2, color2, 1, idx2, []);
        title(ds_labels{2});
        
    	corr_sp = subplot(3,3,[6 9]); % Correlation plot
        
    case 'overlay'
        subplot(3,2,[3 5]); % Cellmap overlay
        plot_boundaries_with_transform(ds1, color1, 2, idx1, []);
        hold on;
        plot_boundaries_with_transform(ds2, color2, 1, idx2, tform);
        
        corr_sp = subplot(3,2,[4 6]); % Correlation plot
end

plot(tr2, tr1, '.k');

switch trace_norm_method
    case 'norm'
        ticks = 0:0.1:1;
    case 'zsc'
        ticks = -50:5:100; % FIXME: Hard-coded   
end
set(corr_sp, 'XTick', ticks);
set(corr_sp, 'YTick', ticks);
grid on;
axis equal tight;
xlabel(sprintf('%s (%s)', ds_labels{2}, trace_norm_method));
ylabel(sprintf('%s (%s)', ds_labels{1}, trace_norm_method));