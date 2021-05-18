function show_corr(ds1, idx1, ds2, idx2, corr_val, varargin)
% A multi-subplot visualization of spatiotemporal correlations between two
% cells, which is used by multiple "applications" (e.g. "match_1p2p",
% "resolve_depths").

% Default settings
display_mode = 'standard';

trace_norm_method = 'norm';
ds_labels = {'ds1', 'ds2'};
frames = [];
tform = [];

% In "overlay" display mode, we can center the cell map at the COM of the
% ds1 cell or the ds2 cell. By default, we target ds2.
zoom_target = 2;

color1 = [0 0.4470 0.7410];
color2 = [0.85 0.325 0.098];

% We use a custom "subplot" command that leaves less unusued space between
% panels
sp = @(m,n,p) subtightplot(m, n, p, 0.05, 0.05, 0.05); % Gap, Margin-X, Margin-Y

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'overlay'
                % Overlay the cellmaps of ds1 and ds2 after applying the
                % spatial transform 'tform' to the ds2 cellmap
                display_mode = 'overlay';
                tform = varargin{k+1};
            case 'zoom_target'
                zoom_target = varargin{k+1};
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
switch trace_norm_method
    case 'norm'
        y_lims = [0 1];
    case 'zsc'
        min_tr_val = min([min(tr1) min(tr2)]);
        max_tr_val = max([max(tr1) max(tr2)]);
        y_lims = [min_tr_val max_tr_val];
        y_lims = y_lims + 1/10*diff(y_lims)*[-1 1];
end


sp(3,1,1);
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
legend(ds_labels, 'Location', 'NorthWest', 'Interpreter', 'None');
xlim([1 length(tr1)]);
xlabel('Frames');
ylim(y_lims);
set(gca, 'TickLength', [0 0]);
title(sprintf('%s cell=%d; %s cell=%d; corr=%.4f',...
      ds_labels{1}, idx1, ds_labels{2}, idx2, corr_val), 'Interpreter', 'none');

% Show cell maps
%------------------------------------------------------------
switch display_mode
    case 'standard'
        sp(3,3,[4 7]); % Cellmap for ds1
        plot_boundaries(ds1, 'Color', color1, 'LineWidth', 1, 'Fill', idx1);
        title(ds_labels{1}, 'Interpreter', 'none');
        
        sp(3,3,[5 8]); % Cellmap for ds2
        plot_boundaries(ds2, 'Color', color2, 'LineWidth', 1, 'Fill', idx2);
        title(ds_labels{2}, 'Interpreter', 'none');
        
    	corr_sp = sp(3,3,[6 9]); % Correlation plot
        
    case 'overlay'
        sp(3,2,[3 5]); % Cellmap overlay
        
        switch zoom_target
            case 1
                zoom_com = ds1.cells(idx1).com;
            case 2
                if ~isempty(tform)
                    zoom_com = transformPointsForward(tform, ds2.cells(idx2).com')';
                else
                    zoom_com = ds2.cells(idx2).com;
                end
        end
        
        % For performance reasons, only show boundaries in the vicinity of
        % the cell under consideration
        plot_boundaries(ds1, 'Color', color1, 'LineWidth', 2, 'Fill', idx1, 'display_center', zoom_com);
        hold on;
        plot_boundaries(ds2, 'Color', color2, 'LineWidth', 1, 'Fill', idx2, 'Transform', tform, 'display_center', zoom_com);
        
        xlim(zoom_com(1) + [-50 50]);
        ylim(zoom_com(2) + [-50 50]);
        
        corr_sp = sp(3,2,[4 6]); % Correlation plot
end

% Show trace correlation
%------------------------------------------------------------
switch trace_norm_method
    case 'norm'
        plot(tr2, tr1, '.k');
        
        ticks = 0:0.1:1;
        set(corr_sp, 'XTick', ticks);
        set(corr_sp, 'YTick', ticks);
        
    case 'zsc'
        plot(tr2, tr1, '.k');
        
        set(corr_sp, 'XTick', 0:5:max(tr2));
        set(corr_sp, 'YTick', 0:5:max(tr1));
        
        hold off;
end

grid on;
axis equal tight;
xlabel(sprintf('%s (%s)', ds_labels{2}, trace_norm_method), 'Interpreter', 'none');
ylabel(sprintf('%s (%s)', ds_labels{1}, trace_norm_method), 'Interpreter', 'none');

end % show_corr