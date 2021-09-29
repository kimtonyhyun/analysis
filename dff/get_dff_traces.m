function [rec_savename, class_savename] = get_dff_traces(ds, fps, varargin)
% Compute DFF traces from the DaySummary 'ds'. The baseline is estimated by
% performing a polynomial fit to "nonactive" periods of the trace, which is
% interactively identified by user (via 'polyfit_nonactive_frames').
%
% The optional input 'F' is the full-field fluorescence of the movie. If
% provided, each DFF trace will be decorrelated with respect to 'F'. The
% user can then choose between the decorrelated and non-decorrelated
% traces. Note: In the future, we could perform multi-variate decorrelation
% with respect to other signals as well.
%
% TODO:
%   - A save/load mechanism?
%

color = [0 0.447 0.741];
color2 = [0.85 0.325 0.098];

F = []; % Full-field fluorescence for decorrelation 

compute_zsc = false;

for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'f'
                fprintf('Full-field fluorescence provided\n');
                F = varargin{k+1};
                F = F - mean(F);
                F = F(:)'; % ds.get_trace returns traces as a row vector. Need to match
            case 'zsc'
                fprintf('Output traces will be z-scored\n');
                compute_zsc = true;
        end
    end
end

% Custom "subplot" command that leaves less unusued space between panels
sp = @(m,n,p) subtightplot(m, n, p, 0.025, [0.1 0.05], [0.05 0.01]); % Gap, Margin-X, Margin-Y

gui.fig_handle = figure;
gui.h_orig = sp(2,1,1); % Original trace
gui.h_dff = sp(2,1,2);  % DFF trace
linkaxes([gui.h_orig, gui.h_dff], 'x');
set(gui.h_orig, 'TickLength', [0 0]);
set(gui.h_orig, 'XTickLabels', []);
set(gui.h_dff, 'TickLength', [0 0]);
grid(gui.h_orig, 'on');
grid(gui.h_dff, 'on');

default_params = struct('threshold', [],...
                        'padding', floor(4*fps),...
                        'order', 1);
params = default_params;

num_cells = ds.num_cells;
num_frames = ds.full_num_frames;
t = 1/fps*(0:num_frames-1);
t_lims = t([1 end]);

% We will populate the following variables
dff_traces = zeros(num_frames, num_cells, 'single');
baseline_fit_infos = cell(num_cells, 1);

% Flag where:
%   - keep_dff{k} == true;  % Keep the k-th cell in export
%   - keep_dff{k} == false; % Omit the k-th cell in export
%   - keep_dff{k} == [];    % k-th cell not yet classified
keep_dff = cell(num_cells, 1);

% Application loop
%------------------------------------------------------------
cell_idx = 1;
while (cell_idx <= num_cells)

    trace = ds.get_trace(cell_idx);
    [baseline, info] = polyfit_nonactive_frames(trace,...
        params.threshold, params.padding, params.order);
    params.threshold = info.threshold;
    dff_trace = (trace - baseline)./baseline;
    if ~isempty(F)
        dff_trace2 = decorrelate_trace(dff_trace, F, info.nonactive_frames);
    end
    
    if compute_zsc
        sigma = std(dff_trace(info.nonactive_frames));
        dff_trace = dff_trace/sigma;
        
        if ~isempty(F)
            sigma2 = std(dff_trace2(info.nonactive_frames));
            dff_trace2 = dff_trace2/sigma;
        end
    end
    
    % Draw results
    %------------------------------------------------------------
    subplot(gui.h_orig);
    cla; hold on;
    nf = info.nonactive_frames;
    af = ~nf; % active_frames
    
    nf_segs = frame_list_to_segments(find(nf));
    for k = 1:size(nf_segs,1)
        frames = nf_segs(k,1):nf_segs(k,2);
        plot(t(frames), trace(frames), 'Color', color);
    end
    
    af_segs = frame_list_to_segments(find(af));
    for k = 1:size(af_segs,1)
        frames = af_segs(k,1):af_segs(k,2);
        plot(t(frames), trace(frames), 'r');
    end
    
    plot(t_lims, info.threshold*[1 1], 'k--');
    plot(t, baseline, 'k-', 'LineWidth', 2);
    hold off;
    ylim(compute_ylims(trace));
    ylabel('Orig. fluorescence');
    title(sprintf('Cell %d: threshold=%.1f, padding=%d, order=%d',...
        cell_idx, params.threshold, params.padding, params.order));
    
    subplot(gui.h_dff);
    cla; hold on;
    plot(t, dff_trace, 'Color', color);
    if ~isempty(F)
        plot(t, dff_trace2, '--', 'Color', color2);
        if ~compute_zsc
            legend('DFF', 'Decorrelated DFF', 'Location', 'NorthWest');
        else
            legend('zsc', 'Decorrelated zsc', 'Location', 'NorthWest');
        end
    end
%     plot(t_lims, [0 0], 'k-', 'LineWidth', 2);
    hold off;
    xlim(t_lims)
    xlabel(sprintf('Time (s); FPS = %.1f Hz', fps));
    ylim(compute_ylims(dff_trace));
    if ~compute_zsc
        ylabel('\DeltaF/F');
    else
        ylabel('z-score');
    end
    
    zoom xon;
    
    % Ask the user to classify the cell candidate
    %------------------------------------------------------------
    prompt = sprintf('Get DFF traces (%d/%d; order=%d) >> ',...
        cell_idx, num_cells, params.order);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid index and jump to it
        if ((1 <= val) && (val <= num_cells))
            cell_idx = val;
            params = default_params;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp)
            go_to_next_cell();
        else
            switch (resp(1))
%                 case 'c' % Accept DFF result
%                     dff_traces(:,cell_idx) = dff_trace;
%                     baseline_fit_infos{cell_idx} = info;
%                     keep_dff{cell_idx} = true;
% 
%                     go_to_next_cell();
                    
                case {'c', 'd'} % Accept decorrelated DFF result
                    if ~isempty(F)
                        fprintf('  Using decorrelated DFF trace for Cell %d\n', cell_idx);
                        dff_traces(:,cell_idx) = dff_trace2;
                    else
                        dff_traces(:,cell_idx) = dff_trace;
                    end
                    baseline_fit_infos{cell_idx} = info;
                    keep_dff{cell_idx} = true;

                    go_to_next_cell();

                case 'n' % Reject DFF
                    fprintf('  Cell %d will be omitted from DFF output\n', cell_idx);
                    keep_dff{cell_idx} = false;
                    go_to_next_cell();

                % Set baseline fitting parameters
                %------------------------------------------------------------
                case 'o' % Set order
                    try
                        val = str2double(resp(2:end));
                        if (~isnan(val)) % Is a number
                            if (val > 0)
                                params.order = val;
                            end
                        end
                    catch
                        fprintf('  Could not parse order command!\n');
                    end
                    
                case 'p' % Set padding
                    try
                        val = str2double(resp(2:end));
                        if (~isnan(val)) % Is a number
                            if (val > 0)
                                params.padding = val;
                            end
                        end
                    catch
                        fprintf('  Could not parse order command!\n');
                    end

                case 't' % Set threshold
                    fprintf('  Plese select a new threshold on the ORIG trace\n');
                    while (1)
                        [~, params.threshold] = ginput(1);
                        if (gca == gui.h_orig)
                            break;
                        else
                            fprintf('  Error! New threshold must be defined on the ORIG trace\n');
                        end
                    end

                case 'q' % Exit
                    close(gui.fig_handle);
                    break;

            end % switch(resp)
        end % isempty(resp)
    end
end

% Save as Rec file
%------------------------------------------------------------

% First, collect ALL spatial filters and traces
[height, width] = size(ds.cell_map_ref_img);
filters = zeros(height, width, num_cells, 'single');
for k = 1:num_cells
    filters(:,:,k) = ds.cells(k).im;
end

% Unlabeled entries in 'keep_dff' will NOT be exported
for k = 1:num_cells
    if isempty(keep_dff{k})
        keep_dff{k} = false;
    end
end
keep_dff = cell2mat(keep_dff);
num_exported_cells = sum(keep_dff);

filters = filters(:,:,keep_dff);
traces = dff_traces(:,keep_dff);

info.type = 'get_dff_traces';
info.num_pairs = num_exported_cells;

info.baseline_infos = baseline_fit_infos(keep_dff);
info.fps = fps;

[rec_savename, timestamp] = save_rec(info, filters, traces);
class_savename = generate_class_file(num_exported_cells, 'timestamp', timestamp);
fprintf('%s: %d filter-trace pairs saved to "%s"\n', datestr(now), num_exported_cells, rec_savename);

    % Auxiliary functions
    %------------------------------------------------------------    
    function go_to_next_cell()
        unprocessed = cellfun(@isempty, keep_dff);
        next_idx = find_next_cell_to_process(cell_idx, unprocessed);
        
        if isempty(next_idx)
            cprintf([0 0.5 0], '  All DFF traces have been computed!\n');
            cell_idx = cell_idx + 1;
            if cell_idx > num_cells
                cell_idx = 1;
            end
        else
            cell_idx = next_idx;
        end
        
        % Reset default fitting method on advance
        params = default_params;
    end % go_to_next_cell
    
end % get_dff_traces

function y_lims = compute_ylims(tr)
    m = min(tr);
    M = max(tr);
    y_lims = [m M] + 0.1*(M-m)*[-1 1];
end % compute_ylims

function tr_dc = decorrelate_trace(tr, mu, frames_for_corr)
    tr2 = tr(frames_for_corr);
    mu2 = mu(frames_for_corr);
    p = polyfit(mu2, tr2, 1); % Linear fit
    tr_dc = tr - polyval(p, mu);
end % decorrelate_trace