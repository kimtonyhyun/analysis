function eval_fit(X, Xest, scores, varargin)

assert(all(size(X)==size(Xest)),...
    'Raw and fitted data do not have the same dimensions!');

[num_neurons, samples_per_trial, num_trials] = size(X);
num_score_types = length(scores);

% Default options and process varargin
%------------------------------------------------------------
selected_cells = 1:num_neurons;
for v = 1:length(varargin)
    vararg = varargin{v};
    if ischar(vararg)
        switch lower(vararg)
            case 'cells'
                selected_cells = varargin{v+1};
        end
    end
end

% Browser loop
%------------------------------------------------------------
h = figure;
score = scores(1);

cidx = 1;

while (1);
    cell_idx = selected_cells(cidx);
    draw_cell(cell_idx);
    
    % Ask user for command
    prompt = sprintf('Evaluator (Cell %d, %s=%.4f) >> ',...
                     cell_idx, score.name, score.vals_n(cell_idx));
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if ~isnan(val) % Is a number
        cidx2 = find(selected_cells == val, 1);
        if ~isempty(cidx2)
            cidx = cidx2;
        else
            fprintf('  Sorry, %d is not a valid cell index\n', val);
        end
    else
        resp = lower(resp);
        if isempty(resp) % Empty string gets mapped to "n"
            resp = 'n';
        end
        
        % Interaction loop
        %------------------------------------------------------------
        switch resp(1)
            case 'n' % Next
                if (cidx < length(selected_cells))
                    cidx = cidx + 1;
                else
                    fprintf('  Already at final cell!\n');
                end
                
            case 'p' % Previous
                if cidx > 1
                    cidx = cidx - 1;
                else
                    fprintf('  Already at first cell!\n');
                end
                
            case 's' % Score
                if length(resp) == 1
                    fprintf('  Available scores:\n');
                    for s = 1:num_score_types
                        fprintf('  %d: %s\n', s, scores(s).name);
                    end
                else
                    val = str2double(resp(2));
                    if ~isnan(val)
                        if (1 <= val) && (val <= num_score_types)
                            score = scores(val);
                        end
                    end
                end
                    
            case 'q'
                close(h);
                break;
                
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end % while 1

    function draw_cell(cell_idx)
        tr = squeeze(X(cell_idx,:,:)); % [Samples_per_trial x Num_trials]
        tr_est = squeeze(Xest(cell_idx,:,:));
        
        tr_unwrap = tr(:);
        tr_est_unwrap = tr_est(:);
        m = min(tr_unwrap);
        M = max(tr_unwrap);
        
        % Plot the full trace with wrapped lines
        %------------------------------------------------------------
        subplot(6,5,[1 2 3 4 6 7 8 9]);
        samples_per_line = 4000;
        [frame_chunks, num_chunks] = make_frame_chunks(length(tr_unwrap), samples_per_line);
        for k = 1:num_chunks
            y_offset = -0.8*(M-m)*(k-1);
            frames = frame_chunks(k,1):frame_chunks(k,2);
            plot(tr_unwrap(frames)+y_offset);
            hold on;
            plot(tr_est_unwrap(frames)+y_offset, 'r');
        end
        hold off;
        set(gca, 'YTickLabel', '');
        grid on;
        ylim([y_offset+m M] + 0.05*(M-m)*[-1 1]);
        zoom xon;
        title(sprintf('Cell %d', cell_idx));
        
        % Plot scatter plot of raw data and fit
        %------------------------------------------------------------
        subplot(6,5,[5 10]);
        plot(tr_est_unwrap, tr_unwrap, '.');
        xlabel('Fit value');
        ylabel('Orig value');
        grid on;
        axis equal;
        xlim([m M]);
        ylim([m M]);
        hold on;
        plot([m M], [m M], 'k--');
        hold off;
        title(sprintf('Full %s: %.4f', score.name, score.vals_n(cell_idx)));
        
        % Plot individual trials
        %------------------------------------------------------------
        scores_per_trial = score.vals(cell_idx,:);
        
        % Plot 10 trials with the BEST fits
        [~, sorted_trials] = sort(scores_per_trial, score.sort_order);
        trials_to_display = sorted_trials(1:10);
        for k = 1:length(trials_to_display)
            subplot(6,5,10+k);
            draw_trial(trials_to_display(k), [0 0.5 0]); % Draw in GREEN
        end
        
        % Plot 10 trials with the WORST fits
        trials_to_display = sorted_trials(end-9:end);
        trials_to_display = fliplr(trials_to_display); % Show worst first
        for k = 1:length(trials_to_display)
            subplot(6,5,20+k);
            draw_trial(trials_to_display(k), [1 0 0]); % Draw in RED
        end
        
        function draw_trial(trial_idx, color)
            plot(tr(:,trial_idx),'.');
            hold on;
            plot(tr_est(:,trial_idx), 'Color', color, 'LineWidth', 2);
            hold off;
            grid on;
            xlim([1 samples_per_trial]);
            ylim([m M]);
            set(gca, 'XTickLabel', '');
            set(gca, 'YTickLabel', '');
            ylabel(sprintf('Trial %d', trial_idx));
        end % draw_trial
        
    end % draw_cell

end % evaluator