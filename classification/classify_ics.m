function classify_ics(sources, ic_dir)
% Perform manual classification of ICA results
%
% Usage:
%   classify_ics(sources, 'ica001')
% where sources is a struct with the following fields:
%   sources.maze: Output from the plus maze (TXT)
%   sources.miniscope: Miniscope recording (TIF or HDF5)
%   sources.fps:  Frame rate of the recording
%   

% Specify ICA
%------------------------------------------------------------
ica_source = get_most_recent_file(ic_dir, 'ica_*.mat');
load(ica_source, 'ica_info', 'ica_filters', 'ica_traces');

% Load data
%------------------------------------------------------------
fprintf('  %s: Loading "%s" to memory...\n', datestr(now), sources.miniscope);
M = load_movie(sources.miniscope);
fprintf('  %s: Done!\n', datestr(now));

fps = sources.fps;
time = 1/fps*((1:size(ica_traces,1))-1); %#ok<*NODEF>
num_ics = ica_info.num_ICs;

% Compute a common scaling for the movie
movie_clim = compute_movie_scale(M);
fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
    datestr(now), movie_clim(1), movie_clim(2));

% Get trial info from maze output (consistent as of Cohort 9)
trial_indices = get_trial_frame_indices(sources.maze);
trial_indices = trial_indices(:,[1 end]); % [Start end]

% Classify ICs
%------------------------------------------------------------
output_name = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));
class = cell(num_ics, 1);

ic_idx = 1;
while (ic_idx <= num_ics)
    % Load IC
    ic_filter = ica_filters(:, :, ic_idx);
    trace = ica_traces(:, ic_idx);
    
    view_trace(time, trace, trial_indices);
    title(sprintf('IC %d of %d', ic_idx, num_ics));
    
    % Ask the user to classify the IC
    prompt = sprintf('Classifier (IC %d/%d) >> ', ...
                        ic_idx, num_ics);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid IC and jump to it
        if ((1 <= val) && (val <= num_ics))
            ic_idx = val;
        else
            fprintf('  Sorry, %d is not a valid trace index\n', val);
        end
    else
        resp = lower(resp);
        switch (resp)
            % Classication options
            %------------------------------------------------------------
            case {'p', 'c'} % Cell
                [~, movie_clim] = view_ic_over_movie_interactively(ic_filter, time, trace, M, fps, movie_clim);
                resp2 = input(sprintf('  Confirm classification ("%s") >> ', resp), 's');
                resp2 = lower(strtrim(resp2));
                if (strcmp(resp, resp2))
                    switch (resp2)
                        case 'p'
                            class{ic_idx} = 'phase-sensitive cell';
                        case 'c'
                            class{ic_idx} = 'cell';
                    end
                    fprintf('  Trace %d classified as %s\n', ic_idx, class{ic_idx});
                    ic_idx = ic_idx + 1;
                end

            case 'n' % Not a cell
                class{ic_idx} = 'not a cell';
                fprintf('  Trace %d classified as %s\n', ic_idx, class{ic_idx});
                ic_idx = ic_idx + 1;
                
            % Application options
            %------------------------------------------------------------
            case 'q' % Exit
                break;
            case 's' % Save classification
                save_classification(class, output_name);
                fprintf('  Saved classification result to %s\n', output_name);
            case 'l' % Load previous classification
                [file, path] = uigetfile('*.txt', 'Select existing classification');
                if (file)
                    full_file = fullfile(path, file);
                    class = load_classification(full_file);
                    fprintf('  Loaded classification from "%s"\n', file);
                end
                ic_idx = find(strcmp(class,''),1); % Go to first unlabeled IC
            case 't' % "Take" screenshot
                screenshot_name = sprintf('ic%03d.png', ic_idx);
                screenshot_name = fullfile(ic_dir, screenshot_name);
                print('-dpng', screenshot_name);
                fprintf('  Plot saved to %s\n', screenshot_name);
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

% Save at end!
save_classification(class, output_name);
