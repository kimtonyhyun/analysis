function classify_ics(sources, ic_dir)

% Usage:
%   classify_ics(sources, 'ica001')

% Specify ICA
%------------------------------------------------------------
ica_source = get_most_recent_file(ic_dir, 'ica_*.mat');

load(ica_source, 'ica_info', 'ica_filters', 'ica_traces');

% Load data
%------------------------------------------------------------
fprintf('  %s: Loading "%s" to memory...\n', datestr(now), sources.miniscope);
movie = load_movie(sources.miniscope);
fprintf('  %s: Done!\n', datestr(now));

time = 1/sources.fps*((1:size(ica_traces,1))-1); %#ok<*NODEF>
num_ics = ica_info.num_ICs;

% Get trial info from maze output (consistent as of Cohort 9)
trial_frame_indices = get_trial_frame_indices(sources.maze);
compressed_indices = compress_frame_indices(trial_frame_indices, sources.trim);
compressed_indices = compressed_indices(:,[1 end]); % [Start end]

% Classify ICs
%------------------------------------------------------------
output_name = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));
class = cell(num_ics, 1);

ic_idx = 1;
while (ic_idx <= num_ics)
    % Load IC
    ic_filter = ica_filters(:, :, ic_idx);
    trace = ica_traces(:, ic_idx);
    
    % First, show the trace
    view_trace(time, trace, compressed_indices);
    title(sprintf('Trace %d of %d', ic_idx, num_ics));
    
    % Ask the user to classify the trace
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
                view_ic_over_movie_interactively(ic_filter, time, trace, movie, 100);
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
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

% Save at end!
save_classification(class, output_name);
