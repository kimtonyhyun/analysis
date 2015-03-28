function M = classify_cells(sources, ic_dir, varargin)
% Perform manual classification of candidate filter/trace pairs
%
% Usage:
%   classify_ics(sources, 'ica001')
%
% where sources is a struct with the following fields:
%   sources.maze: Output from the plus maze (TXT)
%   sources.miniscope: Miniscope recording (TIF or HDF5)
%   sources.fps:  Frame rate of the recording
%
% and 'ica001' is a directory that contains 'ica_*.mat' or 'rec_*.mat'
%   depending on whether one is classifying ICA or reconstructed pairs.
%
% Optional arguments allow for:
%   'reconst': Use reconstructed filter/trace pairs, rather than ICA
%   'movie': Use provided movie, rather than loading from disk
%

movie_provided = 0;
use_reconstruction = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
            case 'reconst'
                use_reconstruction = 1;
            case 'movie'
                movie_provided = 1;
                M = varargin{k+1}; % Note Matlab's lazy eval
        end
    end
end

% Load movie
if movie_provided
    fprintf('  %s: Using provided movie matrix\n', datestr(now));
else
    fprintf('  %s: Loading "%s" to memory...\n', datestr(now), sources.miniscope);
    M = load_movie(sources.miniscope);
    fprintf('  %s: Done!\n', datestr(now));
end

num_frames = size(M, 3);
fps = sources.fps;
time = 1/fps*((1:num_frames)-1); %#ok<*NODEF>

% Compute a common scaling for the movie
movie_clim = compute_movie_scale(M);
fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
    datestr(now), movie_clim(1), movie_clim(2));

% Load filter/trace pairs to be classified
if use_reconstruction
    ds = DaySummary(sources.maze, ic_dir, 'reconst');
else
    ds = DaySummary(sources.maze, ic_dir);
end
num_candidates = ds.num_cells;
fprintf('  %s: Loaded filters/traces from "%s"\n', datestr(now), ic_dir);

trial_indices = ds.trial_indices(:, [1 end]); % [Start end]
assert(num_frames == trial_indices(end,end),...
       'Number of frames in movie does not match trial index table!');

% Begin classification
%------------------------------------------------------------
output_name = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));
class = cell(num_candidates, 1);

ic_idx = 1;
while (ic_idx <= num_candidates)
    % Load IC
    filter = ds.cells(ic_idx).im;
    trace = ds.get_trace(ic_idx);
    
    view_trace(time, trace, trial_indices);
    title(sprintf('Candidate %d of %d', ic_idx, num_candidates));
    
    % Ask the user to classify the IC
    prompt = sprintf('Classifier (%d/%d) >> ', ...
                        ic_idx, num_candidates);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid IC and jump to it
        if ((1 <= val) && (val <= num_candidates))
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
                [~, movie_clim] = view_ic_over_movie_interactively(filter, time, trace, M, fps, movie_clim);
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
                ic_idx = find(strcmp(class,''),1); % Go to first unlabeled pair
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
