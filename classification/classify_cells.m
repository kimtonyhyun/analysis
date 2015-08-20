function M = classify_cells(sources, rec_dir, varargin)
% Perform manual classification of candidate filter/trace pairs
%
% Usage:
%   classify_cells(sources, 'rec001')
%
% where sources is a struct with the following required fields:
%   sources.maze: Output from the plus maze (TXT)
%   sources.miniscope: Miniscope recording (TIF or HDF5)
%   sources.fps:  Frame rate of the recording
%
% and 'rec001' is a directory that contains a 'rec_*.mat' file which has
%   the filters and traces to be classified.
%
% Optional arguments allow for:
%   'movie': Use provided movie, rather than loading from disk
%

movie_provided = 0;
for k = 1:length(varargin)
    if ischar(varargin{k})
        switch lower(varargin{k})
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

% Compute a common scaling for the movie
movie_clim = compute_movie_scale(M);
fprintf('  %s: Movie will be displayed with fixed CLim = [%.3f %.3f]...\n',...
    datestr(now), movie_clim(1), movie_clim(2));

% Load filter/trace pairs to be classified
ds = DaySummary(sources.maze, rec_dir);
num_candidates = ds.num_cells;
fprintf('  %s: Loaded filters/traces from "%s"\n', datestr(now), rec_dir);

trial_indices = ds.trial_indices(:, [1 end]); % [Start end]
assert(size(M,3) == trial_indices(end,end),...
       'Error: Number of frames in movie does not match trial index table!');

% Begin classification
%------------------------------------------------------------
output_name = sprintf('class_%s.txt', datestr(now, 'yymmdd-HHMMSS'));
class = cell(num_candidates, 1);

cell_idx = 1;
while (cell_idx <= num_candidates)
    subplot(3,2,[1 2]);
    ds.plot_trace(cell_idx);
    title(sprintf('Candidate %d of %d', cell_idx, num_candidates));
    
    subplot(3,2,[3 5]);
    ds.plot_superposed_trials(cell_idx);
    
    subplot(3,2,[4 6]);
    ds.plot_cell_raster(cell_idx);
    
    % Ask the user to classify the IC
    prompt = sprintf('Classifier (%d/%d) >> ', ...
                        cell_idx, num_candidates);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number. Check if it is a valid IC and jump to it
        if ((1 <= val) && (val <= num_candidates))
            cell_idx = val;
        else
            fprintf('  Sorry, %d is not a valid trace index\n', val);
        end
    else
        resp = lower(resp);
        switch (resp)
            % Classication options
            %------------------------------------------------------------
            case {'p', 'c'} % Cell
                [~, movie_clim] = view_cell_interactively(ds, cell_idx,...
                                    M, sources.fps, movie_clim);
                resp2 = input(sprintf('  Confirm classification ("%s") >> ', resp), 's');
                resp2 = lower(strtrim(resp2));
                if (strcmp(resp, resp2))
                    switch (resp2)
                        case 'p'
                            class{cell_idx} = 'phase-sensitive cell';
                        case 'c'
                            class{cell_idx} = 'cell';
                    end
                    fprintf('  Trace %d classified as %s\n', cell_idx, class{cell_idx});
                    cell_idx = cell_idx + 1;
                end

            case 'n' % Not a cell
                class{cell_idx} = 'not a cell';
                fprintf('  Trace %d classified as %s\n', cell_idx, class{cell_idx});
                cell_idx = cell_idx + 1;
                
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
                    
                    cell_idx = find(strcmp(class,''),1); % Go to first unlabeled pair
                end
            case 't' % "Take" screenshot
                screenshot_name = sprintf('ic%03d.png', cell_idx);
                screenshot_name = fullfile(rec_dir, screenshot_name);
                print('-dpng', screenshot_name);
                fprintf('  Plot saved to %s\n', screenshot_name);
            otherwise
                fprintf('  Could not parse "%s"\n', resp);
        end
    end
end

% Save at end!
save_classification(class, output_name);
