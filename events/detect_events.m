function detect_events(ds, M, fps, varargin)
% A convenience wrapper around 'detect_events_interactively.m'. Mirrors the
% design of 'classify_cells.m'

state = struct('fig_handle', [],...
               'movie_clim', []);

state.movie_clim = compute_movie_scale(M);

end % detect_events