function backapply_filters(ds, movie_src)

movie_dataset = '/Data/Images';
movie_size = get_dataset_info(movie_src, movie_dataset);
height = movie_size(1);
width = movie_size(2);
num_frames = movie_size(3);

% Collect spatial filters from DaySummary
cell_indices = find(ds.is_cell);
num_cells = length(cell_indices);
filters = zeros(num_cells, height*width);
for k = 1:num_cells
    cell_idx = cell_indices(k);
    filters(k,:) = reshape(ds.cells(cell_idx).im, 1, height*width);
end

% Compute traces by streaming chunks of movie from disk
traces = zeros(num_cells, num_frames);

frame_chunk_size = 2500;
[frame_chunks, num_chunks] = make_frame_chunks(num_frames, frame_chunk_size);

for i = 1:num_chunks
    fprintf('%s: Processing frames %d to %d (out of %d)...\n',...
        datestr(now), frame_chunks(i,1), frame_chunks(i,2), num_frames);
    
    chunk_start = frame_chunks(i,1);
    chunk_count = frame_chunks(i,2) - frame_chunks(i,1) + 1;

    movie_chunk = h5read(movie_src, movie_dataset,...
                         [1 1 chunk_start],...
                         [height width chunk_count]);
                     
    movie_chunk = reshape(movie_chunk, height*width, chunk_count);
    traces(:,frame_chunks(i,1):frame_chunks(i,2)) = ...
        filters * movie_chunk;
end

% Save results into standard Rec format
info.type = 'backapplied_filters';
info.num_pairs = num_cells;
info.backapplied.movie_src = movie_src;

filters = reshape(filters, num_cells, height, width);
filters = permute(filters, [2 3 1]); % [height x width x num_cells]
traces = traces'; % [num_frames x num_cells]

save_rec(info, filters, traces);