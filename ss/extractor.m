function extractor(movie_source,opts)
% Automated cell extraction by source separation algorithm. More
% documentation to come.

if ~exist('opts','var')
    opts = [];
end

% Defaults
if ~isfield(opts,'num_partition_x'),  opts.num_partition_x=3; end
if ~isfield(opts,'num_partition_y'),  opts.num_partition_y=3; end
if ~isfield(opts,'overlap_x'),  opts.overlap_x=30; end
if ~isfield(opts,'overlap_y'),  opts.overlap_y=30; end
if ~isfield(opts,'num_comp'),  opts.num_comp=300; end
if ~isfield(opts,'spat_medfilt_enabled'),  opts.spat_medfilt_enabled=1; end
if ~isfield(opts,'spat_medfilt_halfwidth'),  opts.spat_medfilt_halfwidth=1; end
if ~isfield(opts,'temporal_smooth_enabled'),  opts.temporal_smooth_enabled=1; end
if ~isfield(opts,'temporal_smooth_len'),  opts.temporal_smooth_len=10; end

if opts.spat_medfilt_enabled
    splen = 2*opts.spat_medfilt_halfwidth+1;
else
    splen=0;
end
if opts.temporal_smooth_enabled
    tmplen = opts.temporal_smooth_len;
else
    tmplen=0;
end

svd_savename = sprintf('svd_n%dx%dy%dsp%dtmp%d.mat', opts.num_comp,...
    opts.num_partition_x,opts.num_partition_y,splen,tmplen);

skip_svd = exist(svd_savename,'file');
if ~skip_svd
    dummy_var = '';
    save(svd_savename,'dummy_var'); % Only way to make matlab save an blank .mat
    
    fprintf('%s: Loading %s...\n', datestr(now), movie_source);
    M = load_movie(movie_source);
    [height, width, num_frames] = size(M);

    % Median filter
    if opts.spat_medfilt_enabled
        medfilt_neighborhood = (1+2*opts.spat_medfilt_halfwidth)*[1 1];
        for idx_frame = 1:num_frames
            frame = M(:,:,idx_frame);
            M(:,:,idx_frame) = medfilt2(frame, medfilt_neighborhood);
            if mod(idx_frame,1000)== 0
                fprintf('%s: Median-filtered %d frames (out of %d)...\n',...
                    datestr(now),idx_frame, num_frames);
            end
        end
        fprintf('%s: Finished median filtering!\n', datestr(now));
    end

    % Make M 2D
    num_pixels = height * width;
    M = reshape(M, num_pixels, num_frames);

    % Temporal smoothing with moving average filter
    if opts.temporal_smooth_enabled
        fprintf('%s: Smoothing in time..\n',datestr(now))
        for i = 1:num_pixels
            M(i,:) = smooth(M(i,:),opts.temporal_smooth_len);
        end
        fprintf('%s: Finished smoothing in time!\n',datestr(now));
    end

else
    fprintf('%s: Existing SVD source detected, skipping movie preprocessing and SVD...\n', datestr(now));
    load(svd_savename,'SVD1');
    height = SVD1.info.movie_height;
    width = SVD1.info.movie_width;
    num_frames = SVD1.info.movie_frames;
    num_pixels = height * width;
    clear SVD1;
end

% Divide movie into partitions and run extraction

%shorthand variables
npx = opts.num_partition_x;
npy = opts.num_partition_y;
ox = opts.overlap_x;
oy = opts.overlap_y;

num_partitions = npx*npy;
blocksize_x = ceil((width+(npx-1)*ox)/npx);
blocksize_y = ceil((height+(npy-1)*oy)/npy);  

F = [];
cents = [];
inv_qualities = [];
for idx_partition_x = 1:npx
    for idx_partition_y = 1:npy
        
        part_no = sub2ind([npy,npy],idx_partition_y,idx_partition_x);        
        fprintf('%s: Running partition %d (out of %d)...\n',datestr(now),part_no,num_partitions);
        
        if ~skip_svd % SVD
            x_begin = (idx_partition_x-1)*(blocksize_x-ox)+1;
            x_end = min(x_begin+blocksize_x-1,width);
            x_keep = x_begin:x_end;
            y_begin = (idx_partition_y-1)*(blocksize_y-oy)+1;
            y_end = min(y_begin+blocksize_y-1,height);
            y_keep = y_begin:y_end;
            dum = zeros(height,width);
            dum(y_keep,x_keep)=1;
            idx_kept = find(dum);
            M_small = M(idx_kept,:);

            % Make M_small zero mean in the spatial dimension
            mean_M_small = mean(M_small,1);

            % Do SVD
            M_small = bsxfun(@minus, M_small, mean_M_small);
            [U,S,V] = svds_custom(M_small, opts.num_comp); 
            S = diag(S);
  
            % Save SVD results
            info.movie_height = height;
            info.movie_width  = width;
            info.movie_frames = num_frames;
            info.trim.enabled = 1;
            info.trim.idx_kept = idx_kept;
            info.medfilt.enabled = opts.spat_medfilt_enabled;  %#ok<*STRNU>
            info.medfilt.halfwidth = opts.spat_medfilt_halfwidth;
            svd_source.U = U;
            svd_source.V = V;
            svd_source.S = S;
            svd_source.info = info;
            this_savename = sprintf('SVD%d',part_no);
            savedum.(this_savename) = svd_source;
            save(svd_savename, '-struct','savedum','-append');
            clear M_small;
        end
        
        % Source extraction
        if skip_svd
            svd_source = num2str(part_no);
        end
        [F_this,cents_this,inv_qualities_this] = compute_spatial_weights(svd_source);
        F = [F,F_this];
        cents = [cents,cents_this];
        inv_qualities = [inv_qualities,inv_qualities_this];
        clear F_this cents_this;
    end
end

fprintf('%s: Total of %d cells are found. Removing duplicate cells...\n',datestr(now),size(F,2));
[F,idx_retained] = remove_duplicate_cells(F,cents,[height,width]);
inv_qualities = inv_qualities(idx_retained);

% Sort spatial weights in descending goodness
[~,idx_sort] = sort(inv_qualities);
F = F(:,idx_sort);

fprintf('%s: %d cells were retained after removing duplicates.\n',datestr(now),size(F,2));

clear M;
fprintf('%s: Loading the movie again for extracting temporal traces...\n',datestr(now))
M = load_movie(movie_source);

fprintf('%s: Extracting traces...\n',datestr(now));
traces = extract_traces(M,F);

% Remove baseline from the traces
for k = 1:size(F,2);
    traces(:,k) = fix_baseline(traces(:,k));
end

filters = reshape(F,height,width,size(F,2)); %#ok<*NASGU>

info.type = 'SS';
info.num_pairs = size(F,2);
info.movie_source = movie_source;
info.svd_source = svd_savename; 

% Save the result to mat file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

fprintf('%s: Writing the output to file... \n',datestr(now));
fprintf('%s: Output will be saved in %s \n',datestr(now),rec_savename);

save(rec_savename, 'info', 'filters', 'traces','-v7.3');

fprintf('%s: All done! \n',datestr(now));

%-------------------
% Internal functions
%-------------------

function [filters_out,idx_retained] = remove_duplicate_cells(filters_in,cent_cells,hw)
% Outputs a cell array with indices to be merged

    h = hw(1);w = hw(2);
    thresh_A = 0.7;
    num_cells = size(filters_in,2);
    masks = filters_in>0;  
    cell_sizes = sum(masks,1);
    
    % Mask for closeby cells
    dist_thresh = max(h,w)/40;
    dist_centroids = compute_dist_matrix(cent_cells);
    dist_centroids = dist_centroids < dist_thresh; % Retain only the closeby cells
    
    % Find overlap between closeby cells
    A = zeros(num_cells,num_cells); % Weighted adjacency matrix
    for y = 1:num_cells-1
        for x = y+1:num_cells
            if dist_centroids(y,x)>0
                A(y,x) = comput_overlap(masks(:,y),masks(:,x));
            end
        end
    end
    
    % Apply threshold to binarize A
    A = A > thresh_A;
    A = sparse(double(A+A'+eye(num_cells)));
    
    % do multiple DFSs to get connected components
    idx_elim = [];
    idx_visit = 1:num_cells;
    while true
        idx_this = idx_visit(1);
        [d,~,~,~] = dfs(A,idx_this);
        cc_this = find(d'>-1);
        if length(cc_this)>1 
            these_sizes = cell_sizes(cc_this);
            [~,idx_sort] = sort(these_sizes,'descend');
            
            % Keep the biggest one only
            idx_elim = [idx_elim,cc_this(idx_sort(2:end))];
        end
        idx_visit = setdiff(idx_visit,cc_this);
        if isempty(idx_visit)
            break;
        end
    end
    filters_out = filters_in;
    filters_out(:,idx_elim) = [];
    idx_retained = setdiff(1:num_cells,idx_elim);
end

function overlap = comput_overlap(mask1, mask2)
    intrsct = nnz(mask1 & mask2);
    siz1 = nnz(mask1); siz2 = nnz(mask2);
    overlap = intrsct / min(siz1,siz2);
end

function Dist = compute_dist_matrix(centroids)
% Calculate the matrix of centroid distances

    num_comp = size(centroids,2);
    NN = repmat(sum(centroids.^2,1),num_comp,1);
    Dist = sqrt(max(NN+NN'-2*centroids'*centroids,0));  %#ok<MHERM>
    Dist(Dist<1e-3) = inf; % Set diagonals to infinity
    
end

end