function extractor(movie_source,opts)
% Automated cell extraction by source separation algorithm. More
% documentation to come.

if ~exist('opts','var')
    opts = [];
end

% Defaults
if ~isfield(opts,'num_partition_x'), opts.num_partition_x=1; end
if ~isfield(opts,'num_partition_y'), opts.num_partition_y=1; end
if ~isfield(opts,'overlap_x'), opts.overlap_x=25; end
if ~isfield(opts,'overlap_y'), opts.overlap_y=25; end
if ~isfield(opts,'svd_num_comp'), opts.svd_num_comp=1000; end
if ~isfield(opts,'spat_medfilt_enabled'), opts.spat_medfilt_enabled=1; end
if ~isfield(opts,'spat_linfilt_enabled'), opts.spat_linfilt_enabled=1; end
if ~isfield(opts,'spat_medfilt_halfwidth'), opts.spat_medfilt_halfwidth=1; end
if ~isfield(opts,'spat_linfilt_halfwidth'), opts.spat_linfilt_halfwidth=4; end
if ~isfield(opts,'temporal_smooth_enabled'), opts.temporal_smooth_enabled=1; end
if ~isfield(opts,'temporal_smooth_filter'), opts.temporal_smooth_filter='gauss'; end
if ~isfield(opts,'temporal_smooth_len'), opts.temporal_smooth_len=21; end
if ~isfield(opts,'trim_pixels'), opts.trim_pixels = 0.4; end
if ~isfield(opts,'subsample_time'), opts.subsample_time = 'off'; end
if ~isfield(opts,'subsample_time_ratio'), opts.subsample_time_ratio = 0.5; end
if ~isfield(opts,'verbose'), opts.verbose = 1; end
if ~isfield(opts,'disableGPU'), opts.disableGPU = 0; end
if ~isfield(opts,'mem_occup_ratio_CPU'), opts.mem_occup_ratio_CPU=0.8; end
if ~isfield(opts,'mem_occup_ratio_GPU'), opts.mem_occup_ratio_GPU=0.7; end
if ~isfield(opts,'existing_SVD_file'), opts.existing_SVD_file = ''; end
if ~isfield(opts,'save_SVD'), opts.save_SVD = 1; end
if ~isfield(opts,'ss_num_comp'), opts.ss_num_comp = 0; end
if ~isfield(opts,'ss_max_num_sources'), opts.ss_max_num_sources =1500; end
if ~isfield(opts,'ss_cell_size_threshold'), opts.ss_cell_size_threshold = 0; end
if ~isfield(opts,'remove_duplicate_cells'), opts.remove_duplicate_cells = 0; end

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

if opts.save_SVD
    timestamp = datestr(now, 'yymmdd-HHMMSS');
    svd_savename = sprintf('svd_%s.mat', timestamp);
end

skip_svd = exist(opts.existing_SVD_file,'file');
if ~skip_svd
    if opts.save_SVD
        save(svd_savename,'opts');
    end
    
    dispfun(sprintf('%s: Loading %s...\n', datestr(now), movie_source),opts.verbose~=0);
    M = load_movie(movie_source);
    
    % Mean subtract in time
    mean_M = mean(M,3);
    is_zeromean = std(mean_M(:)) <1e-5;
    if ~is_zeromean
        M = bsxfun(@minus,M,mean_M);
    end
    
    [height, width, num_frames] = size(M);
    num_pixels = height * width;
    
    %Uniform subsampling in time
    if strcmp(opts.subsample_time,'uniform')
        samp_freq = floor(1/opts.subsample_time_ratio);
        idx_end = floor((num_frames-1)/samp_freq)*samp_freq+1;
        idx_keep = 1:samp_freq:idx_end;
        idx_trash = setdiff(1:num_frames,idx_keep);
        M(:,:,idx_trash) = [];
        num_frames = length(idx_keep);
    end
    
    % Spatial smoothing
    if opts.spat_medfilt_enabled || opts.spat_linfilt_enabled
        dispfun(sprintf('%s: Beginning spatial filtering...\n', datestr(now)),opts.verbose~=0);
        medfilt_neighborhood = (1+2*opts.spat_medfilt_halfwidth)*[1 1];
        linfilt_neighborhood = (1+2*opts.spat_linfilt_halfwidth)*[1 1];
        for idx_frame = 1:num_frames
            frame = M(:,:,idx_frame);
            if opts.spat_medfilt_enabled
                frame = medfilt2(frame, medfilt_neighborhood);
            end
            if opts.spat_linfilt_enabled
                frame = wiener2(frame,linfilt_neighborhood);
            end
            M(:,:,idx_frame) = frame;
            if mod(idx_frame,1000)== 0
                dispfun(sprintf('%s: Filtered %d frames (out of %d)...\n',...
                    datestr(now),idx_frame, num_frames),opts.verbose==2);
            end
        end
        dispfun(sprintf('%s: Finished spatial filtering.\n', datestr(now)),opts.verbose~=0);
    end
    
    % Compute pixels to retain in case of trimming
    if opts.trim_pixels
        max_proj = max(M,[],3);
        max_proj = medfilt2(max_proj,[5,5]);
        mask_retained = max_proj>quantile(max_proj(:),opts.trim_pixels);
    end

    % Make M 2D
    M = reshape(M, num_pixels, num_frames);

    % Temporal smoothing
    if opts.temporal_smooth_enabled
        if strcmp(opts.subsample_time,'uniform') % Subsampling was done before
            filt_len = ceil(opts.temporal_smooth_len*opts.subsample_time_ratio);
        else
            filt_len = opts.temporal_smooth_len;
        end
        dispfun(sprintf('%s: Smoothing in time...\n',datestr(now)),opts.verbose~=0);
        if strcmp(opts.temporal_smooth_filter,'gauss')
            idx_filt = -floor(filt_len/2):floor(filt_len/2);
            filt = normpdf(idx_filt,0,2);
            filt = filt/sum(filt);
            mem_factor = 2; % Need this x movie size for gauss filter
        else
            mem_factor = 5; % Need this x movie size for wiener filter
        end
        
        usable_mem_CPU = compute_cpu_memory()*opts.mem_occup_ratio_CPU;
        blocksize_M = floor(usable_mem_CPU/4 / num_frames /mem_factor);
        
        if blocksize_M >= num_pixels
           if strcmp(opts.temporal_smooth_filter,'gauss')
               M= conv2(1,filt,M,'same');
           else
               M = wiener(M,filt_len);
           end
                 
        else
            for idx_begin = 1:blocksize_M:num_pixels 
                idx_end = min(idx_begin+blocksize_M-1,num_pixels);
                if strcmp(opts.temporal_smooth_filter,'gauss')
                    M(idx_begin:idx_end,:)= conv2(1,filt,M(idx_begin:idx_end,:),'same');
                else
                    M(idx_begin:idx_end,:)= wiener(M(idx_begin:idx_end,:),filt_len);
                end
            end
        end
       
        dispfun(sprintf('%s: Finished smoothing in time.\n',datestr(now)),opts.verbose~=0);
    end
   
    % Randomized subsampling in time (from Drineas et.al.)
    if strcmp(opts.subsample_time,'random')
        p = zeros(1,num_frames);
        for i =1:num_frames
            p(i) = sum(M(:,i).^2);
        end
        p = p/sum(p);
        num_keep = ceil(opts.subsample_time_ratio*num_frames);
        idx_keep = randsample(num_frames,num_keep,true,p);
        M = M(:,idx_keep);
        M = bsxfun(@times,M,1./sqrt(num_keep*p(idx_keep)));
        num_frames = num_keep;       
    end

else
    dispfun(sprintf('%s: Using existing SVD source, skipping movie preprocessing and SVD...\n',...
        datestr(now)),opts.verbose~=0);
    load(opts.existing_SVD_file,'SVD1');
    height = SVD1.info.movie_height;
    width = SVD1.info.movie_width;
    num_frames = SVD1.info.movie_frames;
    num_pixels = height * width;
    clear SVD1;
end

%%%%%%%%%%%%%%%%
% Run extraction
%%%%%%%%%%%%%%%%

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
        dispfun(sprintf('%s: Running partition %d (out of %d)...\n',datestr(now),...
            part_no,num_partitions),opts.verbose~=0);
        
        if ~skip_svd % SVD
            x_begin = (idx_partition_x-1)*(blocksize_x-ox)+1;
            x_end = min(x_begin+blocksize_x-1,width);
            x_keep = x_begin:x_end;
            y_begin = (idx_partition_y-1)*(blocksize_y-oy)+1;
            y_end = min(y_begin+blocksize_y-1,height);
            y_keep = y_begin:y_end;
            dum = zeros(height,width);
            dum(y_keep,x_keep)=1;
            if opts.trim_pixels
                dum = dum.*mask_retained;
            end
            idx_kept = find(dum);
            M_small = M(idx_kept,:);
            if num_partitions == 1
                clear M;
            end

            % Make M_small zero mean in the spatial dimension
            mean_M_small = mean(M_small,1);

            dispfun(sprintf('%s: Computing SVD...\n',datestr(now)),opts.verbose~=0);
            % Do SVD
            M_small = bsxfun(@minus, M_small, mean_M_small);
            [U,S,V] = svds_custom(M_small, opts.svd_num_comp,opts.verbose==2);
            S = diag(S);       
            
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
            if opts.save_SVD % Save SVD results
                this_savename = sprintf('SVD%d',part_no);
                savedum.(this_savename) = svd_source;
                save(svd_savename, '-struct','savedum','-append','-v7.3');
            end
            clear M_small;
        end
        
        % Source extraction
        if skip_svd
            svd = load(opts.existing_SVD_file,['SVD',num2str(part_no)]);
            svd_source = eval(['svd.SVD',num2str(part_no)]);
        end
        dispfun(sprintf('%s: Extracting cell spatial weights...\n',datestr(now)),opts.verbose~=0);
        [F_this,cents_this,inv_qualities_this] = compute_spatial_weights(svd_source,opts);
        F = [F,F_this];
        cents = [cents,cents_this];
        inv_qualities = [inv_qualities,inv_qualities_this];
        clear F_this cents_this;
    end
end

dispfun(sprintf('%s: Total of %d cells are found.\n',...
    datestr(now),size(F,2)),opts.verbose~=0);

if opts.remove_duplicate_cells
    dispfun(sprintf('%s: Removing duplicate cells...\n',...
    datestr(now),size(F,2)),opts.verbose~=0);

    [F,idx_retained] = remove_duplicate_cells(F,cents,[height,width]);
    inv_qualities = inv_qualities(idx_retained);

    dispfun(sprintf('%s: %d cells were retained after removing duplicates.\n',...
        datestr(now),size(F,2)),opts.verbose~=0);
end

% Sort spatial weights in descending goodness
[~,idx_sort] = sort(inv_qualities);
F = F(:,idx_sort);
    
clear M;
dispfun(sprintf('%s: Loading the movie again for extracting temporal traces...\n',...
    datestr(now)),opts.verbose~=0);
M = load_movie(movie_source);
M = reshape(M,height*width,num_frames);

dispfun(sprintf('%s: Extracting traces...\n',datestr(now)),opts.verbose~=0);
traces = extract_traces(M,F);


% Remove baseline from the traces
for k = 1:size(F,2);
    traces(:,k) = fix_baseline(traces(:,k));
end

filters = reshape(F,height,width,size(F,2)); %#ok<*NASGU>

% Some statistics
stats.inv_qualities = inv_qualities;
val = goodness_trace(traces);
stats.goodness_traces = val;

info.type = 'SS';
info.num_pairs = size(F,2);
info.movie_source = movie_source;
info.stats = stats;
if opts.save_SVD
    info.svd_source = svd_savename;
end

% Save the result to mat file
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

dispfun(sprintf('%s: Writing the output to file... \n',datestr(now)),opts.verbose~=0);
fprintf('%s: Output will be saved in %s \n',datestr(now),rec_savename);

save(rec_savename, 'info', 'filters', 'traces','-v7.3');

fprintf('%s: All done! \n',datestr(now));

%%%%%%%%%%%%%%%%%%%%
% Internal functions
%%%%%%%%%%%%%%%%%%%%

function F = wiener(G,conv_len)
% Wiener filtering in time, adapted from Matlab built-in wiener2

    localMean = conv2(1,ones(1,conv_len),G,'same') / conv_len;
    localVar = conv2(1,ones(1,conv_len),G.^2,'same') / conv_len - localMean.^2;
    noise = mean(localVar,2);

    F = G - localMean;
    G = bsxfun(@minus,localVar,noise);
    G = max(G, 0);
    localVar = bsxfun(@max,localVar, noise);
    F = F ./ localVar;
    F = F .* G;
    F = F + localMean;
end

function [filters_out,idx_retained] = remove_duplicate_cells(filters_in,cent_cells,hw)
% Outputs a cell array with indices to be merged

    h = hw(1);w = hw(2);
    thresh_A = 0.65;
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

function val = goodness_trace(traces)
%Determine the "spiky-ness" by comparing distribution to gaussian  

    num_cells = size(traces,2);
    val = zeros(2,num_cells);
    for idx_cell =1:num_cells
        tr = traces(:,idx_cell);
        tr(tr<0) = [];
        tr = [tr;-tr]; %make distribution symmetric
        tr = tr/std(tr); % standardize
        [c,x] = hist(tr,100);
        f_g = 1/sqrt(2*pi)*exp(-x.^2/2); % gaussian pdf
        f_r = c/length(tr)/(x(2)-x(1)); % empirical pdf
        val(1,idx_cell) = norm(f_g-f_r,2);
        val(2,idx_cell) = norm(f_g-f_r,1);
    end
end

function dispfun(str,state)
    if state == 1
        fprintf(str);
    end
end

end