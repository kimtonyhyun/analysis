function extract_cells(movie_source,pca_source,varargin)

%   Cell extraction method that is based on an algorithm that essentially
%   solves the sparse source separation problem to find cell filters.
%   Traces are then found by using least squares.
%
%   Inputs:
%
%       movie_source: path to the movie matrix (height x width x frames)
%
%       pca_sources: path to the pca .mat file. It is assumed that there is
%       a variable named pca_filters that has size [# of PCs x # of pixels]
%       in the .mat file.
%       
%
%   Outputs:
%
%       Script saves the traces and filters in a .mat file.
%
%       inv_quality: Inverse quality metric that indicates how good each
%       found potential is as far as the extraction algorithm is concerned.
%       It is a [1 x # of extracted cells] array whose entries are
%       normalized to [(1/# of pixels),1] range. This array is also plotted
%       on the screen at the end of execution.
%
%   Hakan Inan (Aug 2015)
%

% Some parameters
mem_occup_scale_CPU = 0.5; % Occupy this fraction of available memory on RAM
mem_occup_scale_GPU = 0.7; % Occupy this fraction of available memory on GPU
corr_thresh = 0.2;
filter_thresh = 0.25;
max_num = 3000;
do_baseline_fix = 1;
good_cell_limit = 0.4; % threshold to call cell (wrt a custom metric, see code)
bad_cell_limit = 0.2; % threshold to call non-cell

if ~isempty(varargin)
    for k = 1:length(varargin)
        switch varargin{k}
            case 'fix_baseline'
                do_baseline_fix = varargin{k+1};
                if do_baseline_fix~=0 || do_baseline_fix~=1
                    error('fix_baseline must be 0 or 1');
                end
        end
    end
end

fprintf('%s: Loading movie singular vectors...\n',datestr(now));

pca = load(pca_source);
is_trimmed = pca.info.trim.enabled;
h = pca.info.movie_height;
w = pca.info.movie_width;
t = pca.info.movie_frames;
idx_kept = pca.info.trim.idx_kept;

fprintf('%s: Extracting filters...\n',datestr(now));

U = (pca.filters');
[N,num_pcs] = size(U);
norms_U = sqrt(sum(U.^2,2));
U_norm = bsxfun(@times,U,1./norms_U);
[Q_fine,idx_cells] = extract_filters(U,max_num);
num_extracted = size(Q_fine,2);
F = U*Q_fine;

if is_trimmed % Untrim the pixels
    F_temp = zeros(h*w,num_extracted);
    F_temp(idx_kept,:) = F;
    F = F_temp;
    clear F_temp;
    idx_cells_trimmed = idx_cells;
    idx_cells = idx_kept(idx_cells);
end

% Cell centroids
cent_cells = zeros(2,num_extracted);
[cent_cells(1,:),cent_cells(2,:)] = ind2sub([h,w],idx_cells); % [vert;horiz] format

fprintf('\t \t \t Extracted %d potential cells \n',num_extracted);
fprintf('%s: Cleaning filters and removing duplicates...\n',datestr(now));
[F,cent_cells] = modify_filters(F,filter_thresh);
F = bsxfun(@times,F,1./sum(F,1)); % Make each filter sum up to 1

fprintf('%s: Loading movie for trace extraction...\n',datestr(now));
M = load_movie(movie_source);

% Extract time traces
fprintf('%s: Extracting traces...\n',datestr(now));
M = reshape(M,h*w,t);
traces = extract_traces(M,F);

% Remove baseline from the traces
if do_baseline_fix
    for k = 1:size(F,2);
        traces(:,k) = fix_baseline(traces(:,k));
    end
end

% % Identify split cells by looking at traces
% fprintf('%s: Identifying possible split cells...\n',datestr(now));
% cc = similar_traces(traces,cent_cells);
% 
% fprintf('\t\t\t Identified %d split cells (total of %d objects), merging now...\n',...
%     length(cc),length(cell2mat(cc)));
% % Merge filters for split cells
% idx_elim = [];
% for i = 1:length(cc)
%     idx_to_merge = cc{i};
%     
%     if is_trimmed
%         filters_to_merge = U*U_norm(idx_cells_trimmed(idx_to_merge),:)';
%         F_temp = zeros(h*w,length(idx_to_merge));
%         F_temp(idx_kept,:) = filters_to_merge;
%         filters_to_merge = F_temp;
%         clear F_temp;
%     else
%         filters_to_merge = U*U_norm(idx_cells(idx_to_merge),:)';
%     end
%     
%     [~,merged_filter,~,~] = cleanup_filters(sum(filters_to_merge,2),filter_thresh);
%     merged_filter = merged_filter / sum(merged_filter);
% %     merged_filter = sum(F(:,idx_to_merge),2)/length(idx_to_merge);
%     F(:,idx_to_merge(1)) = merged_filter;
%     idx_elim = [idx_elim,idx_to_merge(2:end)];
% end
% F(:,idx_elim) = [];
% cent_cells(:,idx_elim) = [];
% 
% 
% % Extract traces again
% fprintf('%s: Extracting traces with updated filters...\n',datestr(now));
% traces = extract_traces(M,F);
clear M;

filters = reshape(F,h,w,size(F,2)); %#ok<*NASGU>

% % Sort cells wrt a metric that determines how good a trace looks
% quality_traces = goodness_trace(traces);
% [ss,idx] = sort(quality_traces,'descend');
% idx_stationary = sort(idx(ss>good_cell_limit));
% idx_move = idx(ss<=good_cell_limit);
% filters = filters(:,:,[idx_stationary,idx_move]);
% traces = traces(:,[idx_stationary,idx_move]);
% cent_cells = cent_cells(:,[idx_stationary,idx_move]);
% save('cent_cells','cent_cells');
% s1 = sprintf('Total of %d potential cells are extracted.',...
%     length(ss));
% s2 = sprintf('First %d extracted cells are likely to be cells.',...
%     sum(ss>good_cell_limit));
% s3 = sprintf('It is likely that there are few or no cells from  #%d onwards. ',...
%     idx(find(ss<=bad_cell_limit,1)));
% fprintf('%s: Cell Statistics:\n \t\t\t %s\n \t\t\t %s\n \t\t\t %s\n',...
%     datestr(now),s1,s2,s3);

% Reconstruction settings
info.type = 'PCA+SSS';
info.num_pairs = size(F,2);
info.movie_source = movie_source;
info.pca_source = pca_source; %#ok<*STRNU>

% Save the result to mat file
%------------------------------------------------------------
timestamp = datestr(now, 'yymmdd-HHMMSS');
rec_savename = sprintf('rec_%s.mat', timestamp);

fprintf('%s: Writing the output to file... \n',datestr(now));
fprintf('%s: Output will be saved in %s \n',datestr(now),rec_savename);

save(rec_savename, 'info', 'filters', 'traces','-v7.3');

fprintf('%s: All done! \n',datestr(now));



%-------------------
% Internal functions
%-------------------

function avail_mem = compute_gpu_memory()
% Available memory on GPU in bytes

    D = gpuDevice;
    % Handle 2 different versions of gpuDevice
    if isprop(D,'AvailableMemory')
        avail_mem = D.AvailableMemory;
    elseif isprop(D,'FreeMemory')
        avail_mem = D.FreeMemory;
    else 
        avail_mem = 0;
    end
end

function avail_mem = compute_cpu_memory()
% Available memory on CPU in bytes

    if ispc % Windows
        [~,sys] = memory;
        avail_mem = (sys.PhysicalMemory.Available);
    elseif ismac %Mac
        % TODO
    else % Linux
         [~,meminfo] = system('cat /proc/meminfo');
         tokens = regexpi(meminfo,'^*MemAvailable:\s*(\d+)\s','tokens');
         avail_mem = str2double(tokens{1}{1})*1024;
    end
end

function [Q_fine,idx_cells] = extract_filters(U,max_num)
% Extract filters recursively using the 1-norm criterion

    [N,num_pcs] = size(U);
    norms_U = sqrt(sum(U.^2,2));
    U_norm = bsxfun(@times,U,1./norms_U);
    
    % Brute-force approach for 1-norm computation, use GPU if available
    use_gpu = 0;
    gpu_exists = gpuDeviceCount>0;
    if gpu_exists
        avail_mem = compute_gpu_memory();    
        f = max(avail_mem/4 - N*(num_pcs+2),0); % maximum available element size 
        f = f*mem_occup_scale_GPU;
        chunk_size = floor(f/ (2*N+num_pcs)); %# of pixels at once
        use_gpu = chunk_size>=1;
    end    

    if use_gpu %GPU    
        fprintf('\t \t \t Using GPU to accelarate processing\n');
        U = gpuArray(U);
        one_norms = gpuArray(zeros(1,N,'single'));
        norms_U = gpuArray(norms_U);

        for i = 1:chunk_size:N
            fin = min(chunk_size-1,N-i);
            Q = bsxfun(@times,U(i:i+fin,:),1./norms_U(i:i+fin))';
            one_norms(i:i+fin) = sum(abs(U*Q),1);
            clear Q;
        end

        % Transfer variables back to CPU
        U = gather(U);
        one_norms = gather(one_norms);
        norms_U = gather(norms_U);

    else % CPU
        avail_mem = compute_cpu_memory();
        f = avail_mem/4 ; % maximum available element size
        f = f*mem_occup_scale_CPU;
        chunk_size = floor(f/ (2*N+num_pcs)); %# of pixels at once

        one_norms = zeros(1,N);
        for i = 1:chunk_size:N
            fin = min(chunk_size-1,N-i);
            Q = U_norm(i:i+fin,:)';
            one_norms(i:i+fin) = sum(abs(U*Q),1);
            clear Q;
        end

    end

    idx_possible = 1:N;
    F = zeros(N,max_num);
    inv_quality_filters = zeros(1,max_num);
    idx_cells = zeros(1,max_num);
    acc=0;

    while true
        acc = acc+1;
        [val,idx_this] = min(one_norms(idx_possible));
        idx_this = idx_possible(idx_this);
        sig_this = U*U_norm(idx_this,:)';
        correl = sig_this.*(1./norms_U);
        idx_possible = intersect(find(abs(correl)<corr_thresh),idx_possible);
        F(:,acc) = sig_this;
        inv_quality_filters(acc) = val;
        idx_cells(acc) = idx_this;
        % Termination condition
        if isempty(idx_possible) || acc==max_num
            break;
        end
    end

    % Truncate F in case there are less components than max_num
    F = F(:,1:acc);
    inv_quality_filters = inv_quality_filters(1:acc)/sqrt(N);
    idx_cells = idx_cells(1:acc);

    % Q is the coarse linear mapping from U to filters (UQ=F)
    Q = U_norm(idx_cells,:)';

    % Fine tuning: estimate the optimum # of pcs to use for each filter
    tic
    num_sweep = 25;
    rat_sweep = linspace(0.2,1,num_sweep);
    
    one_norms_F = zeros(acc,num_sweep);    
    for rat = rat_sweep % forward sweep
        idx_stop = ceil(num_pcs*rat);
        Q_trunc = Q(1:idx_stop,:);
        norms_Q_trunc = sqrt(sum(Q_trunc.^2,1));
        Q_trunc = bsxfun(@times,Q_trunc,1./norms_Q_trunc);
        F_trunc = U(:,1:idx_stop)*Q_trunc;
        one_norms_F(:,rat_sweep==rat) = sum(abs(F_trunc),1)';
        toc
    end

    % Reconstruct Q (Q_fine)
    [best_one_norms,idx_best_rat] = min(one_norms_F,[],2);
    idx_uncertain = find(idx_best_rat==1);
    
    % Deal with indices where its not certain which ratio is the best 
%     idx_best_rat(idx_best_rat==1) = num_sweep/2;
    for ii = idx_uncertain'
        one_norms_vs_rat = one_norms_F(ii,:)/sqrt(N);
        x = smooth(one_norms_vs_rat,floor(num_sweep/4));
        [xmax,imax,xmin,imin] = extrema(x);
        % Remove the first minimum (first index)
        xmin = xmin(2:end);
        imin = imin(2:end);
        if ~isempty(imin)
            % Find local minimum with big enough a basin
            for iii = 1:length(imin) 
                idx_candid = imin(iii);
                idx_max_before_this = find((idx_candid-imax)>0,1);
                if xmax(idx_max_before_this)- xmin(iii) > 3e-3
                    idx_best_rat(ii) = idx_candid;
                    best_one_norms(ii) = one_norms_F(ii,idx_candid);
                    break;
                end
            end
        end
    end
    
    best_rat = rat_sweep(idx_best_rat);
    idx_stop_best = ceil(num_pcs*best_rat);
    Q_fine = Q;
    for ii = 1:acc
        Q_fine(idx_stop_best(ii)+1:end,ii) = 0;
    end
    Q_fine = bsxfun(@times,Q_fine,1./sqrt(sum(Q_fine.^2,1)));

    % Eliminate possible duplicate cells by a 2nd step
    num_cells_step1 = acc;
    Corr_cells = Q_fine'*Q_fine;
    
    idx_possible = 1:num_cells_step1;
    idx_cells = zeros(1,num_cells_step1);
    acc=0;

    while true
        acc = acc+1;
        [~,idx_this] = min(best_one_norms(idx_possible));
        idx_this = idx_possible(idx_this);
        sig_this = U*Q_fine(:,idx_this);
        correl = Corr_cells(idx_this,:);
        idx_possible = intersect(find(abs(correl)<corr_thresh+0.2),idx_possible);
        idx_cells(acc) = idx_this;
        % Termination condition
        if isempty(idx_possible) || acc==max_num
            break;
        end
    end
    idx_cells = idx_cells(1:acc);
    Q_fine = Q_fine(:,idx_cells);
    elim = num_cells_step1-acc;
   
end

function [filters_out,cent_out] = modify_filters(filters_in,filter_thresh)
    
    % Clean up the filters
    [masks,filters_int,cent_int] = cleanup_filters(filters_in,filter_thresh);
    num_filters = size(filters_int,2);

    % Identify duplicate cells
    new_assignments = compute_reassignments(masks,cent_int);

    % Write to output variable (merging the duplicates)
    filters_out = zeros(h*w,num_filters);%length(unique(new_assignments)));
    cent_out = zeros(2,num_filters);%length(unique(new_assignments)));
    
    acc = 0;
    idx_merged = [];
    for idx_filt = 1:num_filters
        cells_this = find(new_assignments==idx_filt);
        len = length(cells_this);
        if true%len==1 % not an overlapping cell
            acc = acc+1;
            filters_out(:,acc) = filters_int(:,idx_filt);
            cent_out(:,acc) = cent_int(:,idx_filt);
        elseif len>1 % overlapping cells, merge
            idx_merged = [idx_merged;setdiff(cells_this,idx_filt)];
            merged_filter = sum(filters_in(:,idx_keep(cells_this)),2);
            merged_filter = merged_filter/norm(merged_filter);
            [~,merged_filter,cent_dum] = cleanup_filters(merged_filter,filter_thresh);
            if sum(merged_filter)>0
                acc = acc+1;
                filters_out(:,acc) = merged_filter;
                cent_out(:,acc) = cent_dum;
            else
                idx_merged(end+1) = idx_filt;
            end
        end
    end
    
     % Truncate output variables
    filters_out = filters_out(:,1:acc);
    cent_out = cent_out(:,1:acc);
    
    if ~isempty(idx_merged)
        fprintf('\t \t \t %d cells were merged into others \n',length(idx_merged));
    end
    
    filters_out(:,acc+1:end) = []; % Remove the 0's in the end (if any)

end

function [masks_out,filters_out,cent_out] = cleanup_filters(filters_in,filter_thresh)
    
    size_thresh = 30;
    num_filters = size(filters_in,2);
    
    % Initialize output variables
    cent_out = zeros(2,num_filters);
    filters_out = zeros(h,w,num_filters*2); % Be generous in size
    masks_out = zeros(h,w,num_filters);
    
    % Remove the filter baselines
    meds = median(filters_in,1);
    filters_in = bsxfun(@minus,filters_in,meds);

    % Reshape filters_out into 2D
    filters_in = reshape(filters_in,h,w,num_filters);

    % Threshold filters + deal with multiple boundaries    
    acc = 0;
    num_elim_objects = 0;
    for idx_filt = 1:num_filters
        this_filter = filters_in(:,:,idx_filt);
        [boundaries, ~] = compute_ic_boundary(this_filter, filter_thresh);
        num_objects = length(boundaries);
        if num_objects <20
            for idx_obj = 1:min(1,num_objects)
                mask_candid = poly2mask(boundaries{idx_obj}(:,1), boundaries{idx_obj}(:,2), h, w);
                if nnz(mask_candid)>=size_thresh
                    acc = acc+1;
                    s = regionprops(mask_candid,'centroid');
                    filters_out(:,:,acc) = mask_candid.*this_filter;
                    masks_out(:,:,acc) = mask_candid;
                    cent_out(:,acc) = s(1).Centroid;
                else
                    num_elim_objects = num_elim_objects + 1;
                end
            end
        end
    end
    
    % Truncate output variables
    filters_out = filters_out(:,:,1:acc);
    masks_out = masks_out(:,:,1:acc);
    cent_out = cent_out(:,1:acc);
    
    % Reshape filters_out back into 1D
    filters_out = reshape(filters_out,h*w,acc);
    
    % num_filters + num_extra_objects = acc
    num_extra_objects = acc-num_filters;
    if num_extra_objects>0
        fprintf('\t \t \t %d cells were created by splitting existing cells \n',num_extra_objects);
    end
%     if ~isempty(num_elim_objects)
%         fprintf('\t \t \t %d cells were eliminated due to minimum size requirement(<%d pixels)\n',...
%             num_elim_objects,size_thresh);
%     end
end

function new_assignments = compute_reassignments(masks,cent)
    
    num_filters = size(masks,3);
    dist_thresh = max(h,w)/30;
    Dist = compute_dist_matrix(cent);
    Dist = Dist < dist_thresh; % Retain only the closeby cells
    
    % Calculate overlap between closeby cells
    J_sim = cell(num_filters,1);
    for idx_filt = 1:num_filters
        closeby_cells = find(Dist(idx_filt,:)==1);
        if ~isempty(closeby_cells)
            sim = zeros(length(closeby_cells),2);
            for j = 1:length(closeby_cells)
                sim(j,:) = [closeby_cells(j),compute_overlap(masks(:,:,idx_filt),masks(:,:,closeby_cells(j)))];
            end
            J_sim{idx_filt} = sim;
        end
    end
    
    % Merge cells that have similarity above a threshold;
    sim_thresh = 0.8;
    new_assignments = (1:num_filters)';
    for idx_filt = 1:num_filters
        if ~isempty(J_sim{idx_filt}) 
            for j = 1:size(J_sim{idx_filt},1) % for all cells that overlap with this
                this_sim = J_sim{idx_filt}(j,2);
                if this_sim > sim_thresh % similarity is above threshold
                    that_cell = J_sim{idx_filt}(j,1);
                    new_assignments(new_assignments==new_assignments(that_cell)) = idx_filt;
                    J_sim{that_cell} = [];
                end
            end
        end
    end
    
end

function overlap = compute_overlap(mask1, mask2)
    intrsct = nnz(mask1 & mask2);
    siz1 = nnz(mask1); siz2 = nnz(mask2);
    overlap = intrsct / min(siz1,siz2);
end

function val = goodness_trace(traces)
%Determine the "spiky-ness" by comparing distribution to gaussian  

    num_cells = size(traces,2);
    val = zeros(1,num_cells);
    for idx_cell =1:num_cells
        tr = traces(:,idx_cell);
        tr(tr<0) = [];
        tr = [tr;-tr]; %make distribution symmetric
        tr = tr/std(tr); % standardize
        [c,x] = hist(tr,100);
        f_g = 1/sqrt(2*pi)*exp(-x.^2/2); % gaussian pdf
        f_r = c/length(tr)/(x(2)-x(1)); % empirical pdf
        val(idx_cell) = norm(f_g-f_r,2);
    end
end

function cc = similar_traces(traces,cent_cells)
% Outputs a cell array with indices to be merged

    % Mask for closeby cells
    dist_thresh = max(h,w)/40;
    dist_centroids = compute_dist_matrix(cent_cells);
    dist_centroids = dist_centroids < dist_thresh; % Retain only the closeby cells
    
    % Smooth traces
    traces = medfilt1(double(traces),5);
    
    % Normalize traces
    [num_frames,num_cells] = size(traces);
    mean_traces = sum(traces,1)/num_frames;
    centered_traces = bsxfun(@minus,traces,mean_traces);
    normed_traces = bsxfun(@times,centered_traces,1./sqrt(sum(centered_traces.^2,1)) );
    
%     %Smooth traces
%     normed_traces = medfilt1(double(normed_traces),5);
    
    % Find distance between traces
    dist_traces = compute_dist_matrix(normed_traces);
    dist_traces(dist_traces==inf)=0;
    dist_traces = dist_traces<0.8;
    
    % Connectivity matrix
    A = sparse(double(dist_traces.*dist_centroids)+eye(num_cells));
    
    % do multiple DFSs to get connected components
    cc = {};
    idx_visit = 1:num_cells;
    while true
        idx_this = idx_visit(1);
        [d,~,~,~] = dfs(A,idx_this);
        cc_this = find(d'>-1);
        if length(cc_this)>1
            cc{end+1} = cc_this;          
        end
        idx_visit = setdiff(idx_visit,cc_this);
        if isempty(idx_visit)
            break;
        end
    end
end

function Dist = compute_dist_matrix(centroids)
% Calculate the matrix of centroid distances

    num_comp = size(centroids,2);
    NN = repmat(sum(centroids.^2,1),num_comp,1);
    Dist = sqrt(max(NN+NN'-2*centroids'*centroids,0));  %#ok<MHERM>
    Dist(Dist<1e-3) = inf; % Set diagonals to infinity
    
end

end
