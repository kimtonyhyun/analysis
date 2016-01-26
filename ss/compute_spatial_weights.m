function [F,cents,inv_qualities] = compute_spatial_weights(svd_source,opts)
%   Cell extraction method that is based on an algorithm that essentially
%   solves the sparse source separation problem to find cell spatial weights.

% Some parameters
mag_thresh = 0.05;

verb1=1;
verb2=0;
if opts.verbose==0
    verb1=0;
elseif opts.verbose==2
    verb2=1;
end

U = svd_source.U;
info = svd_source.info;

if opts.ss_num_comp>0 % number of components are set
    num_comp = opts.ss_num_comp;
else % set to #PCs
   num_comp = opts.svd_num_comp;
end
U = U(:,1:num_comp);

h = info.movie_height;
w = info.movie_width;
idx_kept = info.trim.idx_kept;

[Q_fine,inv_qualities] = extract_sources(U,opts.ss_max_num_sources);
num_extracted = size(Q_fine,2);
F = U*Q_fine;

% Untrim pixels
F_temp = zeros(h*w,num_extracted);
F_temp(idx_kept,:) = F;
F = F_temp;
clear F_temp;
% save('F','F');

dispfun(sprintf('%s: Checking cells and doing some cleaning...\n',datestr(now)),verb2);
[F,cents,idx_retained] = cleanup_sources(F,mag_thresh);
inv_qualities = inv_qualities(idx_retained);

F = bsxfun(@times,F,1./max(F,[],1)); % Normalize to recover df/f at traces

dispfun(sprintf('%s: Done. Source count: %d.\n',datestr(now),size(F,2)),verb1);

%%%%%%%%%%%%%%%%%%%%
% Internal functions
%%%%%%%%%%%%%%%%%%%%

function [Q,goodness_sigs] = extract_sources(U,max_num)
% Extract sources recursively using the 1-norm criterion

    [N,num_pcs] = size(U);
    norms_U = sqrt(sum(U.^2,2));
    U_norm = bsxfun(@times,U,1./norms_U);
    
    % Brute-force approach for 1-norm computation, use GPU if available
    gpu_exists = gpuDeviceCount>0;
    if gpu_exists && (~opts.disableGPU)
        try
            avail_mem = compute_gpu_memory();
        catch
            warning('GPU was detected but available memory could not be retrieved. Using CPU...');
            avail_mem = 0;
        end
        f = max(avail_mem/4 - N*(num_pcs+2),0); % maximum available element size 
        f = f*opts.mem_occup_ratio_GPU;
        
        chunk_size = floor(f/ (2*N+num_pcs)/2); %# of pixels at once
        use_gpu = chunk_size>=1;
    else
        use_gpu = 0;
    end    
    
    if use_gpu %GPU
        dispfun(sprintf('\t\t\t GPU detected. Using it for computation...\n'),verb2);
        U = gpuArray(U);
        one_norms = gpuArray(zeros(1,N,'single'));
        norms_U = gpuArray(norms_U);

        for i = 1:chunk_size:N
            fin = min(chunk_size-1,N-i);
            Q = bsxfun(@times,U(i:i+fin,:),1./norms_U(i:i+fin))';
            %new
%             dum = U*Q;
%             stds = 3*std(dum,0,1);
%             dum = bsxfun(@minus,dum,stds);
%             dum = dum>0;
%             one_norms(i:i+fin) = sum(dum,1);
            %old
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
        f = f*opts.mem_occup_ratio_CPU;
        chunk_size = floor(f/ (2*N+num_pcs)); %# of pixels at once

        one_norms = zeros(1,N);
        for i = 1:chunk_size:N
            fin = min(chunk_size-1,N-i);
            Q = U_norm(i:i+fin,:)';
            %new
%             dum = U*Q;
%             stds = 3*std(dum,0,1);
%             dum = bsxfun(@minus,dum,stds);
%             dum = dum>0;
%             one_norms(i:i+fin) = sum(dum,1);
            %old
            one_norms(i:i+fin) = (sum(abs(U*Q),1));
            clear Q;
        end
    end

    dispfun(sprintf('\t\t\t Beginning source separation..\n'),verb2);
    idx_possible = 1:N;
    goodness_sigs = zeros(1,max_num);
    Q = zeros(num_pcs,max_num,'single');
    S = sparse(N,max_num);
    
    acc=0;
    w_count = 0;
    num_sweep = 25;
    idx_stop = ceil(num_pcs*linspace(0.2,1,num_sweep));
    num_reject_overlap = 0;
    num_reject_small = 0;
    
    % Statistics for debugging 
    % 1:is_good , 2:one_norm , 3-4:idx , 5-6:max_loc , 7-8:centroid
    % 9:cell_size , 10:#pixels left to search 11:neighbor_one_norms
    debug_stats = zeros(max_num,11);
    
    while true        
        w_count = w_count+1;
        
        [~,idx_this] = min(one_norms(idx_possible));
        idx_this = idx_possible(idx_this);
        
        [sub_this_y,sub_this_x] = ind2sub([h,w],idx_kept(idx_this));
        debug_stats(w_count,[3,4]) = [sub_this_x,sub_this_y];        
        
        % Find the best #pcs for the chosen signal index
        q = U_norm(idx_this,:)';
        q_sweep = repmat(q,1,num_sweep);
        for i = 1:num_sweep
            q_sweep(idx_stop(i)+1:end,i) = 0;
        end
        norms_q_sweep = sqrt(sum(q_sweep.^2,1));
        q_sweep = bsxfun(@times,q_sweep,1./norms_q_sweep);
        sigs_sweep = U*q_sweep;
        %new
%         stds = 3*std(sigs_sweep,0,1);
%         sigs_sweep = bsxfun(@minus,sigs_sweep,stds);
%         sigs_sweep = sigs_sweep>0;
%         goodness_sweep = sum(sigs_sweep,1);

        %old
        goodness_sweep = sum(abs(sigs_sweep),1);
        
        [goodness_sig,best_idx] = min(goodness_sweep);
        
        debug_stats(w_count,2) = goodness_sig/sqrt(N);
        
        sig_this = U*q_sweep(:,best_idx); 
        std_sig = double(std(sig_this));
        
        [~,idx] = max(sig_this);
        [sub_max_y,sub_max_x] = ind2sub([h,w],idx_kept(idx));
        debug_stats(w_count,[5,6]) = [sub_max_x,sub_max_y];
        
        % Calculate the signal goodness if neighbors of max-location were
        % selected
        idx_neighbors = [];
        for x = -1:1
            for y = -1:1
                if  (x~=0) || (y~=0)
                    neighbor = [sub_max_y+y,sub_max_x+x];
                    if ~any(neighbor<=0) && ~any(neighbor>h) && ~any(neighbor>w) %check if in range
                        idx = sub2ind([h,w],neighbor(1),neighbor(2));
                        idx_trimmed = find(idx_kept==idx);
                        if idx_trimmed>0
                            idx_neighbors(end+1) = idx_trimmed;
                        end
                    end
                end
            end
        end
       if ~isempty(idx_neighbors)
            q_neighbor = U(idx_neighbors,1:idx_stop(best_idx))';
            norm_q_neighbor = sqrt(sum(q_neighbor.^2,1));
            q_neighbor = bsxfun(@times,q_neighbor,1./norm_q_neighbor);
            dum = U(:,1:idx_stop(best_idx))*q_neighbor;
            one_norms_neighbors = sum(abs(dum),1)/sqrt(N);
            debug_stats(w_count,11) = max(one_norms_neighbors);
            debug_stats(w_count,12) = mean(one_norms_neighbors);
       end
        
        % Check if it overlaps significantly with any other signal
        sig_temp = zeros(h*w,1);
        sig_temp(idx_kept) = sig_this;
        [sig_temp,cent] = cleanup_source(sig_temp,0.2);
        pixels_this = find(sig_temp);%sig_this>3*std_sig
        
        debug_stats(w_count,[7,8]) = [cent(1),cent(2)];
        debug_stats(w_count,9) = length(pixels_this);
        
        if ~isempty(pixels_this) && length(pixels_this)< h*w/100
            debug_stats(w_count,1) = 1;
            sig_temp = sig_temp(idx_kept);
            pixels_this = find(sig_temp);%sig_this>3*std_sig
            sig_temp = sparse(double(sig_temp));
            sig_temp = sig_temp/norm(sig_temp);
            overlaps = S'*sig_temp;

            % Add if new signal
            if ~any(overlaps>0.65)
                debug_stats(w_count,1) = 2;
                acc = acc+1;
                S(sig_temp>0,acc) = sig_temp(sig_temp>0);
                Q(:,acc) = q_sweep(:,best_idx);        
                goodness_sigs(acc) = goodness_sig;           
                idx_possible = setdiff(idx_possible,pixels_this);
                
                debug_stats(w_count,10) = length(idx_possible);
            else
                idx_possible = setdiff(idx_possible,pixels_this);
                num_reject_overlap = num_reject_overlap+1;
            end
            
        elseif length(pixels_this)>= h*w/100
            debug_stats(w_count,1) = 10;           
        end
        
        idx_possible = setdiff(idx_possible,idx_this);

        if isempty(idx_possible) || w_count==max_num
            break;
        end
        fprintf('iter:%d, num_sig: %d \n',w_count,acc);
    end
    % Truncate in case there are less components than max_num
    S = S(:,1:acc);
    Q = Q(:,1:acc);
    norms_Q = sqrt(sum(Q.^2,1));
    Q = bsxfun(@times,Q,1./norms_Q);
    goodness_sigs = goodness_sigs(1:acc);
    debug_stats = debug_stats(1:w_count,:);
    %sort debug_stats by their inv_qualities
    [~,idx] = sort(debug_stats(:,2),'ascend');
    debug_stats = debug_stats(idx,:);
    save('debug_stats','debug_stats');
    dispfun(sprintf('\t\t\t Done. Extracted %d potential cells in total. \n',acc),verb2);
end


function [F_out,cent_out,idx_retained] = cleanup_sources(F_in,mag_thresh)
    
    size_thresh = opts.ss_cell_size_threshold;
    num_cells = size(F_in,2);
    
   % Initialize output variables
    cent_out = zeros(2,num_cells);
    F_out = zeros(h,w,num_cells,'single'); 
    
%     % Remove baselines
%     meds = median(F_in,1);
%     F_in = bsxfun(@minus,F_in,meds);

    % Reshape into 2D
    F_in = reshape(F_in,h,w,num_cells);
 
    acc = 0;
    idx_elim = [];
    for idx_cell = 1:num_cells
        this_cell = F_in(:,:,idx_cell);
        [mx,idx_mx] = max(this_cell(:));
        this_cell(this_cell<mx*mag_thresh) = 0; % Threshold
        this_mask = this_cell>0;
        CC = bwconncomp(this_mask);
        
        % Find the connected component that has idx_mx
        lens = cellfun(@length, CC.PixelIdxList);
        [~,idx_sort] = sort(lens,'descend');
        CC.PixelIdxList = CC.PixelIdxList(idx_sort); 
        
        acc2 = 0;
        while 1
            acc2 = acc2+1;
            if ~isempty(find(CC.PixelIdxList{acc2}==idx_mx,1))
                break;
            end
        end
        this_mask = zeros(h,w,'single');
        this_mask(CC.PixelIdxList{acc2})=1;
        if nnz(this_mask) >= size_thresh
            acc = acc+1;
            F_out(:,:,acc) = this_mask.*this_cell;
            s = regionprops(this_mask,'centroid');
            cent_out(:,acc) = s(1).Centroid;
        else
            idx_elim = [idx_elim,idx_cell];
        end        
    end
        
    % Truncate output variables
    F_out = F_out(:,:,1:acc);
    cent_out = cent_out(:,1:acc);
    
    % Reshape back into 1D
    F_out = reshape(F_out,h*w,acc);
    idx_retained = setdiff(1:num_cells,idx_elim);
end

function [F_out,cent] = cleanup_source(F_in,mag_thresh)

    size_thresh = opts.ss_cell_size_threshold;
   % Initialize output variable
    F_out = zeros(h,w,'single'); 

    % Reshape into 2D
    F_in = reshape(F_in,h,w);

    [mx,idx_mx] = max(F_in(:));
    F_in(F_in<mx*mag_thresh) = 0; % Threshold
    this_mask = F_in>0;
    CC = bwconncomp(this_mask);

    % Find the connected component that has idx_mx
    lens = cellfun(@length, CC.PixelIdxList);
    [~,idx_sort] = sort(lens,'descend');
    CC.PixelIdxList = CC.PixelIdxList(idx_sort); 

    acc2 = 0;
    while 1
        acc2 = acc2+1;
        if ~isempty(find(CC.PixelIdxList{acc2}==idx_mx,1))
            break;
        end
    end
    this_mask = zeros(h,w,'single');
    this_mask(CC.PixelIdxList{acc2})=1;

    F_out = this_mask.*F_in;
    cent_x = sum(F_in*(1:w)')/sum(F_in(:));
    cent_y = sum((1:h)*F_in)/sum(F_in(:));
    cent = [cent_x,cent_y];
    % Reshape back into 1D
    if sum(F_out(:)>mag_thresh*mx) >  size_thresh
        
        F_out = reshape(F_out,h*w,1);
    else
        F_out = [];
    end
    
end


function dispfun(str,state)
    if state == 1
        fprintf(str);
    end
end

end