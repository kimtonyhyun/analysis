function [F,cents,inv_qualities] = compute_spatial_weights(svd_source,opts)
%   Cell extraction method that is based on an algorithm that essentially
%   solves the sparse source separation problem to find cell spatial weights.

% Some parameters
corr_thresh = 0.3;
mag_thresh = 0.1;

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

function [Q_fine,inv_qualities] = extract_sources(U,max_num)
% Extract sources recursively using the 1-norm criterion

    [N,num_pcs] = size(U);
    norms_U = sqrt(sum(U.^2,2));
    U_norm = bsxfun(@times,U,1./norms_U);
    
    % Brute-force approach for 1-norm computation, use GPU if available
    gpu_exists = gpuDeviceCount>0;
    if gpu_exists && (~opts.disableGPU)
        avail_mem = compute_gpu_memory();    
        f = max(avail_mem/4 - N*(num_pcs+2),0); % maximum available element size 
        f = f*opts.mem_occup_ratio_GPU;
        
        chunk_size = floor(f/ (2*N+num_pcs)); %# of pixels at once
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
            one_norms(i:i+fin) = sum(abs(U*Q),1);
            clear Q;
        end
    end

    dispfun(sprintf('\t\t\t Beginning source separation..\n'),verb2);
    idx_possible = 1:N;
    idx_cells = zeros(1,max_num);
    acc=0;

    while true
        acc = acc+1;
        [~,idx_this] = min(one_norms(idx_possible));
        idx_this = idx_possible(idx_this);
        sig_this = U*U(idx_this,:)'/norms_U(idx_this);
        correl = sig_this.*(1./norms_U);
        idx_possible = intersect(find(abs(correl)<corr_thresh),idx_possible);
        idx_cells(acc) = idx_this;
        % Termination condition
        if isempty(idx_possible) || acc==max_num
            break;
        end
    end

    % Truncate in case there are less components than max_num
    idx_cells = idx_cells(1:acc);

    % Q is the coarse linear mapping from U to weights (UQ=F)
    Q = U_norm(idx_cells,:)';
    
    dispfun(sprintf('\t\t\t Found %d sources in the initial pass.\n',acc),verb2);
    
    dispfun(sprintf('\t\t\t Fine tuning the sources...\n'),verb2);
    % Fine tuning: estimate the optimum # of singular vectors to use for each source
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
    end
    
    % Reconstruct Q (Q_fine)
    [best_one_norms,idx_best_rat] = min(one_norms_F,[],2);
    idx_uncertain = find(idx_best_rat==1);
%     idx_best_rat(idx_best_rat==1) = num_sweep/2;

    % Deal with indices where its not certain which ratio is the best 
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

    dispfun(sprintf('\t\t\t Doing a 2nd pass for removing duplicate sources...\n'),verb2);
    
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
        correl = Corr_cells(idx_this,:);
        idx_possible = intersect(find(abs(correl)<corr_thresh+0.1),idx_possible);
        idx_cells(acc) = idx_this;
        % Termination condition
        if isempty(idx_possible) || acc==max_num
            break;
        end
    end
    idx_cells = idx_cells(1:acc);
    Q_fine = Q_fine(:,idx_cells);
    inv_qualities = best_one_norms(idx_cells)'/sqrt(N);
    size(inv_qualities,1);
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

    % Threshold sources + deal with multiple boundaries    
    acc = 0;
    idx_elim = [];
    for idx_cell = 1:num_cells
        this_cell = F_in(:,:,idx_cell);
        [boundaries, ~] = compute_ic_boundary(this_cell, mag_thresh);
        num_objects = length(boundaries);
        if num_objects <100
            for idx_obj = 1:min(1,num_objects)
                mask_candid = poly2mask(boundaries{idx_obj}(:,1), boundaries{idx_obj}(:,2), h, w);
                if nnz(mask_candid)>=size_thresh
                    s = regionprops(mask_candid,'centroid');
                    if ~isempty(s)
                        acc = acc+1;
                        F_out(:,:,acc) = mask_candid.*this_cell;
                        cent_out(:,acc) = s(1).Centroid;
                    else
                        idx_elim = [idx_elim,idx_cell];
                    end
                else
                    idx_elim = [idx_elim,idx_cell];
                end
            end
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

function dispfun(str,state)
    if state == 1
        fprintf(str);
    end
end

end