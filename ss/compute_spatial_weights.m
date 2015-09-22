function [F,cents,inv_qualities] = compute_spatial_weights(svd_source)

%   Cell extraction method that is based on an algorithm that essentially
%   solves the sparse source separation problem to find cell spatial weights.

% Some parameters
mem_occup_scale_CPU = 0.5; % Occupy this fraction of available memory on RAM
mem_occup_scale_GPU = 0.7; % Occupy this fraction of available memory on GPU
corr_thresh = 0.2;
mag_thresh = 0.25;
max_num = 3000;

if ischar(svd_source)
    str = svd_source;
    svd_source = get_most_recent_file('','svd_*.mat');
    svd = load(svd_source,['SVD',str]);
    U = eval(['svd.SVD',str,'.U']);
    info = eval(['svd.SVD',str,'.info']);
else
    U = svd_source.U;
    info = svd_source.info;
end
% U = U(:,1:200);
is_trimmed = info.trim.enabled;
h = info.movie_height;
w = info.movie_width;
idx_kept = info.trim.idx_kept;

fprintf('%s: Extracting cell spatial weights...\n',datestr(now));

[Q_fine,inv_qualities] = extract_sources(U,max_num);
num_extracted = size(Q_fine,2);
F = U*Q_fine;

if is_trimmed % Untrim the pixels
    F_temp = zeros(h*w,num_extracted);
    F_temp(idx_kept,:) = F;
    F = F_temp;
    clear F_temp;
end

fprintf('%s: Checking cells and doing some cleaning...\n',datestr(now));
[F,cents,idx_retained] = cleanup_sources(F,mag_thresh);
inv_qualities = inv_qualities(idx_retained);

F = bsxfun(@times,F,1./sum(F,1)); % Make each weight vector sum up to 1

fprintf('%s: Done. \n',datestr(now));

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

function [Q_fine,inv_qualities] = extract_sources(U,max_num)
% Extract sources recursively using the 1-norm criterion

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
        fprintf('\t\t\t GPU detected. Using it for computation...\n');
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

    fprintf('\t\t\t Beginning source separation..\n');
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
    
    fprintf('\t\t\t Found %d sources in the initial pass.\n',acc);
    
    fprintf('\t\t\t Fine tuning the sources...\n');
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

    fprintf('\t\t\t Doing a 2nd pass for removing duplicate sources...\n');
    
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
    fprintf('\t\t\t Done. Extracted %d potential cells in total. \n',acc);
end


function [F_out,cent_out,idx_retained] = cleanup_sources(F_in,mag_thresh)
    
    size_thresh = 40;
    num_cells = size(F_in,2);
    
    % Initialize output variables
    cent_out = zeros(2,num_cells);
    F_out = zeros(h,w,num_cells*2); % Be generous in size
    
    % Remove baselines
    meds = median(F_in,1);
    F_in = bsxfun(@minus,F_in,meds);

    % Reshape into 2D
    F_in = reshape(F_in,h,w,num_cells);

    % Threshold sources + deal with multiple boundaries    
    acc = 0;
    idx_elim = [];
    for idx_cell = 1:num_cells
        this_cell = F_in(:,:,idx_cell);
        [boundaries, ~] = compute_ic_boundary(this_cell, mag_thresh);
        num_objects = length(boundaries);
        if num_objects <20
            for idx_obj = 1:min(1,num_objects)
                mask_candid = poly2mask(boundaries{idx_obj}(:,1), boundaries{idx_obj}(:,2), h, w);
                if nnz(mask_candid)>=size_thresh
                    acc = acc+1;
                    s = regionprops(mask_candid,'centroid');
                    F_out(:,:,acc) = mask_candid.*this_cell;
                    cent_out(:,acc) = s(1).Centroid;
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
    
    fprintf('\t\t\t %d Cells were retained after cleanup.\n',acc);
end

end