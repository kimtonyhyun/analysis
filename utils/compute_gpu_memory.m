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