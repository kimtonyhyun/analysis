function avail_mem = compute_cpu_memory()
% Available memory on CPU in bytes

    if ispc % Windows
         [~,sys] = memory;
         avail_mem = (sys.PhysicalMemory.Available);
    elseif ismac %Mac
         [~,m]=unix('vm_stat | grep free');
         spaces=strfind(m,' ');
         avail_mem = str2double(m(spaces(end):end))*4096;
    else % Linux
         [~,meminfo] = system('cat /proc/meminfo');
         tokens = regexpi(meminfo,'^*MemAvailable:\s*(\d+)\s','tokens');
         avail_mem = str2double(tokens{1}{1})*1024;
    end
end