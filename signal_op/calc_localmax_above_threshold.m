function localmax = calc_localmax_above_threshold(threshold,data)

% Find the local maxima in a sequence of points provided that the maxima are 
%above a specified threshold

% Make data a row vector
if size(data,2)==1, data = data'; end

thresh_dist = 5; % minimum distance required between events
data_binary = data>threshold; % Binarize

localmax = [];

k=1;
while k<length(data_binary)
    if data_binary(k)==1 % Start a sequence of 1's
        
        begin_seq = k; 
        while (data_binary(k)==1)&&(k<length(data_binary)), k=k+1; end
        end_seq = k-1;
        
        % After the end of sequence calculate local maxima in it
        
        seq = data(begin_seq:end_seq); % Magnitudes of the events in the sequence
        idx_localmax = find(conv(seq,[-1,2,-1],'valid'))+1; % Indices of local maxima
        
        for i = idx_localmax % Check the quality of each local maxima
            dists = abs(find(seq>seq(i))-i);
            dists = setdiff(dists,0); % Exclude self
            min_dist = min(dists);
            
            if isempty(dists)
                localmax = [localmax,i+begin_seq-1]; %#ok
            elseif min_dist>thresh_dist % Impose min_dist requirement
                localmax = [localmax,i+begin_seq-1]; %#ok
            end
            
        end

    else %skip through 0's
        k = k+1;
    end
end