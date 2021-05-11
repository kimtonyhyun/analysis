function matched_inds = browse_corrlist(corrlist, ds1, ds2, varargin)

% Default settings
app_name = 'Browse corrlist';
num_pairs = size(corrlist, 1);
matched_inds = size(num_pairs,2);

% Interactive loop
%------------------------------------------------------------
hf = figure;

idx = 1;
num_matched = 0;
while (1)
    corrdata = corrlist(idx, :); % [cell1_idx cell2_idx corr]
    show_corr(ds1, corrdata(1), ds2, corrdata(2), corrdata(3), varargin{:});
    
    prompt = sprintf('%s (%d of %d) >> ', app_name, idx, num_pairs);
    resp = strtrim(input(prompt, 's'));
    
    val = str2double(resp);
    if (~isnan(val)) % Is a number
        if (1 <= val) && (val <= num_pairs)
            idx = val;
        end
    else
        resp = lower(resp);
        if isempty(resp)
            idx = idx + 1;
            idx = min(num_pairs, idx);
        else
            switch resp(1)
                case 'c' % Save to inds (used in 1P:2P)
                    num_matched = num_matched + 1;
                    matched_inds(num_matched,1) = corrdata(1);
                    matched_inds(num_matched,2) = corrdata(2);
                    fprintf('  Saved to matched_inds\n');
                    
                    idx = idx + 1;
                    idx = min(num_pairs, idx);
                    
                case 's' % Save image
                    image_filename = sprintf('cell%d.png', idx);
                    print(image_filename, '-dpng');
                    fprintf('  Imaged saved to "%s"\n', image_filename);

                case 'p' % Previous
                    idx = idx - 1;
                    idx = max(1, idx);

                case 'q' % Exit
                    close(hf);
                    break;

                otherwise
                    fprintf('  Could not parse "%s"\n', resp);
            end
        end
    end
end % while (1)

matched_inds = matched_inds(1:num_matched,:);