function browse_corrlist(corrlist, ds1, ds2, varargin)

% Default settings
app_name = 'Browse corrlist';
num_pairs = size(corrlist, 1);

% Interactive loop
%------------------------------------------------------------
hf = figure;

idx = 1; 
while (1)
    corrdata = corrlist(idx, :);
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
                case 's' % Save image
                    print(sprintf('cell%d.png', idx), '-dpng');

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

