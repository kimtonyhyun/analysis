function [good_fits, okay_fits, bad_fits] = group_fits(score)

num_cells = length(score.vals_n);

switch score.name
    case 'Rsq'
        bad_threshold = 0.4;
        okay_threshold = 0.7;
        
        vals = score.vals_n;
        
        bad_fits = find(vals < bad_threshold);
        
        okay_fits = find(vals < okay_threshold);
        okay_fits = setdiff(okay_fits, bad_fits);
        
        good_fits = setdiff((1:num_cells)', [okay_fits; bad_fits]);
        
    otherwise
        error('Score name "%s" is not supported!', score.name);
end