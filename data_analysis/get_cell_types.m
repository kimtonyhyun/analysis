function p_mat = get_cell_types(ds,visualize)
%
% Finds cell types by testing population mean difference for trials in and out of
% each defined category.
% There are 14 categories (se sw en es c i se-c se-i sw-c sw-i en-c en-i
% es-c es-i).
%
% Inputs:
%   ds: DaySummary object
%   visualize (optional): Set to 1 for plotting a detailed visualization
%   for each cell that passes the test at level 1e-9 for at least one
%   category.
%
% Output:
%   p_mat : [# of classified cells]x[# of categories] array containing
%   p-values for each cell and for each category.
%

    idx_cells = find(ds.is_cell);
    num_cells = length(idx_cells);
    num_trials = ds.num_trials;
    num_categ = 14; 
    feature_mat = zeros(num_trials,num_cells);
    categ_mat = zeros(num_trials,num_categ);

    % Construct the two input matrices
    idx_probes=[];
    for idx_trial = 1:num_trials
        
        trial = ds.trials(idx_trial);
        feature_mat(idx_trial,:) = max(trial.traces(idx_cells,:),[],2)';

        val_start = strcmp(trial.start,'west')-strcmp(trial.start,'east');
        if val_start==0 % probe trial, skip
            idx_probes(end+1) = idx_trial;
            categ_mat(idx_trial,:) = -ones(1,num_categ);
            continue
        end
        val_finish = strcmp(trial.end,'south')-strcmp(trial.end,'north');
        val_correct = trial.correct;
        categ_mat(idx_trial,:) = compute_categories(val_start,val_finish,val_correct);
        
    end
    
    % Decide cell type
    p_mat = zeros(num_cells,num_categ);
    for idx_cell = 1:num_cells
        means_cell = feature_mat(:,idx_cell);
        for idx_categ = 1:num_categ
            categ = categ_mat(:,idx_categ);
            [~,p] = ttest2(means_cell(categ==1),means_cell(categ==0),'Tail','right','Vartype','equal' );
            p_mat(idx_cell,idx_categ) = p;
        end
    end
    
    if exist('visualize','var')
        if visualize
            visualize_cell_types;
        end
    end

    function categ = compute_categories(val_start,val_finish,val_correct)
        categ = zeros(1,14);
        if val_start==-1 % east-start
            categ(1) = 1;
            if val_correct
                categ(7) = 1;
            else
                categ(8) = 1;
            end
        elseif val_start==1 % west-start
            categ(2) = 1;
            if val_correct
                categ(9) = 1;
            else
                categ(10) = 1;
            end
        end

        if val_finish==-1 % north-finish
            categ(3) = 1;
            if val_correct
                categ(11) = 1;
            else
                categ(12) = 1;
            end
        elseif val_finish==1 % south-finish
            categ(4) = 1;
            if val_correct
                categ(13) = 1;
            else
                categ(14) = 1;
            end
        end
        
        if val_correct %correct
            categ(5)=1;
        else
            categ(6)=1;
        end
    end

    function visualize_cell_types
        categ_strs = {'se','sw','en','es','c','i','se-c','se-i','sw-c','sw-i','en-c','en-i','es-c','es-i'};
        [ps,idx_categs] = min(p_mat,[],2);

        alph = 1e-9; % Significance threshold (test level)

        for i = 1:length(ps)
            if ps(i) < alph

                categ_str = categ_strs{idx_categs(i)};
                is_in_categ = categ_mat(:,idx_categs(i))==1;
                is_out_categ = categ_mat(:,idx_categs(i))==0;
                features_in = feature_mat(is_in_categ,i);
                features_out = feature_mat(is_out_categ,i);

                subplot(3,3,[1,4,7])
                whole_raster = ds.plot_cell_raster(idx_cells(i));
                minval = min(min(whole_raster));
                maxval = max(max(whole_raster));
                str = sprintf('cell %d, cell type = %s \n p=%d',idx_cells(i),categ_str,ps(i));
                title(str);
                subplot(3,3,[2,3])
                plot(features_in);
                hold on
                plot(features_out);
                hold off
                xlabel('Trial index')
                title(sprintf('Max values of traces in each trial \n Trials in and outside of the category are separated'))
                legend('in the category','not in the category')
                subplot(3,3,[5,8]);
                x_axis = linspace(0, 1, 1000);
                imagesc(x_axis,1:sum(is_in_categ), whole_raster(is_in_categ,:),[minval,maxval]);
                xlabel('Trial phase [a.u.]');
                ylabel('Trial index');
                title(sprintf('Trials in the category  \"%s\" ',categ_str));
                subplot(3,3,[6,9]);
                imagesc(x_axis,1:sum(is_out_categ), whole_raster(is_out_categ,:),[minval,maxval]);
                xlabel('Trial phase [a.u.]');
                ylabel('Trial index');
                title(sprintf('Trials not in the category  \"%s\" ',categ_str));
                pause;
            end
        end
    
    end
end