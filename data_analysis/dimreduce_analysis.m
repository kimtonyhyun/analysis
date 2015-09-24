function dimreduce_analysis(md,method,training_set,test_set,opts)
% Do population analysis
%
% Inputs: 
%
%   md: multiday object
%
%   method: 'pca' or 'pls'
%
%   training_set: struct with following fields: 'day', 'attr','label',
%   'disp_label'. 'attr' is a cell array specifying the desired trial types in
%   a day. An example 'attr' is {'start','east','correct'}. This only
%   retains the trials that start in the east arm and are correct. 'Label'
%   can be (in theory) any real numeric value, but recommended type is
%   integer. 'label' is not mandatory for method='pca'. 'displabel' is an
%   optional field that specifies how the element is desired to be referred
%   to in the plots. If left blank, it's set to be the aggregation of
%   'attr' elements and the corresponding day.
%
%   test_set: Optional input for testing under the model constructed with 
%   training_set. Type is same as training_set. 'label' field is optional for
%   test_set.
%
%   opts: Optional struct input whose fields overwrites defaults set in the
%   script. These are:
%       'trial_phase': Default is 1(pre-run). Can be a 1D array composed
%       of numbers 1(pre-run), 2(run), and 3(post-run).
%       'ratio_PCs': PCA is computed prior to PLS to reduce noise. This
%       parameter controls what ratio of PCs is retained for the PLS step.
%       Default is 0.5.
%       pca_dims: 2x1 array with PC dimensions to use. Default is [1,2]
%       pls_dims: same as above for pls
%
% Output: Produces a plot.
%

if ~exist('opts','var')
    opts = [];
end

% Defaults
do_test = 0;
if exist('test_set','var') && ~isempty(test_set)
    do_test = 1;
end

if ~isfield(opts,'trial_phase'),  opts.trial_phase=1; end
if ~isfield(opts,'ratio_PCs'),  opts.ratio_PCs=0.5; end

num_cells = md.num_cells;

days_training = cell2mat({training_set.day});
num_training = length(days_training);
attr_training = {training_set.attr};

if do_test
    days_test = cell2mat({test_set.day});
    num_test = length(days_test);
    attr_test = {test_set.attr};

    days = [days_training,days_test];
    num_days = num_training+num_test;
    attr = [attr_training,attr_test];
else
    days = days_training;
    num_days = num_training;
    attr = attr_training;
end

%----------------------------
% Construct features and labels
%----------------------------
data = [];
num_trials_each_day = [];
for i = 1:num_days
    day = days(i);
    
    % Check attributes and retrieve trials
    start_arm = '';
    end_arm='';
    turn_dir='';
    start_spec = find(strcmp(attr{i},'start'));
    if start_spec>0
        start_arm = attr{i}{start_spec+1};
    end
    end_spec = find(strcmp(attr{i},'end'));
    if end_spec>0
        end_arm = attr{i}{end_spec+1};
    end
    turn_spec = find(strcmp(attr{i},'end'));
    if turn_spec>0
        turn_dir = attr{i}{turn_spec+1};
    end
    tr = '';
    trials_spec = find(strcmp(attr{i},'trials'));
    if trials_spec>0
        tr = attr{i}{trials_spec+1};
    end
    if any(strcmp(attr{i},'correct'))
        chunk_day = md.get_trials(day,'trials',tr,'start',start_arm,'end',end_arm,...
            'turn',turn_dir,'correct','split_trial_phases');
    elseif any(strcmp(attr{i},'incorrect'))
        chunk_day = md.get_trials(day,'trials',tr,'start',start_arm,'end',end_arm,...
            'turn',turn_dir,'incorrect','split_trial_phases');
    else
        chunk_day = md.get_trials(day,'trials',tr,'start',start_arm,'end',end_arm,...
            'turn',turn_dir,'split_trial_phases');
    end
    
    % Convert fluorescence into features
    num_trials = length(chunk_day);
    for j = 1:num_trials % Add 1 sample at a time to the feature matrix
        chunk_trial = chunk_day(j).traces; % 1x3 cell 
        data_trial = zeros(num_cells,1);
        for trial_phase = opts.trial_phase
            data_trial = data_trial+mean(chunk_trial{trial_phase},2);
        end
        data_trial = data_trial / length(opts.trial_phase); % Average
        data = [data;data_trial'];
    end
    num_trials_each_day = [num_trials_each_day,num_trials];
end

% Split data into 2 training and test
idx_end_training = sum(num_trials_each_day(1:num_training));
data_training = data(1:idx_end_training,:);
if do_test
    data_test = data(idx_end_training+1:end,:);
end

% Class and display labels
if strcmp(method,'pls')    
    labels_training_days = cell2mat({training_set.label});
    labels_training = [];
    for i = 1:num_training
        labels_training = [labels_training;...
            repmat(labels_training_days(i),num_trials_each_day(i),1)];
    end
end

displabels_training = {};
for i = 1:num_training
    if isfield(training_set,'displabel')
        displabel_str = training_set(i).displabel;
    else
        displabel_str = num2str(days_training(i));
        for j = 1:length(attr{i})
            this_attr = attr_training{i}{j};
            if isnumeric(this_attr) % Trials:indicate begin and end
                str_start = num2str(this_attr(1));
                str_end = num2str(this_attr(end));
                displabel_str = [displabel_str,'-',str_start,'-',str_end];
            else
                displabel_str = [displabel_str,'-',this_attr];
            end            
        end
    end
    for j = 1:num_trials_each_day(i)
        displabels_training = [displabels_training;displabel_str];
    end
end

if do_test
    if isfield(test_set,'label')
        labels_test_days = cell2mat({test_set.label});
        labels_test = [];
        for i = 1:num_test
            labels_test = [labels_test;...
                repmat(labels_test_days(i),num_trials_each_day(i+num_training),1)];
        end
    end

    displabels_test = {};
    for i = 1:num_test
        if isfield(training_set,'displabel')
            displabel_str = test_set(i).displabel;
        else
            displabel_str = num2str(days_test(i));
            for j = 1:length(attr_test{i})
                this_attr = attr_test{i}{j};
                if isnumeric(this_attr) % Trials:indicate begin and end
                    str_start = num2str(this_attr(1));
                    str_end = num2str(this_attr(end));
                    displabel_str = [displabel_str,'-',str_start,'-',str_end];
                else
                    displabel_str = [displabel_str,'-',this_attr];
                end
            end
        end
        for j = 1:num_trials_each_day(i+num_training)
            displabels_test = [displabels_test;displabel_str];
        end
    end
end

%-----------------
% Dimension reduction
%-----------------

mean_data = mean(data_training,1);
[coeff,score,~,~,explained,~] = pca(data_training);

if strcmp(method,'pca') % PCA
    if isfield(opts,'pca_dims')
        dim1 = opts.pca_dims(1);
        dim2 = opts.pca_dims(2);
    else
        dim1=1;dim2=2;
    end
    
    data_1 = score(:,dim1);
    data_2 = score(:,dim2);
    displabels = displabels_training;
    markers = repmat('.',1,num_training);
    title_str = 'PCA : Projection onto 2 dimensions';
    
    if do_test
        centered_test_data = bsxfun(@minus,data_test,mean_data);
        dr_filter = coeff(:,[dim1,dim2]);
        dum =  centered_test_data * dr_filter ;
        test_points = dum(:,[dim1,dim2]);
        
        data_1  = [data_1;test_points(:,1)];
        data_2  = [data_2;test_points(:,2)];
        displabels = [displabels;displabels_test];
        markers = [markers,repmat('o',1,num_test)];
        title_str = [title_str,' - Filled: Training, Hollow: Test'];
    end
    gscatter(data_1,data_2,displabels,[],markers);
    
    xlabel(['PC Dimension ',num2str(dim1)],'FontSize',15)
    ylabel(['PC Dimension ',num2str(dim2)],'FontSize',15)
    prcnt_explained = sum(explained([dim1,dim2]));
    title(sprintf('%s\n Percent explained: %.2f',title_str,prcnt_explained),'FontSize',18);
    if isfield(opts,'pca_xylim'),
        xlim(opts.pca_xylim(1,:));
        ylim(opts.pca_xylim(2,:));
    end
    
else % PLS
    num_PCs =size(score,2)*opts.ratio_PCs;
    % Partial Least squares
    data_reduced = score(:,1:num_PCs);
    [XL,~,XS] = plsregress(data_reduced,labels_training);
    if isfield(opts,'pls_dims')
        dim1 = opts.pca_dims(1);
        dim2 = opts.pca_dims(2);
    else
        dim1=1;dim2=2;
    end
    
    % Fit 2D glm (suppress warnings when perfectly separated)
    warning('off', 'stats:glmfit:IterationLimit');
    warning('off', 'stats:glmfit:PerfectSeparation');
    w = glmfit(XS(:,[dim1,dim2]),labels_training,'binomial');
    l_tr = mean((XS(:,[dim1,dim2])*w(2:3)+w(1) > 0)==labels_training);   
    
    data_1 = XS(:,dim1);
    data_2 = XS(:,dim2);
    displabels = displabels_training;
    markers = repmat('.',1,num_training);
    accuracy_str = sprintf('Training Accuracy: %d percent',round(100*l_tr));
    title_str = 'PLS : Projection onto 2 dimensions';
    
    if do_test
        centered_test_data = bsxfun(@minus,data_test,mean_data);
        dr_filter = ( XL \coeff(:,1:num_PCs)')';
        dum =  centered_test_data * dr_filter ;
        test_points = dum(:,[dim1,dim2]);  
        
        data_1  = [data_1;test_points(:,1)];
        data_2  = [data_2;test_points(:,2)];
        displabels = [displabels;displabels_test];
        markers = [markers,repmat('o',1,num_test)];
        title_str = [title_str,' - Filled: Training, Hollow: Test'];
        if isfield(test_set,'label')
            l_te = mean((test_points*w(2:3)+w(1) > 0)==labels_test);
            accuracy_str = {accuracy_str,sprintf('Test Accuracy: %d percent',...
                round(100*l_te))};
        end
    else
        accuracy_str = {accuracy_str}; % Still make a cell array
    end
    
    gscatter(data_1,data_2,displabels,[],markers);
    
    hold on;
    %Plot Decision Boundary
    dx = linspace(min(XS(:,dim1)),max(XS(:,dim1)),5);
    dy = (-dx*w(2)-w(1))/w(3);
    plot(dx,dy,'--k','LineWidth',2); 

    if isfield(opts,'pls_xylim'),
        xlim(opts.pls_xylim(1,:));
        ylim(opts.pls_xylim(2,:));
    end
    xlabel(['PLS Dimension ',num2str(dim1)],'FontSize',15)
    ylabel(['PLS Dimension ',num2str(dim2)],'FontSize',15)
    title(title_str,'FontSize',18);
    
    textloc = [min(dx),min(dy)];
    text(textloc(1),textloc(2),accuracy_str{1},'Fontsize',15,'Color','r');
    if do_test
        text(textloc(1),textloc(2)-0.05,accuracy_str{2},'Fontsize',15,'Color',[0,0.5,0]);
    end
    legend('boxoff')
    hold off
end
