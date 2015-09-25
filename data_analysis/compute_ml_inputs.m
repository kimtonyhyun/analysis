function data = compute_ml_inputs(md,training_set,test_set,opts)
% Computes features and labels from neuron traces to be used in subsequent 
% analysis methods.
% Meant to be called from data analysis functions, refer to one of them for
% documentation on the inputs.
%

% Defaults
do_test = 0;
if exist('test_set','var') && ~isempty(test_set)
    do_test = 1;
end

if ~isfield(opts,'trial_phase'),  opts.trial_phase=1; end

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
features = [];
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
        features_trial = zeros(num_cells,1);
        for trial_phase = opts.trial_phase
            features_trial = features_trial+mean(chunk_trial{trial_phase},2);
        end
        features_trial = features_trial / length(opts.trial_phase); % Average
        features = [features;features_trial'];
    end
    num_trials_each_day = [num_trials_each_day,num_trials];
end

% Split data into training and test
idx_end_training = sum(num_trials_each_day(1:num_training));
data.features_training = features(1:idx_end_training,:);
if do_test
    data.features_test = features(idx_end_training+1:end,:);
end

% Class and display labels
if isfield(training_set,'label')    
    labels_training_days = cell2mat({training_set.label});
    labels_training = [];
    for i = 1:num_training
        labels_training = [labels_training;...
            repmat(labels_training_days(i),num_trials_each_day(i),1)];
    end
    data.labels_training = labels_training;
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
data.displabels_training = displabels_training;

if do_test
    if isfield(test_set,'label')
        labels_test_days = cell2mat({test_set.label});
        labels_test = [];
        for i = 1:num_test
            labels_test = [labels_test;...
                repmat(labels_test_days(i),num_trials_each_day(i+num_training),1)];
        end
        data.labels_test = labels_test;
    end

    displabels_test = {};
    for i = 1:num_test
        if isfield(training_set,'displabel')
            displabel_str = test_set(i).displabel;
        else
            displabel_str = num2str(days_test(i));
            for j = 1:length(attr_test{i})
                this_attr = attr_test{i}{j};
                if isnumeric(this_attr) % Trials: indicate begin and end
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
    data.displabels_test = displabels_test;
end

end