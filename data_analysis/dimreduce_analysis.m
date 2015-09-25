function dimreduce_analysis(md,method,training_set,test_set,opts)
% Reduce dimensionality from <number of cells> to 2 and visualize.
%
% Inputs: 
%
%   md: multiday object.
%
%   method: 'pca' or 'pls'.
%
%   training_set: struct with following fields: 'day', 'attr','label',
%   'disp_label'. 'attr' is a cell array specifying the desired trial types in
%   a day. An example 'attr' is {'start','east','correct'}. This only
%   retains the trials that start in the east arm and are correct. 'Label'
%   can be (in theory) any real numeric value, but recommended type is
%   integer. 'label' is not checked for method='pca'. 'displabel' is an
%   optional field that specifies how the element is desired to be referred
%   to in the plots. If left blank, it's set to be the aggregation of
%   'attr' elements and the corresponding day.
%
%   test_set: Optional input for testing under the model constructed with 
%   training_set. Type is same as training_set. 'label' field is optional for
%   test_set.
%
%   opts: Optional struct input whose fields overwrite defaults set in the
%   script. These are:
%       'trial_phase': Default is 1(pre-run). Can be a 1D array composed
%       of numbers 1(pre-run), 2(run), and 3(post-run).
%       'ratio_PCs': PCA is computed prior to PLS to reduce noise. This
%       parameter controls what ratio of PCs is retained for the PLS step.
%       Default is 0.5.
%       pca_dims: 2x1 array with PC dimensions to use. Default is [1,2].
%       pls_dims: same as above but for PLS.
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
else
    test_set = '';
end

if ~isfield(opts,'ratio_PCs'),  opts.ratio_PCs=0.5; end

num_training = length(training_set);
if do_test
    num_test = length(test_set);
end

data = compute_ml_inputs(md,training_set,test_set,opts);
features_training = data.features_training;
displabels_training = data.displabels_training;
if isfield(data,'labels_training')
    labels_training = data.labels_training;
end
if do_test
    features_test = data.features_test;
    displabels_test = data.displabels_test;
    if isfield(data,'labels_test'),
        labels_test = data.labels_test;
    end
end

%-----------------
% Dimension reduction
%-----------------

mean_features_training = mean(features_training,1);
[coeff,score,~,~,explained,~] = pca(features_training);

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
        centered_test_data = bsxfun(@minus,features_test,mean_features_training);
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
        centered_test_data = bsxfun(@minus,features_test,mean_features_training);
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
    hold off
end
legend('boxoff')
end
