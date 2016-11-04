classdef Population < handle
    properties
        num_cells
        num_trials
        num_embedding_dims
        data
        var_explained
    end

    properties(Access = private)
        W
        b
        intercept
    end

    methods
        function obj = Population(data)
        % Initializer method.
        %
        % Inputs:
        %   data: 2D array with shape [num_trials,num_cells]
        %

            obj.data = data;
            [obj.num_trials,obj.num_cells] = size(data);
            obj.b = 0;
            obj.intercept = 0;
            obj.var_explained = [];
            obj.W = eye(obj.num_cells);
            obj.num_embedding_dims = obj.num_cells;
        end
        
        
        function dim_reduce(obj,varargin)    
        % Dimension reduction (currently PCA and PLS are supported).
        %
        % Variable input name-value pairs:
        %   'num_comp': Followed by an integer denoting the number of
        %   desired reduced dimensions. If set to 0 (default), number is
        %   automatically selected to capture >90% of the variance
        %   
        %   'idx': Followed by a 1D array with the desired trial indices to
        %   use for dimension reduction. Default is all the trials
        %
        %   'method': Followed by 'PCA'(default) or 'PLS'. If PLS, then
        %   labels have to be provided (see below).
        %
        %   'labels': Followed by 1D array with the desired labels for each
        %   trial. Needed only for supervised dimension reduction.
        %

            % Defaults
            num_comp = 0;  % Get as many components as needed for 90% var
            idx = 1:obj.num_trials;  % All trials
            method = 'PCA';
            for i =1:length(varargin)
                if ischar(varargin{i})
                    switch varargin{i}
                        case 'num_comp'
                            num_comp = varargin{i+1};
                        case 'idx'
                            idx = varargin{i+1};
                        case 'method'
                            method = varargin{i+1};
                            % Clear var_explained everytime a method is
                            % called
                            obj.var_explained = [];
                        case 'labels'
                            labels = varargin{i+1};
                    end
                end
            end
            if strcmp(method,'PCA')
                [obj.b,obj.W,obj.var_explained] = do_PCA(obj.data,num_comp,idx);
            elseif strcmp(method,'PLS')
                if ~exist('labels','var')
                    error('Need labels for PLS')
                end
                [obj.b,obj.W] = do_PLS(obj.data,num_comp,idx,labels);
            elseif strcmp(method,'logistic_reg')
                if ~exist('labels','var')
                    error('Need labels for classification')
                end
                [obj.b,obj.W,obj.intercept] = do_logistic_lasso(obj.data,idx,labels);
            end
            obj.num_embedding_dims = size(obj.W,2);
        end
        
        function X = represent(obj,A,dims)
        % Representation in the reduced dimension space.
        %
        % Inputs:
        %   A: 2D array with shape [_,num_cells].
        % Optional inputs:
        %   dims = 1D array of dimensions to use for computing the
        %       reduced dimension representation.
        % Outputs:
        %   X: Reduced dimension representation of the input with shape
            if exist('dims','var')
                dims(dims>obj.num_embedding_dims) = [];  % get rid of invalid idx
            else
                dims = 1:obj.num_embedding_dims;
            end
            if  size(A,2)~=obj.num_cells
                error('Inputs need to have row dimension = num_cells')
            end
            X = bsxfun(@minus,A,obj.b)*obj.W(:,dims)+obj.intercept;
        end
        
        function d = dist(obj,A,B,dims)
        % Calculate distance between the rows of A and mean of B in the
        % representation space.
        %
        % Inputs:
        %   A: 2D array with shape [_,num_cells]
        %
        %   B: 2D array with shape [_,num_cells]
        % Optional inputs:
        %   dims = 1D array of dimensions to use for computing the
        %       reduced dimension representation.
        % Outputs:
        %   d: distances, 1D array with same length as size(A,1)
        %
            if exist('dims','var')
                dims(dims>obj.num_embedding_dims) = [];  % get rid of invalid idx
            else
                dims = 1:obj.num_embedding_dims;
            end
            if size(B,2)~= obj.num_cells || size(A,2)~=obj.num_cells
                error('Inputs need to have row dimension = num_cells')
            end
            A_reduced = obj.represent(A,dims);
            B_reduced = obj.represent(mean(B,1),dims);
            diff =  bsxfun(@minus,A_reduced,B_reduced);
            d = sum(diff.^2,2);
        end

    end


end

%------------------
% Utility functions
%------------------

function [mu,coeff,var_explained] = do_PCA(data,num_comp,idx)
    [coeff,~,latent,~,var_explained,mu] = pca(data(idx,:));
    if num_comp==0
        cum_explained = cumsum(var_explained);
        num_comp = find(cum_explained>90,1);
    end
    coeff = coeff(:,1:num_comp)*diag(sqrt(1./latent(1:num_comp)));
end

function [mu,coeff] = do_PLS(data,num_comp,idx,labels)
    % Use PCA to retain 90% of variance
    [mu,coeff,~] = do_PCA(data,0,idx);
    % If pca returns less comps than num_comp, override num_comp
    if size(coeff,2)<num_comp || num_comp==0
        num_comp = size(coeff,2);
    end
    % Get PCA scores
    PCA_scores = bsxfun(@minus,data(idx,:),mu)*coeff;
    % Do PLS
    [XL,~] = plsregress(PCA_scores,labels);
    coeff = (XL \ coeff')';
    coeff = coeff(:,1:num_comp);
end

function [mu,coeff,intercept] = do_logistic_lasso(data,idx,labels)
    % Use PCA to retain 90% of variance
    [mu,coeff] = do_PCA(data,0,idx);
    % Get PCA scores
    PCA_scores = bsxfun(@minus,data(idx,:),mu)*coeff;
    % Do lassoglm
    [B,fit_info] = lassoglm(PCA_scores,labels,'binomial','lambda',0.1,'NumLambda',1,'alpha',0.5);
    coeff = coeff*B;
    intercept = fit_info.Intercept;
end

