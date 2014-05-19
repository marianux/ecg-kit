function ydata = tsne(X, no_dims, initial_dims, perplexity)
%TSNE Performs symmetric t-SNE on dataset this_X
%
%   mappedX = tsne(this_X, labels, no_dims, initial_dims, perplexity)
%   mappedX = tsne(this_X, labels, initial_solution, perplexity)
%
% The function performs symmetric t-SNE on the NxD dataset this_X to reduce its 
% dimensionality to no_dims dimensions (default = 2). The data is 
% preprocessed using PCA, reducing the dimensionality to initial_dims 
% dimensions (default = 30). Alternatively, an initial solution obtained 
% from an other dimensionality reduction technique may be specified in 
% initial_solution. The perplexity of the Gaussian kernel that is employed 
% can be specified through perplexity (default = 30). The labels of the
% data are not used by t-SNE itself, however, they are used to color
% intermediate plots. Please provide an empty labels matrix [] if you
% don't want to plot results during the optimization.
% The low-dimensional data representation is returned in mappedX.
%
%

% This file is part of the Matlab Toolbox for Dimensionality Reduction.
% The toolbox can be obtained from http://homepage.tudelft.nl/19j49
% You are free to use, change, or redistribute this code in any way you
% want for non-commercial purposes. However, it is appreciated if you 
% maintain the name of the original author.
%
% (C) Laurens van der Maaten, Delft University of Technology
% adapted by MEC

    [ cant_data cant_features ] = getsize(X);

    if nargin < 2 || isempty(no_dims)
        no_dims = 2;
    end
    
    if nargin < 3 || isempty(initial_dims)
        initial_dims = min(50, cant_features);
    end
    
    if nargin < 4 || isempty(perplexity)
        perplexity = 30;
    end
    
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
        perplexity = initial_dims;
    else
        initial_solution = false;
    end
    
%     file_index = getident(X, 'file_index');
%     files_present = unique(file_index);
    
    ydata = nan(cant_data,no_dims);
    n_iter = 5000;
    
%     for ii = rowvec(files_present)
%     
%         aux_idx = find(file_index == ii);

        aux_idx = 1:cant_data;
        laux_idx = length(aux_idx);
        
        iterations = 1:n_iter:laux_idx;
        cant_iterations = length(iterations);
        
        if( cant_iterations == 1)
            iterations = [ 1  laux_idx];
        else
            cant_iterations = cant_iterations - 1;
            iterations = [ [1 iterations(2:end-1)+1 ]; iterations(2:end)  ]';
        end
    
        for jj = cant_iterations
            
            this_aux_idx = aux_idx(iterations(jj,1):iterations(jj,2));
            this_X = +(X(this_aux_idx,:));
            
            % Normalize input data
            this_X = this_X - min(this_X(:));
            this_X = this_X / max(this_X(:));
            this_X = bsxfun(@minus, this_X, mean(this_X, 1));

            % Perform preprocessing using PCA
            if ~initial_solution
                disp('Preprocessing data using PCA...');
                if size(this_X, 2) < size(this_X, 1)
                    C = this_X' * this_X;
                else
                    C = (1 / size(this_X, 1)) * (this_X * this_X');
                end
                [M, lambda] = eig(C);
                [lambda, ind] = sort(diag(lambda), 'descend');
                M = M(:,ind(1:initial_dims));
                lambda = lambda(1:initial_dims);
                if ~(size(this_X, 2) < size(this_X, 1))
                    M = bsxfun(@times, this_X' * M, (1 ./ sqrt(size(this_X, 1) .* lambda))');
                end
                this_X = this_X * M;
                clear M lambda ind
            end

            % Compute pairwise distance matrix
            sum_X = sum(this_X .^ 2, 2);
            D = bsxfun(@plus, sum_X, bsxfun(@plus, sum_X', -2 * (this_X * this_X')));

            % Compute joint probabilities
            P = d2p(D, perplexity, 1e-5);                                           % compute affinities using fixed perplexity
            clear D

            % Run t-SNE
            if initial_solution
                ydata(this_aux_idx,:) = tsne_p(P, ydata);
            else
                ydata(this_aux_idx,:) = tsne_p(P, no_dims);
            end
            
        end
        
%     end
    
    
    ydata = setdata(X, ydata);

    ydata = setfeatlab(ydata, (1:no_dims)');
    
    ui_data = getuser(ydata);
    ui_data.featlabs = char({'tSNE1' 'tSNE2'});
    ui_data.feature_properties.BeatFeaturesLabels = ui_data.featlabs;
    ui_data.feature_properties.featureTransformation = cell(no_dims,1);
    ui_data.feature_properties.debug_script_handlers = repmat({@debug_Show_ECG},no_dims,1);
    ui_data.feature_properties.Group_of_features = [];
    ui_data.feature_properties.Group_of_features_labels = [];
    
    ydata = setuser(ydata, ui_data);

    