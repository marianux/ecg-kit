function ydata = tsne_p(P, no_dims)
%TSNE_P Performs symmetric t-SNE on affinity matrix P
%
%   mappedX = tsne_p(P, labels, no_dims)
%
% The function performs symmetric t-SNE on pairwise similarity matrix P 
% to create a low-dimensional map of no_dims dimensions (default = 2).
% The matrix P is assumed to be symmetric, sum up to 1, and have zeros
% on the diagonal.
% The labels of the data are not used by t-SNE itself, however, they 
% are used to color intermediate plots. Please provide an empty labels
% matrix [] if you don't want to plot results during the optimization.
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



%     if ~exist('labels', 'var')
%         labels = [];
%     end

    if nargin < 2 || isempty(no_dims)
        no_dims = 2;
    end
    
    % First check whether we already have an initial solution
    if numel(no_dims) > 1
        initial_solution = true;
        ydata = no_dims;
        no_dims = size(ydata, 2);
    else
        initial_solution = false;
    end
    
    % Initialize some variables
    n = size(P, 1);                                     % number of instances
    momentum = 0.5;                                     % initial momentum
    final_momentum = 0.8;                               % value to which momentum is changed
    mom_switch_iter = 250;                              % iteration at which momentum is changed
    stop_lying_iter = 100;                              % iteration at which lying about P-values is stopped
    max_iter = 400;                                    % maximum number of iterations
    epsilon = 500;                                      % initial learning rate
    min_gain = .01;                                     % minimum gain for delta-bar-delta
    
    % Make sure P-vals are set properly
    P(1:n + 1:end) = 0;                                 % set diagonal to zero
    P = 0.5 * (P + P');                                 % symmetrize P-values
    P = max(P ./ sum(P(:)), realmin);                   % make sure P-values sum to one
    const = sum(P(:) .* log(P(:)));                     % constant in KL divergence
    if ~initial_solution
        P = P * 4;                                      % lie about the P-vals to find better local minima
    end
    
    % Initialize the solution
    if ~initial_solution
        ydata = .0001 * randn(n, no_dims);
    end
    y_incs  = zeros(size(ydata));
    gains = ones(size(ydata));
    
    % Run the iterations
%     for iter=1:max_iter
    iter = 1;
    cost = 1;
    initial_cost = 0;
    k_cost = 0.1;
    while iter <= max_iter && cost > (k_cost*initial_cost)
        
        % Compute joint probability that point i and j are neighbors
        sum_ydata = sum(ydata .^ 2, 2);
        num = 1 ./ (1 + bsxfun(@plus, sum_ydata, bsxfun(@plus, sum_ydata', -2 * (ydata * ydata')))); % Student-t distribution
        num(1:n+1:end) = 0;                                                 % set diagonal to zero
        Q = max(num ./ sum(num(:)), realmin);                               % normalize to get probabilities
        
        % Compute the gradients (faster implementation)
        L = (P - Q) .* num;
        y_grads = 4 * (diag(sum(L, 1)) - L) * ydata;
            
        % Update the solution
        gains = (gains + .2) .* (sign(y_grads) ~= sign(y_incs)) ...         % note that the y_grads are actually -y_grads
              + (gains * .8) .* (sign(y_grads) == sign(y_incs));
        gains(gains < min_gain) = min_gain;
        y_incs = momentum * y_incs - epsilon * (gains .* y_grads);
        ydata = ydata + y_incs;
        ydata = bsxfun(@minus, ydata, mean(ydata, 1));
        
        % Update the momentum if necessary
        if iter == mom_switch_iter
            momentum = final_momentum;
        end
        if iter == stop_lying_iter && ~initial_solution
            P = P ./ 4;
        end
        
        % Print out progress
        if ~rem(iter, 10)
            cost = const - sum(P(:) .* log(Q(:)));
            disp(['Iteration ' num2str(iter) ': error is ' num2str(cost)]);
        end
        
        if initial_cost == 0
            cost = const - sum(P(:) .* log(Q(:)));
            initial_cost = cost;
        end

        iter = iter + 1;
        
%         % Display scatter plot (maximally first three dimensions)
%         if ~isempty(labels)
%             if no_dims == 1
%                 scatter(ydata, ydata, 9, labels, 'filled');
%             elseif no_dims == 2
%                 scatter(ydata(:,1), ydata(:,2), 9, labels, 'filled');
%             else
%                 scatter3(ydata(:,1), ydata(:,2), ydata(:,3), 40, labels, 'filled');
%             end
%             axis equal tight
% %             axis off
%             drawnow
%         end
        
        
    end
    