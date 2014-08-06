function B = vpc(A, W)
%VPC Voted perceptron classifier
%
%   W = VPC(A)
%   W = VPC(A, N)
%
% INPUT 
%   A   Dataset
%   N   Number of sweeps
%
% OUTPUT 
%   W   Voted perceptron classifier
%
% DESCRIPTION 
% The classifier trains an ensemble of perceptrons on dataset A. The
% training procedure performs N full sweeps through the training data. If N
% is not specified, 10 sweeps are performed. The total number of
% perceptrons in the ensemble is equal to N x GETSIZE(A, 1).
%
% The classification of new objects is performed by allowing the ensemble  
% of perceptrons to vote on the label of each test point.
%
% REFERENCE
% Y. Freund and R.E. Schapire. Large Margin Classification Using the 
% Perceptron Algorithm. Machine Learning 37(3):277-296, 1999.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PERLC, ADABOOSTC

% (C) Laurens van der Maaten, 2010
% University of California, San Diego

    name = 'Voted perceptron';

    % Handle untrained calls like W = vpc([]);
    if nargin == 0 || isempty(A)
        B = prmapping(mfilename); 
        B = setname(B, name); 
        return;
        
    % Handle training on dataset A (use A * vpc, A * vpc([]), and vpc(A))
    elseif (nargin == 1 && isa(A, 'prdataset')) || (isa(A, 'prdataset') && isa(W, 'double'))
        if nargin == 1
            W = 10;
        end
        islabtype(A, 'crisp');
        isvaldfile(A, 1, 2);
        A = testdatasize(A, 'features');
        A = setprior(A, getprior(A)); 
        [m, k, c] = getsize(A);
        [data.v, data.c] = train_voted_perceptron(+A, getnlab(A), W);
        B = prmapping(mfilename, 'trained', data, getlablist(A), k, c);
        B = setname(B, name);
        
    % Handle evaluation of a trained voted perceptron W for a dataset A 
    elseif (isa(A, 'prdataset') || isdouble(A)) && isa(W, 'prmapping')        
        test_post = predict_voted_perceptron(+A, W.data.v, W.data.c);
        A = prdataset(A); 
        B = setdata(A, test_post, getlabels(W));
    
    % This should not happen
    else
        error('Illegal call');
    end
end

function [v, c] = train_voted_perceptron(train_X, train_labels, max_iter)
%TRAIN_VOTED_PERCEPTRON Trains a voted perceptron on the specified data set
%
%   [v, c] = train_voted_perceptron(train_X, train_labels, max_iter)
%
% The function trains a voted perceptron on the data specified by train_X
% and train_labels. For a problem with K classes, the functions trains K
% voted perceptrons in a one-vs-all manner. The cell-array v contains the
% perceptrons for each class, whereas the cell-array c contains the
% corresponding weights of the perceptrons.
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego

    % Add the bias term to the data
    train_X = [train_X ones(size(train_X, 1), 1)];

    % Initialize some variables
    [n, d] = size(train_X);
    lablist = unique(train_labels);
    no_labels = length(lablist);
    
    % Shuffle the data
    ind = randperm(n);
    train_X = train_X(ind,:);
    train_labels = train_labels(ind);
    
    % Initialize voted perceptron cell arrays
    c = cell(no_labels, 1);
    v = cell(no_labels, 1);
    
    % Loop over all classes
    for k=1:no_labels
        
        % Initialize voted perceptron
        v{k} = zeros(d, ceil(n * max_iter / 4));          % perceptrons
        c{k} = zeros(1, ceil(n * max_iter / 4));          % survival times
        index = 1;
        
        % Construct appropriate labelling
        labels = train_labels;
        labels(labels ~= k) = -1;
        labels(labels == k) =  1;

        % Perform learning iterations
        for iter=1:max_iter 

            % Loop over all data points
            for i=1:n

                % Perform prediction
                y = +sign(train_X(i,:) * v{k}(:,index));

                % Train new perceptron if misclassified, or increase weight
                if y == labels(i)
                    c{k}(index) = c{k}(index) + 1;
                else
                    v{k}(:,index + 1) = v{k}(:,index) + labels(i) .* train_X(i,:)';
                    c{k}(  index + 1) = 1;
                    index = index + 1;
                end
            end
        end
        
        % Make sure we do not have 'empty' perceptrons left
        v{k} = v{k}(:,1:index);
        c{k} = c{k}(  1:index);
    end
end

function test_post = predict_voted_perceptron(test_X, v, c)
%PREDICT_VOTED_PERCEPTRON Classify test data using a voted perceptron
%
%   test_post = predict_voted_perceptron(test_X, v, c, lablist)
%
% The function classifies the data points in test_X using the voted
% perceptron specified by v and c. The voted perceptron can be trained
% using the TRAIN_VOTED_PERCEPTRON function.
%
%
% (C) Laurens van der Maaten, 2010
% University of California, San Diego

    
    % Add the bias term to the data
    test_X = [test_X ones(size(test_X, 1), 1)];

    % Prediction is maximum over all voted perceptrons
    test_post = zeros(size(test_X, 1), length(v));
    for k=1:length(v)
        test_post(:,k) = sum(bsxfun(@times, c{k}, +sign(test_X * v{k})), 2);
    end
    
    % Estimate posterior if requested
    for k=1:length(c)
        test_post(:,k) = test_post(:,k) + sum(c{k});
    end
    test_post = bsxfun(@rdivide, test_post, sum(test_post, 2));
end
    
