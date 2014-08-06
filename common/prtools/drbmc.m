function [B, ll] = drbmc(varargin)
%DRBMC Trainable classifier by Discriminative Restricted Boltzmann Machine
%
%   W = DRBMC(A,N,L)
%   W = A*DRBMC([],N,L)
%   W = A*DRBMC(N,L)
%
% INPUT 
%   A   Dataset
%   N   Number of hidden units
%   L   Regularization parameter (L2)
%
% OUTPUT 
%   W   Discriminative Restricted Boltzmann Machine classifier
%
% DESCRIPTION 
% The classifier trains a discriminative Restricted Boltzmann Machine (RBM)
% on dataset A. The discriminative RBM can be viewed as a logistic
% regressor with hidden units. The discriminative RBM has N hidden units
% (default = 5). It is trained with L2 regularization using regularization
% parameter L (default = 0).
%
% New objects are classified in the same way as in a logistic regressor,
% using a softmax function over the labels.
%
% REFERENCE
% H. Larochelle and Y. Bengio. Classification using Discriminative Restricted 
% Boltzmann Machines. Proceedings of the 25th International Conference on 
% Machine Learning (ICML), pages 536?543, 2008.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, LOGLC

% (C) Laurens van der Maaten, 2010
% University of California, San Diego

% 28 Jan 2013, debugged, default N = 5, PRTools call tests

    name = 'Discr. RBM';
    argin = shiftargin(varargin,'scalar');
    argin = setdefaults(argin,[],5,0,1e-3);

    % Handle untrained calls like W = drbmc([]);
    if mapping_task(argin,'definition')
        B = define_mapping(argin,'untrained',name);
        %B = prmapping(mfilename); 
        %B = setname(B, name); 
        return;
        
    % Handle training on dataset A (use A * drbmc, A * drbmc([]), and drbmc(A))
    elseif mapping_task(argin,'training')
      [A,W,L,I] = deal(argin{:});
        if ~isa(L, 'double')
          error('Illegal call');
        end
        if ~isa(I, 'double')
            error('Illegal call');
        end
        islabtype(A, 'crisp');
        isvaldfile(A, 1, 2);
        A = testdatasize(A, 'features');
        A = setprior(A, getprior(A)); 
        [m, k, c] = getsize(A);
        
        % Normalize data and add extra feature
        train_X = +A;
        data.min_X = min(train_X(:));
        train_X = train_X - data.min_X;
        data.max_X = max(train_X(:));
        train_X = train_X / data.max_X;
        train_X = (train_X * 2) - 1;
        data.max_radius = max(sum(train_X .^ 2, 2));
        train_X = [train_X sqrt(data.max_radius - sum(train_X .^ 2, 2))];
        
        % Train the dRBM
        [data.machine, ll] = train_drbm(train_X, getnlab(A), W, 100, L, I);
        B = prmapping(mfilename, 'trained', data, getlablist(A), k, c);
        B = setname(B, name);
        
    % Handle evaluation of a trained dRBM W for a dataset A 
    elseif mapping_task(argin,'execution')
        [A,W] = deal(argin{1:2});
        % Normalize data and add extra feature
        test_X = +A;
        test_X = test_X - W.data.min_X;
        test_X = test_X / W.data.max_X;
        test_X = (test_X * 2) - 1;
        test_X = [test_X sqrt(W.data.max_radius - sum(test_X .^ 2, 2))];
        
        % Evaluate dRBM
        [test_labels, test_post] = predict_ff(W.data.machine, test_X, length(W.labels));
        A = prdataset(A); 
        B = setdata(A, real(test_post), getlabels(W));
        ll = [];
    
    % This should not happen
    else
        error('Illegal call');
    end
end

function [machine, ll] = train_drbm(X, labels, h, max_iter, weight_cost, variance)
%TRAIN_DRBM Trains a discrimative RBM using gradient descent
%
%   machine = train_drbm(X, labels, h, max_iter, weight_cost, variance)
%
% Trains a discriminative Restricted Boltzmann Machine on dataset X. The RBM
% has h hidden nodes (default = 100). The training is performed by means of
% the gradient descent. The visible and hidden units are binary stochastic, 
% whereas the label layer has softmax units. The maximum number of 
% iterations can be specified through max_iter (default = 30).
% The trained RBM is returned in the machine struct.
%
%
% (C) Laurens van der Maaten
% University of California San Diego, 2010

    % Process inputs
    if ischar(X)
        load(X);
    end
    if ~exist('h', 'var') || isempty(h)
        h = 100;
    end
    if ~exist('max_iter', 'var') || isempty(max_iter)
        max_iter = 15;        
    end
    if ~exist('weight_cost', 'var') || isempty(weight_cost)
        weight_cost = 0;
    end 
    
    % Randomly permute data
    ind = randperm(size(X, 1));
    X = X(ind,:);
    labels = labels(ind);
    
    % Convert labels to binary matrix
    no_labels = length(unique(labels));
    label_matrix = zeros(length(labels), no_labels);
    label_matrix(sub2ind(size(label_matrix), (1:length(labels))', labels)) = 1;
    
    % Initialize some variables
    [n, v] = size(X);
    machine.W = randn(v, h) * variance;
    machine.labW = randn(no_labels, h) * variance;
    machine.bias_upW = zeros(1, h);
    machine.bias_labW = zeros(1, no_labels); 
    
    % Main loop
    for iter=1:max_iter        
        
        % Encode current solution
        x = [machine.W(:); machine.labW(:); machine.bias_upW(:); machine.bias_labW(:)];

        % Run conjugate gradients using five linesearches
%         checkgrad('drbm_grad', x, 1e-7, machine, cur_X, cur_lab, cur_labmat, weight_cost)
        x = minimize(x, 'drbm_grad', 9, machine, X, labels, label_matrix, weight_cost);

        % Process the updated solution
        ii = 1;
        machine.W = reshape(x(ii:ii - 1 + numel(machine.W)), size(machine.W));
        ii = ii + numel(machine.W);
        machine.labW = reshape(x(ii:ii - 1 + numel(machine.labW)), size(machine.labW));
        ii = ii + numel(machine.labW);                
        machine.bias_upW = reshape(x(ii:ii - 1 + numel(machine.bias_upW)), size(machine.bias_upW));
        ii = ii + numel(machine.bias_upW);
        machine.bias_labW = reshape(x(ii:ii - 1 + numel(machine.bias_labW)), size(machine.bias_labW));
    end
    machine.type = 'drbm';
    ll = drbm_grad(x, machine, X, labels, label_matrix, weight_cost);
end

function [test_labels, test_post] = predict_ff(machine, test_X, no_labels)
%PREDICT_FF Perform label predictions in a feed-forward neural network
%
%   [test_labels, test_post] = predict_ff(network, test_X, no_labels)
%
% The function performs label prediction in the feed-forward neural net
% specified in network (trained using TRAIN_DBN) on dataset test_X. The
% predicted labels are returned in test_labels. The class posterior is
% returned in test_post.
%
%
% (C) Laurens van der Maaten
% University of California San Diego, 2010

    % Precompute W dot X + bias
    Wx = bsxfun(@plus, test_X * machine.W, machine.bias_upW);

    % Loop over all labels to compute p(y|x)
    n = size(test_X, 1);
    energy = zeros(n, size(machine.labW, 2), no_labels);
    label_logpost = zeros(n, no_labels);
    for i=1:no_labels

        % Clamp current label
        lab = zeros(1, no_labels);
        lab(i) = 1;

        % Compute log-posterior probability for current class
        energy(:,:,i) = bsxfun(@plus, lab * machine.labW, Wx);
        label_logpost(:,i) = machine.bias_labW(i) + sum(log(1 + exp(energy(:,:,i))), 2);
    end
    test_post = exp(bsxfun(@minus, label_logpost, max(label_logpost, [], 2)));
    test_post = bsxfun(@rdivide, test_post, sum(test_post, 2));
    [foo, test_labels] = max(test_post, [], 2);
end

function [C, dC] = drbm_grad(x, machine, X, labels, label_matrix, weight_cost)
%DRBM_GRAD Computes gradient for training a discriminative RBM
%
%   [C, dC] = drbm_grad(x, machine, X, labels, label_matrix)
%
% Computes gradient for training a discriminative RBM.
%
%
% (C) Laurens van der Maaten
% University of California San Diego, 2010

    % Decode the current solution
    if ~isempty(x)
        ii = 1;
        machine.W = reshape(x(ii:ii - 1 + numel(machine.W)), size(machine.W));
        ii = ii + numel(machine.W);
        machine.labW = reshape(x(ii:ii - 1 + numel(machine.labW)), size(machine.labW));
        ii = ii + numel(machine.labW);                
        machine.bias_upW = reshape(x(ii:ii - 1 + numel(machine.bias_upW)), size(machine.bias_upW));
        ii = ii + numel(machine.bias_upW);
        machine.bias_labW = reshape(x(ii:ii - 1 + numel(machine.bias_labW)), size(machine.bias_labW));
    end
    no_labels = size(label_matrix, 2);
    n = size(X, 1);

    % Precompute W dot X + bias
    Wx = bsxfun(@plus, X * machine.W, machine.bias_upW);

    % Loop over all labels to compute p(y|x)
    energy = zeros(n, size(machine.labW, 2), no_labels);
    label_logpost = zeros(n, no_labels);
    for i=1:no_labels

        % Clamp current label
        lab = zeros(1, no_labels);
        lab(i) = 1;

        % Compute log-posterior probability for current class
        energy(:,:,i) = bsxfun(@plus, lab * machine.labW, Wx);
        label_logpost(:,i) = machine.bias_labW(i) + sum(log(1 + exp(energy(:,:,i))), 2);
    end
    label_post = exp(bsxfun(@minus, label_logpost, max(label_logpost, [], 2)));
    label_post = bsxfun(@rdivide, label_post, sum(label_post, 2));

    % Compute cost function -sum(log p(y|x))
    Pyx = log(sum(label_matrix .* label_post, 2));
    C = -sum(Pyx) + weight_cost * sum(machine.W(:) .^ 2);
    
    % Only compute gradients if they are requested
    if nargout > 1
        
        % Precompute sigmoids of energies, and products with label posterior
        sigmoids = 1 ./ (1 + exp(-energy));
        sigmoids_label_post = sigmoids;
        for i=1:no_labels
            sigmoids_label_post(:,:,i) = bsxfun(@times, sigmoids_label_post(:,:,i), label_post(:,i));
        end

        % Initialize gradients
        grad_W          = zeros(size(machine.W));
        grad_labW       = zeros(size(machine.labW));
        grad_bias_upW   = zeros(size(machine.bias_upW));

        % Sum gradient with respect to the data-hidden weights
        for i=1:no_labels
            ind = find(labels == i);
            grad_W = grad_W + X(ind,:)' * sigmoids(ind,:,i) - X' * sigmoids_label_post(:,:,i);
        end
        grad_W = grad_W - 2 * weight_cost * machine.W;

        % Compute gradient with respect to bias on label-hidden weights
        for i=1:no_labels
            grad_labW(i,:) = sum(sigmoids(labels == i,:,i), 1) - sum(sigmoids_label_post(:,:,i), 1);
        end

        % Compute gradient with respect to bias on hidden units
        for i=1:no_labels
            grad_bias_upW = grad_bias_upW + sum(sigmoids(labels == i,:,i), 1) - sum(sigmoids_label_post(:,:,i), 1);
        end

        % Compute the gradient with respect to bias on labels
        grad_bias_labW = sum(label_matrix - label_post, 1);

        % Encode the gradient
        dC = -[grad_W(:); grad_labW(:); grad_bias_upW(:); grad_bias_labW(:)];
    end
end 
   

function [X, fX, i] = minimize(X, f, length, P1, P2, P3, P4, P5);

% Minimize a continuous differentialble multivariate function. Starting point
% is given by "X" (D by 1), and the function named in the string "f", must
% return a function value and a vector of partial derivatives. The Polack-
% Ribiere flavour of conjugate gradients is used to compute search directions,
% and a line search using quadratic and cubic polynomial approximations and the
% Wolfe-Powell stopping criteria is used together with the slope ratio method
% for guessing initial step sizes. Additionally a bunch of checks are made to
% make sure that exploration is taking place and that extrapolation will not
% be unboundedly large. The "length" gives the length of the run: if it is
% positive, it gives the maximum number of line searches, if negative its
% absolute gives the maximum allowed number of function evaluations. You can
% (optionally) give "length" a second component, which will indicate the
% reduction in function value to be expected in the first line-search (defaults
% to 1.0). The function returns when either its length is up, or if no further
% progress can be made (ie, we are at a minimum, or so close that due to
% numerical problems, we cannot get any closer). If the function terminates
% within a few iterations, it could be an indication that the function value
% and derivatives are not consistent (ie, there may be a bug in the
% implementation of your "f" function). The function returns the found
% solution "X", a vector of function values "fX" indicating the progress made
% and "i" the number of iterations (line searches or function evaluations,
% depending on the sign of "length") used.
%
% Usage: [X, fX, i] = minimize(X, f, length, P1, P2, P3, P4, P5)
%
% Copyright (C) 2001 and 2002 by Carl Edward Rasmussen. Date 2002-02-13

    RHO = 0.01;                            % a bunch of constants for line searches
    SIG = 0.5;       % RHO and SIG are the constants in the Wolfe-Powell conditions
    INT = 0.1;    % don't reevaluate within 0.1 of the limit of the current bracket
    EXT = 3.0;                    % extrapolate maximum 3 times the current bracket
    MAX = 20;                         % max 20 function evaluations per line search
    RATIO = 100;                                      % maximum allowed slope ratio

    argstr = [f, '(X'];                      % compose string used to call function
    for i = 1:(nargin - 3)
    argstr = [argstr, ',P', int2str(i)];
    end
    argstr = [argstr, ')'];

    if max(size(length)) == 2, red=length(2); length=length(1); else red=1; end
    if length>0, S=['Linesearch']; else S=['Function evaluation']; end

    i = 0;                                            % zero the run length counter
    ls_failed = 0;                             % no previous line search has failed
    fX = [];
    [f1 df1] = eval(argstr);                      % get function value and gradient
    i = i + (length<0);                                            % count epochs?!
    s = -df1;                                        % search direction is steepest
    d1 = -s'*s;                                                 % this is the slope
    z1 = red/(1-d1);                                  % initial step is red/(|s|+1)

    while i < abs(length)                                      % while not finished
    i = i + (length>0);                                      % count iterations?!

    X0 = X; f0 = f1; df0 = df1;                   % make a copy of current values
    X = X + z1*s;                                             % begin line search
    [f2 df2] = eval(argstr);
    i = i + (length<0);                                          % count epochs?!
    d2 = df2'*s;
    f3 = f1; d3 = d1; z3 = -z1;             % initialize point 3 equal to point 1
    if length>0, M = MAX; else M = min(MAX, -length-i); end
    success = 0; limit = -1;                     % initialize quanteties
    while 1
        while ((f2 > f1+z1*RHO*d1) | (d2 > -SIG*d1)) & (M > 0)
            limit = z1;                                         % tighten the bracket
            if f2 > f1
                z2 = z3 - (0.5*d3*z3*z3)/(d3*z3+f2-f3);                 % quadratic fit
            else
                A = 6*(f2-f3)/z3+3*(d2+d3);                                 % cubic fit
                B = 3*(f3-f2)-z3*(d3+2*d2);
                z2 = (sqrt(B*B-A*d2*z3*z3)-B)/A;       % numerical error possible - ok!
            end
            if isnan(z2) | isinf(z2)
                z2 = z3/2;                  % if we had a numerical problem then bisect
            end
            z2 = max(min(z2, INT*z3),(1-INT)*z3);  % don't accept too close to limits
            z1 = z1 + z2;                                           % update the step
            X = X + z2*s;
            [f2 df2] = eval(argstr);
            M = M - 1; i = i + (length<0);                           % count epochs?!
            d2 = df2'*s;
            z3 = z3-z2;                    % z3 is now relative to the location of z2
        end
        if f2 > f1+z1*RHO*d1 | d2 > -SIG*d1
            break;                                                % this is a failure
        elseif d2 > SIG*d1
            success = 1; break;                                             % success
        elseif M == 0
            break;                                                          % failure
        end
        A = 6*(f2-f3)/z3+3*(d2+d3);                      % make cubic extrapolation
        B = 3*(f3-f2)-z3*(d3+2*d2);
        z2 = -d2*z3*z3/(B+sqrt(B*B-A*d2*z3*z3));        % num. error possible - ok!
        if ~isreal(z2) | isnan(z2) | isinf(z2) | z2 < 0   % num prob or wrong sign?
            if limit < -0.5                               % if we have no upper limit
                z2 = z1 * (EXT-1);                 % the extrapolate the maximum amount
            else
                z2 = (limit-z1)/2;                                   % otherwise bisect
            end
        elseif (limit > -0.5) & (z2+z1 > limit)          % extraplation beyond max?
            z2 = (limit-z1)/2;                                               % bisect
        elseif (limit < -0.5) & (z2+z1 > z1*EXT)       % extrapolation beyond limit
            z2 = z1*(EXT-1.0);                           % set to extrapolation limit
        elseif z2 < -z3*INT
            z2 = -z3*INT;
        elseif (limit > -0.5) & (z2 < (limit-z1)*(1.0-INT))   % too close to limit?
            z2 = (limit-z1)*(1.0-INT);
        end
        f3 = f2; d3 = d2; z3 = -z2;                  % set point 3 equal to point 2
        z1 = z1 + z2; X = X + z2*s;                      % update current estimates
        [f2 df2] = eval(argstr);
        M = M - 1; i = i + (length<0);                             % count epochs?!
        d2 = df2'*s;
    end                                                      % end of line search

    if success                                         % if line search succeeded
        f1 = f2; fX = [fX' f1]';
        s = (df2'*df2-df1'*df2)/(df1'*df1)*s - df2;      % Polack-Ribiere direction
        tmp = df1; df1 = df2; df2 = tmp;                         % swap derivatives
        d2 = df1'*s;
        if d2 > 0                                      % new slope must be negative
            s = -df1;                              % otherwise use steepest direction
            d2 = -s'*s;
        end
        z1 = z1 * min(RATIO, d1/(d2-realmin));          % slope ratio but max RATIO
        d1 = d2;
        ls_failed = 0;                              % this line search did not fail
    else
        X = X0; f1 = f0; df1 = df0;  % restore point from before failed line search
        if ls_failed | i > abs(length)          % line search failed twice in a row
            break;                             % or we ran out of time, so we give up
        end
        tmp = df1; df1 = df2; df2 = tmp;                         % swap derivatives
        s = -df1;                                                    % try steepest
        d1 = -s'*s;
        z1 = 1/(1-d1);
        ls_failed = 1;                                    % this line search failed
    end
    end
end

function d = checkgrad(f, X, e, P1, P2, P3, P4, P5, P6);

% checkgrad checks the derivatives in a function, by comparing them to finite
% differences approximations. The partial derivatives and the approximation
% are printed and the norm of the diffrence divided by the norm of the sum is
% returned as an indication of accuracy.
%
% usage: checkgrad('f', X, e, P1, P2, ...)
%
% where X is the argument and e is the small perturbation used for the finite
% differences. and the P1, P2, ... are optional additional parameters which
% get passed to f. The function f should be of the type 
%
% [fX, dfX] = f(X, P1, P2, ...)
%
% where fX is the function value and dfX is a vector of partial derivatives.
%
% Carl Edward Rasmussen, 2001-08-01.

argstr = [f, '(X'];                            % assemble function call strings
argstrd = [f, '(X+dx'];
for i = 1:(nargin - 3)
  argstr = [argstr, ',P', int2str(i)];
  argstrd = [argstrd, ',P', int2str(i)];
end
argstr = [argstr, ')'];
argstrd = [argstrd, ')'];

[y dy] = eval(argstr);                         % get the partial derivatives dy

dh = zeros(length(X),1) ;
for j = 1:length(X)
  dx = zeros(length(X),1);
  dx(j) = dx(j) + e;                               % perturb a single dimension
  y2 = eval(argstrd);
  dx = -dx ;
  y1 = eval(argstrd);
  dh(j) = (y2 - y1)/(2*e);
end

disp([dy dh])                                        % print the two vectors
d = norm(dh-dy)/norm(dh+dy);       % return norm of diff divided by norm of sum
end
    