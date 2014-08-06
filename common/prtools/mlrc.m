% MLRC Muli-response Linear Regression Combiner
% 
%   W = A*(WU*MLRC)
%   W = WT*MLRC(B*WT)
%   D = C*W
%     
% INPUT
%   A   Dataset used for training base classifiers as well as combiner
%   B   Dataset used for training combiner of trained base classifiers
%   C   Dataset used for testing (executing) the combiner
%   WU  Set of untrained base classifiers, see STACKED
%   WT  Set of trained base classifiers, see STACKED
% 
% OUTPUT
%   W   Trained Muli-response Linear Regression Combiner
%   D   Dataset with prob. products (over base classifiers) per class
% 
% DESCRIPTION
% Using dataset A that contains the posterior probabilities of each instance 
% belonging to each class predicted by the base classifiers to train a
% multi-response linear regression combiner.
% If the original classification problem has K classes, it is converted
% into K seperate regression problems, where the problem for class c has
% instances with responses equal to 1 when they have label c and zero
% otherwise. Put in another way, this function establish a multi-response 
% linear regression model for each class and utilize these models to estimate 
% the probability that the instances belong to each class.
% Note that in the model for class c, only the probabilities of class c 
% predicted by the set of base classifiers are used. 
% 
% REFERENCE
% 1. Ting, KM, Witten IH. Issues in stacked generalization, Journal of  
% Artificial Intelligent Research, 1999, 10: 271-289.
% 2. Dzeroski S, Zenko B. Is combining classifiers with stacking better
% than selecting the best one? Machine Learning, 2004, 54(3): 255-273.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, STACKED, CLASSC, TESTD, LABELD

% Copyright: Chunxia Zhang, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function W = mlrc(A)

        name = 'MLR combiner';

    % If there are no inputs, return an untrained mapping.
    % Handle untrained calls like MLRC([])
    if nargin < 1 || isempty(A)
        W = prmapping(mfilename);
        W = setname(W,name);
        return
    end

    islabtype(A,'crisp');       % allow crisp labels only
    isvaldfile(A,1,2);          % at least one object per class and 2 classes

    A = testdatasize(A,'features'); % test whether they fit
    A = setprior(A,getprior(A));    % avoid many warnings
    [m,k,c] = getsize(A);           % size of training set; (m objects; k features; c classes)
    L = k/c;                        % compute the number of classifiers
    A = setfeatlab(A,repmat([1:c]',L,1)); % reset the feature labels of dataset A such that the first c features correspond to
                                          % the first classifier, the next c features correspond to the second classifier, ect.
    C = zeros(c*L,c);                     % register the coefficients of each model
    options = optimset('ToLX',1e-4);
    for i = 1:c                         % run over all classes
        Res = zeros(m,1);
        Index = find(A.featlab == i);   % find the indices correspond to the jth class
        B = seldat(A,[],Index);         % select the data corresponding to the jth class(m x L matrix)
        I = A.nlab == i;
        Res(I) = 1;
        [x,resnorm,residual,exitflag] = lsqnonneg(B.data,Res,[],options); % compute the nonnegative coefficients
        if exitflag == 0
            resnorm
        end
        
        for j = 1:L            
            C((j-1)*c+i,i) = x(j);
        end      
    end

    W = affine(C,[],A,getlablist(A),k,c);
    W = setname(W,name);

return 

