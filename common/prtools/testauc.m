%TESTAUC Multiclass error area under the ROC
%
%   E = TESTAUC(A*W)
%   E = TESTAUC(A,W)
%   E = A*W*TESTAUC
%
% INPUT
%   A  Dataset to be classified
%   W  Classifier
%
% OUTPUT
%   E  Error, Area under the ROC
%
% DESCRIPTION
% The area under the ROC is computed for the datset A w.r.t. the
% classifer W. The estimator is based on a rank analysis of the classifier
% outcomes. Ties are broken by a two-way sorting and averaging. 
%
% The multiclass situation is solved by averaging over all outcomes of
% the one-against-rest ROCs.
%
% Note that E is an error and not a performance measure like the AUC often
% used in literature.
%
% SEE ALSO
% DATASETS, MAPPINGS, TESTC, ROC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function e = testauc(a,w)

	
	if (nargin == 0) | (isempty(a))
		% No input arguments given: return mapping information.
		e = prmapping(mfilename,'fixed');
        return
	elseif (nargin == 1)
		% Classification matrix already computed
        d = a;
    else
        % Compute classification matrix now
        d = a*w;
    end

    [m,k,c] = getsize(d);
    s = classsizes(d);

    if k == 1 % classifier with a single class outcome, make two for consistency
        d = [d 1-d];
        k = 2;
    end

    class_names = getfeatlab(d);   % class names
    e = zeros(1,c);
    for j=1:c                      % run over all classes
        % forward sorting
        [ds,J1] = sort(-d(:,j));  [j1,K1] = sort(J1);
        % backward sorting to solve ties
        [ds,J2] = sort(flipud(-d(:,j))); [j2,K2] = sort(J2); K2 = flipud(K2);
        % get all object indices for this class
        K = findlabels(d,class_names(j,:));
        % retrieve number of wrong pairs
        e(j) = (sum(K1(K)) + sum(K2(K))-(s(j)*(s(j)+1)))/2;
        % error contribution
        e(j) = e(j) / ((m-s(j))*s(j));
    end
    % average over all classes
    e = mean(e);
