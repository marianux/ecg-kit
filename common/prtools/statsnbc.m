%STATSNBC Stats Naive Bayes Classifier (Matlab Stats Toolbox)
%
%   W = STATSNBC(A,'PARAM1',val1,'PARAM2',val2,...)
%   W = A*STATSNBC([],'PARAM1',val1,'PARAM2',val2,...)
%   D = B*W
%
% INPUT
%   A          Dataset used for training
%   PARAM1     Optional parameter, see NAIVEBAYES.FIT
%   B          Dataset used for evaluation
%
% OUTPUT
%   W          Naive Bayes classifier  
%   D          Classification matrix, dataset with posteriors
%
% DESCRIPTION
% This is the PRTools interface to the NaiveBayes classifier of the Matlab
% Stats toolbox. See there for more information. It is assumed that objects
% labels, feature labels and class priors are included in the dataset A.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, NAIVEBAYES, NAIVEBC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function W = statsnbc(varargin)

  name = 'Stats NaiveBayes';
  
  if mapping_task(varargin,'definition')
    W = define_mapping(varargin,[],name);
	elseif mapping_task(varargin,'training')
    A       = varargin{1};
    data    = +A;
    labels  = getlabels(A);
    prior   = getprior(A);
    tree    = NaiveBayes.fit(data,labels,'prior',prior, ...
    varargin{2:end});
    W = trained_mapping(A,tree);
  else % evaluation
    [A,W]    = deal(varargin{:});
    res      = getdata(W);
    post     = posterior(res,+A);
    W        = setdat(A,post,W);
  end
  
return