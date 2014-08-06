%STATSDTC Stats Decision tree Classifier (Matlab Stats Toolbox)
%
%   W = STATSDTC(A,'PARAM1',val1,'PARAM2',val2,...)
%   W = A*STATSDTC([],'PARAM1',val1,'PARAM2',val2,...)
%   D = B*W
%
% INPUT
%   A          Dataset used for training
%   PARAM1     Optional parameter, see CLASSIFICATIONTREE.FIT
%   B          Dataset used for evaluation
%
% OUTPUT
%   W          Decision tree classifier  
%   D          Classification matrix, dataset with posteriors
%
% DESCRIPTION
% This is the PRTools interface to the CLASSIFICATIONTREE of the Matlab
% Stats toolbox. See there for more information. It is assumed that objects
% labels, feature labels and class priors are included in the dataset A.
%
% The decision tree is stored in W and can be retrieved by T = +W or by
% T = getdata(W). The Stats toolbox command VIEW can be used to visualize
% it, either in the command window (default) or graphically setting the
% 'mode' options to 'graph'.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, DTC, TREEC, CLASSIFICATIONTREE, VIEW

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function W = statsdtc(varargin)

  name = 'Stats DecTree';
  
  if mapping_task(varargin,'definition')
    W = define_mapping(varargin,[],name);
  elseif mapping_task(varargin,'training')
    A       = varargin{1};
    data    = +A;
    labels  = getlabels(A);
    prior   = getprior(A);
    featlab = getfeatlab(A);
    if ischar(featlab)
      featlab = cellstr(featlab);
    end
    tree    = ClassificationTree.fit(data,labels,'prior',prior, ...
    'PredictorNames',featlab,varargin{2:end});
    W = trained_mapping(A,tree);
  else % evaluation
    [A,W]    = deal(varargin{:});
    tree     = getdata(W);
    [dummy,post] = predict(tree,+A);
    W        = setdat(A,post,W);
  end
  
return