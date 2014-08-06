%STATSLINC Stats Linear Classifier (Matlab Stats Toolbox)
%
%   W = STATSLINC(A,'PARAM1',val1,'PARAM2',val2,...)
%   W = A*STATSLINC([],'PARAM1',val1,'PARAM2',val2,...)
%   D = B*W
%
% INPUT
%   A          Dataset used for training
%   PARAM1     Optional parameter, see CLASSIFICATIONDISCRIMINANT.FIT
%   B          Dataset used for evaluation
%
% OUTPUT
%   W          Linear classifier  
%   D          Classification matrix, dataset with posteriors
%
% DESCRIPTION
% This is the PRTools interface to the CLASSIFICATIONDISCRIMINANT of the 
% Matlab Stats toolbox. See there for more information. It is assumed that 
% objects labels, feature labels and class priors are included in the 
% dataset A.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, CLASSIFICATIONDISCRIMINANT, FISHERC, LDC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function W = statslinc(varargin)

  name = 'Stats LinClassf';
  
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
    tree    = ClassificationDiscriminant.fit(data,labels,'prior',prior, ...
    'PredictorNames',featlab,varargin{2:end});
    W = trained_mapping(A,tree);
  else % evaluation
    [A,W]    = deal(varargin{:});
    tree     = getdata(W);
    [dummy,post] = predict(tree,+A);
    W        = setdat(A,post,W);
  end
  
return