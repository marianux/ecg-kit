%STATSKNNC Stats KNN Classifier (Matlab Stats Toolbox)
%
%   W = STATSKNNC(A,'PARAM1',val1,'PARAM2',val2,...)
%   W = A*STATSKNNC([],'PARAM1',val1,'PARAM2',val2,...)
%   D = B*W
%
% INPUT
%   A          Dataset used for training
%   PARAM1     Optional parameter, see CLASSIFICATIONKNN.FIT
%   B          Dataset used for evaluation
%
% OUTPUT
%   W          KNN classifier  
%   D          Classification matrix, dataset with posteriors (0-1)
%
% DESCRIPTION
% This is the PRTools interface to the KNN classifier of the Matlab
% Stats toolbox. See there for more information. It is assumed that objects
% labels, feature labels and class priors are included in the dataset A.
% The classification matrix D is for this classifier a 0-1 matrix with just
% a 1 in the column of the assigned class.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, CLASSIFICATIONKNN, KNNC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function W = statsknnc(varargin)

  name = 'Stats KNN';
  
  if mapping_task(varargin,'definition')
    W = define_mapping(varargin,[],name);
	elseif mapping_task(varargin,'training')
    A       = varargin{1};
    data    = +A;
    labels  = getlabels(A);
    res    = ClassificationKNN.fit(data,labels,varargin{2:end});
    W       = trained_mapping(A,res);
  else % evaluation
    [A,W]    = deal(varargin{:});
    res      = getdata(W);
    [dummy,post] = predict(res,+A);
    W        = setdat(A,post,W);
  end
  
return