%LINEARR Trainable mapping for linear regression
%
%     Y = LINEARR(X,LAMBDA,N)
%     Y = X*LINEARR([],LAMBDA,N)
%     Y = X*LINEARR(LAMBDA,N)
%
% INPUT
%   X       Dataset
%   LAMBDA  Regularization parameter (default: no regularization)
%   N       Order of polynomial (default: 1)
%
% OUTPUT
%   Y       Linear (or higher order) regression
%
% DESCRIPTION
% Perform a linear regression on dataset X, with regularization
% parameter LAMBDA. When N is supplied, also higher order polynomials
% are possible.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% RIDGER, TESTR, PLOTR, VANDERMONDEM

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = linearr(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],[],1);
  
  if mapping_task(argin,'definition')
    y = define_mapping(argin,'untrained');
    y = setname(y,'Linear regression');
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [x,lambda,p] = deal(argin{:});
    [n,d] = size(x);
    X = +vandermondem(x,p);
    if isempty(lambda)
      beta = prinv(X'*X)*X'*gettargets(x);
    else
      dimp = size(X,2);
      beta = prinv(X'*X + lambda*eye(dimp))*X'*gettargets(x);
    end
    W.beta = beta;
    W.n = p;
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'Linear regression');
  else % Evaluation
    [x,v] = deal(argin{1:2});
    w = getdata(v);
    out = vandermondem(x,w.n)*w.beta;
    y = setdat(x,out);

  end
