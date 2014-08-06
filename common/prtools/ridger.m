%RIDGER Trainable mapping for ridge regression
%
%     W = RIDGER(X,LAMBDA)
%     W = X*RIDGER([],LAMBDA)
%     W = X*RIDGER(LAMBDA)
%
% INPUT
%   X        Regression dataset
%   LAMBDA   Regularization parameter (default LAMBDA=1)
%
% OUTPUT
%   W        Ridge regression mapping
%
% DESCRIPTION
% Perform a ridge regression on dataset X, with the regularization
% parameter LAMBDA.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% LASSOR, PLOTR, LINEARR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = ridger(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    y = define_mapping(argin,'fixed');
    y = setname(y,'Ridge regression');
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [x,lambda] = deal(argin{:});
    [n,d] = size(x);
    y = gettargets(x);
    X = +x;
    beta = prpinv(X'*X + diag(repmat(lambda,d,1)))*X'*(y-mean(y));
    W = [mean(y); beta];   % don't forget the offset
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'Ridge regression');
    
  else	% Evaluation
    
    [x,v] = deal(argin{1:2});
    [n,d] = size(x);
    w = getdata(v);
    out = [ones(n,1) +x]*w;
    y = setdat(x,out);
    
  end
	
end
