%LASSOR Trainable mapping for LASSO regression
%
%    W = LASSOR(X,LAMBDA)
%    W = X*LASSOR([],LAMBDA)
%    W = X*LASSOR(LAMBDA)
%
% INPUT
%   X       Regression dataset
%   LAMBDA  Regularization parameter
%
% OUTPUT
%   W       LASSO regression mapping
%
% DESCRIPTION
%  The 'Least Absolute Shrinkage and Selection Operator' regression,
%  using the regularization parameter LAMBDA.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  RIDGER, LINEARR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = lassor(varargin)

  mapname = 'LASSO regression';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    
    y = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.
 
    [x,lambda] = deal(argin{:});
    [n,d] = size(x);
    y = gettargets(x);
    W = arrfit(+x,(y-mean(y)),lambda);
    W = [mean(y); W];
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,mapname);
    
  else                                     % Evaluation
    
    [x,v] = deal(argin{1:2});
    w = getdata(v);
    [n,d] = size(x);
    out = [ones(n,1) +x]*w;
    y = setdat(x,out);

  end
  
