%TESTR MSE for regression
%
%      E = TESTR(X,W,TYPE)
%      E = TESTR(X*W,TYPE)
%      E = X*W*TESTR([],TYPE)
%      E = X*W*TESTR(TYPE)
%
% INPUT
%   X    Regression dataset
%   W    Regression mapping
%   TYPE Type of error measure, default: mean squared error
%
% OUTPUT
%   E    Mean squared error
%
% DESCRIPTION
% Compute the error of regression W on dataset X. The following error
% measures have been defined for TYPE:
% 'mse'    mean squared error (default)
% 'mad'    mean absolute deviation
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  RSQUARED, TESTC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function e = testr(varargin)
  
  argin = shiftargin(varargin,'char');
  argin = setdefaults(argin,[],[],'mse');
  
  if mapping_task(argin,'definition')
    e = define_mapping(argin,'fixed');
    
  else	% Evaluate.
  
    [x,w,type] = deal(argin{:}); 

    if (ismapping(w) & istrained(w))
      x = x*w;
    end
    if ischar(w)
      type = w;
    end
    switch type
      case 'mse'
        e = mean((+x(:,1) - gettargets(x)).^2);
      case 'mad'
        e = mean(abs(+x(:,1) - gettargets(x)));
      otherwise
        error('Error %s is not implemented.',type);
    end

    if nargout==0
      %display results on the screen:
      fprintf('Error on %d objects: %f.\n',...
        size(x,1), e);
      clear e;
    end
    
  end

