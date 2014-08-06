%KSMOOTHR Kernel smoother
%
%      W = KSMOOTHR(X,H)
%      W = X*KSMOOTHR([],H)
%      W = X*KSMOOTHR(X,H)
%
% INPUT
%   X    Regression dataset
%   H    Width parameter (default H=1)
%
% OUTPUT
%   W    Kernel smoother mapping
%
% DESCRIPTION
% Train a kernel smoothing W on data X, with width parameter H.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  KNNR, TESTR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = ksmoothr(varargin)

  mapname = 'Kernel smoother';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    
    y = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.

    [x,h] = deal(argin{:});
    [n,d] = size(x);
    W.x = +x;
    W.y = gettargets(x);
    W.h = h;
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'Kernel smoother');
    
  else
    
    % evaluation
    [x,v] = deal(argin{:});
    W = getdata(v);
    [n,d] = size(x);
    m = size(W.x,1);
    xtst = +x;
    gamma = -1/(W.h*W.h); % tiny speedup
    % now go through all test data:
    y = zeros(n,1);
    K = exp(gamma*distm(xtst,W.x));
    y = (K*W.y)./sum(K,2);
    y = setdat(x,y);

  end
