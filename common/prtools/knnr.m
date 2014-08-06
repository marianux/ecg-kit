%KNNR Trainable Nearest Neighbor Regression
%
%    Y = KNNR(X,K)
%    Y = X*KNNR([],K)
%    Y = X*KNNR(K)
%
% INPUT
%   X    Regression dataset, used for training
%   K    number of neighbors (default K=3)
%
% OUTPUT
%   Y    k-nearest neighbor regression
%
% DESCRIPTION
% Define a k-Nearest neighbor regression on dataset X.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% LINEARR, TESTR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = knnr(varargin)

  mapname = 'KNN regression';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],3);
  
  if mapping_task(argin,'definition')
    
    y = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.

		[x,k] = deal(argin{:});
    [n,d] = size(x);
    W.x = +x;
    W.y = gettargets(x);
    W.k = k;
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'k-nearest neighbor regression');
    
  else                                      % Evaluation
    
    [x,v] = deal(argin{1:2});
    w = getdata(v);
    [n,d] = size(x);
    D = distm(+x,w.x);
    [sD,I] = sort(D,2);
    if n==1
      out = mean(w.y(I(:,1:w.k)));
    else
      out = mean(w.y(I(:,1:w.k)),2);
    end
    y = setdat(x,out);

  end
