% [X,centering,scaling] = pls_prepro(X,centering,scaling, flag)
function [X,centering,scaling] = pls_prepro(X,centering,scaling, flag)

% Copyright: S.Verzakov, serguei@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

if nargin<4
  flag = 1;
end

[N,d] = size(X);

centering = centering(:).';
scaling = scaling(:).';

if flag >= 0
  if length(centering) == 1
    if isnan(centering)
      centering = mean(X,1);
      X = X - repmat(centering, [N,1]);
    else
      X = X - centering;
    end

  elseif length(centering) == d
    idx = find(isnan(centering));
    centering(idx) = mean(X(:,idx),1);
    X = X - repmat(centering, [N,1]);
  end

  if length(scaling) == 1
    if isnan(scaling)
      scaling = std(X,0,1);
      idx0 = find(scaling == 0);
      scaling(idx0) = 1;
      warning(['features ' num2str(idx(:)) ' have std = 0 and are not scaled']);
      X = X ./ repmat(scaling, [N,1]);
    else
      X = X / scaling;
    end

  elseif length(scaling) == d
    idx = find(isnan(scaling));
    scaling(idx) = std(X(:,idx),0,1);
    idx0 = find(scaling(idx) == 0);
    scaling(idx(idx0)) = 1;  
    warning(['features ' num2str(idx(idx0(:))) ' have std = 0 and are not scaled']);
    X = X ./ repmat(scaling, [N,1]);
  end

else
  if length(centering) > 0 & all(~isnan(centering))
    if length(centering) == 1 
      X = X + centering;
    elseif length(centering) == d 
      X = X + repmat(centering, [N,1]);
    end
  end

  if length(scaling) > 0 & all(~isnan(scaling))
    if length(scaling) == 1
      X = X * scaling;
    elseif length(scaling) == d
      X = X .* repmat(scaling, [N,1]);
    end
  end
end

return

