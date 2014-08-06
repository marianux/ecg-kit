%PARZENM Trainable mapping for estimating Parzen densities
%
%   W = PARZENM(A,H)
%   W = A*PARZENM([],H)
%   W = A*PARZENM(H)
%
%   D = B*W
%
% INPUT
%   A  Input dataset
%   H  Smoothing parameters (scalar, vector)
%
% OUTPUT
%   W  output mapping
%
% DESCRIPTION
% A Parzen distribution is estimated for the labeled objects in A. Unlabeled
% objects are neglected, unless A is entirely unlabeled or double. Then all
% objects are used. If A is a multi-class dataset the densities are estimated
% class by class and then weighted and combined according their prior
% probabilities. In all cases, just single density estimator W is computed.
% 
% The mapping W may be applied to a new dataset B using DENSITY = B*W.
%
% The smoothing parameter H is estimated by PARZENML if not supplied. It
% can be a scalar or a vector with as many components as A has features.
%
% EXAMPLE
% See PREX_DENSITY for densities and PREX_PARZEN for differences between
% PARZENC, PARZENDC and PARZENM.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, KNNM, GAUSSM, PARZENML, PARZENDC, PARZENC, KNNM,
% PREX_PARZEN

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = parzenm(varargin)
  
	mapname = 'Parzen Density';
  argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[],1);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,h,n] = deal(argin{:});
    a = prdataset(a); a = remclass(a);
    labname = getname(a);
    islabtype(a,'crisp','soft');
    isvaldfile(a,2,1); % at least 2 objects per class, 1 class
    a = testdatasize(a);
    a = testdatasize(a,'objects');

    if (getsize(a,3) ~= 1)
      w = mclassm(a,prmapping(mfilename,h),'weight');
      w = setlabels(w,labname);
      w = setname(w,mapname);
      return
    end

    [m,k] = size(a);

    % Scale A such that its mean is shifted to the origin and 
    % the variances of all features are scaled to 1. 

    if isempty(h) % if no smoothing parameter given, we have to estimate
                  % it later, lets scale first
      ws = scalem(a,'variance');
    else
      ws = affine(ones(1,k),zeros(1,k),a);
    end
    b = a*ws;

    % SCALE is basically [1/mean(A) 1/STD(A)] based on the properties of SCALEM.
    scale = ws.data.rot;				
    if (size(scale,1) ~= 1) % formally ws.data.rot stores a rotation matrix 
      scale = diag(scale)'; % extract the diagonal if it does,
    end                     % otherwise we already have the diagonal
    if isempty(h)
      if n==1
        h = repmat(parzenml(b),1,k)./scale;
      else
        h = repmat(emparzenml(b,n),1,k)./repmat(scale,n,1);
      end
    end
    w = prmapping('parzen_map','trained',{a,h},labname,k,1);
    w = setname(w,mapname);
    
  end
return
