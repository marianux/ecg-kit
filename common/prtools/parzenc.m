%PARZENC Optimisation of the Parzen classifier
% 
%  [W,H] = PARZENC(A,H)
%  [W,H] = A*PARZENC([],H)
%  [W,H] = A*PARZENC(H)
% 
% INPUT
%  A    dataset
%  H    smoothing parameter (may be scalar, vector of per-class
%       parameters, or matrix with parameters for each class (rows) and
%       dimension (columns))
%
% OUTPUT
%  W    trained mapping
%  H    estimated smoothing (scalar value)
%
% DESCRIPTION
% Computation of the optimum smoothing parameter H for the Parzen 
% classifier between the classes in the dataset A. The leave-one-out 
% Lissack & Fu estimate is used for the classification error E. The 
% final classifier is stored as a mapping in W. It may be converted
% into a classifier by W*CLASSC. PARZENC cannot be used for density
% estimation.
% 
% In case smoothing H is specified, no learning is performed, just the
% discriminant W is produced for the given smoothing parameters H.
% Smoothing parameters may be scalar, vector of per-class parameters, or 
% a matrix with individual smoothing for each class (rows) and feature
% directions (columns)
%
% EXAMPLES
% See PREX_DENSITY for densities and PREX_PARZEN for differences between
% PARZENC, PARZENDC and PARZENM.
%
% REFERENCES
% T. Lissack and K.S. Fu, Error estimation in pattern recognition via
% L-distance between posterior density functions, IEEE Trans. Inform. 
% Theory, vol. 22, pp. 34-45, 1976.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PARZEN_MAP, PARZENML, PARZENDC, PARZENM, CLASSC, 
% PREX_PARZEN 
 
% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: parzenc.m,v 1.6 2008/07/03 09:11:44 duin Exp $

function [W,h] = parzenc(varargin)
  
	mapname = 'ParzenC';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],[]);
  
  if mapping_task(argin,'definition')
    W = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,h] = deal(argin{:});
    islabtype(a,'crisp','soft');
    isvaldfile(a,2,2); % at least 2 objects per class, 2 classes
    a = testdatasize(a);
    a = testdatasize(a,'objects');

    [m,k,c] = getsize(a);
    nlab = getnlab(a);

    if ~isempty(h)       % take user setting for smoothing parameter

      if size(h,1) == 1, h = repmat(h,c,1); end
      if size(h,2) == 1, h = repmat(h,1,k); end
      if any(size(h) ~= [c,k])
        error('Array with smoothing parameters has wrong size');
      end
      W = prmapping('parzen_map','trained',{a,h},getlablist(a),k,c);
      W = setname(W,'Parzen Classifier');
      return

    end

    % compute all object distances
    % make diagonal inf to exclude objects own contribution
    D = +distm(a) + diag(inf*ones(1,m));

    % find object frequencies
    if islabtype(a,'crisp')
      csize = classsizes(a);
      of = csize(nlab);
    else
      csize = sum(gettargets(a),1);
    end

    % find object weights q
    p = getprior(a);
    a = setprior(a,p);
    q = p(nlab)./csize(nlab);

    % initialise
    h = max(std(a)); % for sure a too high value
    L = -inf;
    Ln = 0;
    z = 0.1^(1/k); % initial step size

    % iterate

    iter = 0;
    prwaitbar(100,'parzenc: Optimizing smoothing parameter',m > 100);
    while abs(Ln-L) > 0.001 & z < 1

      % In L we store the best performance estimate found so far.
      % Ln is the actual performance (for the actual h)
      % If Ln > L we improve the bound L, and so we rest it.

      if Ln > L, L = Ln; end
      iter = iter+1;
      prwaitbar(100,100-100*exp(-iter/10));

      r = -0.5/(h^2);
      F = q(ones(1,m),:)'.*exp(D*r);           % density contributions
      FS = sum(F)*((m-1)/m); IFS = find(FS>0); % joint density distribution
      if islabtype(a,'crisp');
        G = sum(F .* (nlab(:,ones(1,m)) == nlab(:,ones(1,m))'));
        G = G.*(of-1)./of;                     % true-class densities
      else
        % here we are for soft labels (stored in targets)
        G = zeros(1,m);
        for j=1:c
          G = G + sum(F .* (a.targets(:,j) * a.targets(:,j)'));
          % to be corrected for bias?
        end
      end

      % performance estimate
      en = max(p)*ones(1,m);
      en(IFS) = (G(IFS))./FS(IFS);
      Ln = exp(sum(log(en))/m);

      if Ln < L            % compute next estimate
        z = sqrt(z);       % adjust stepsize up (recall: 0 < z < 1)
        h = h / z;         % if we don't improve, increase h (approach h_opt from below)
      else
        h = h * z;         % if we improve, decrease h (approach h_opt from above)
      end
    end
    prwaitbar(0);
    W = prmapping('parzen_map','trained',{a,repmat(h,c,k);},getlablist(a),k,c);
    W = setname(W,mapname);
    W = setcost(W,a);
    
  end

return
