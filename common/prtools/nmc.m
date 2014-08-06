%NMC Nearest Mean Classifier
% 
%   W = NMC(A)
%   W = A*NMC
%
% INPUT
%   A    Dataset
%
% OUTPUT
%   W    Nearest Mean Classifier   
%
% DESCRIPTION
% Computation of the nearest mean classifier between the classes in the
% dataset A.  The use of soft labels is supported. Prior probabilities are
% not used.
%
% The difference with NMSC is that NMSC is based on an assumption of normal
% distributions and thereby automatically scales the features and is
% sensitive to class priors. NMC is a plain nearest mean classifier for
% which the assigned classes are are sensitive to feature scaling and 
% unsensitive to class priors.
%
% The estimated class confidences by B*NMC(A), however, are based on
% assumed spherical gaussian distributions of the same size.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, NMSC, LDC, FISHERC, QDC, UDC 

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: nmc.m,v 1.14 2009/12/09 15:49:54 duin Exp $

function W = nmc(varargin)
  
	mapname = 'NearestMean';
  argin = setdefaults(varargin,[],false);
  
  if mapping_task(argin,'definition')
    W = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,flag] = deal(argin{:});
    % flag = true forces execution by knn_map instead of normal_map
    % flag = false (defaults) finds posteriors by spherical normal
    % distributions

    if nargin < 2, flag = false; end
    % flag=true forces to skip computation of posteriors based on normal
    % densities


    if nargin < 1 | isempty(a)
      W = prmapping(mfilename);
      W = setname(W,'NearestMean');
      return
    end

    islabtype(a,'crisp','soft');
    isvaldfile(a,1,2); % at least 1 object per class, 2 classes

    [m,k,c] = getsize(a);
    if isdatafile(a), a = setfeatsize(a,k); end

    if ~flag & c == 2      % 2-class case: store linear classifier
      u = meancov(a);
      u1 = +u(1,:);
      u2 = +u(2,:);
      R = [u1-u2]';
      offset =(u2*u2' - u1*u1')/2;
      W = affine([R -R],[offset -offset],a,getlablist(a));
      W = cnormc(W,a);
      W = setname(W,'Nearest Mean');
    else 
      if ~flag & all(classsizes(a) > 1)
        a = setprior(a,0);  % NMC should be independent of priors: make them equal
        p = getprior(a);
        U = zeros(c,k);
        V = zeros(c,k);
        s = sprintf('NMC: Processing %i classes: ',c);
        prwaitbar(c,s);
        for j=1:c
          prwaitbar(c,j,[s int2str(j)]);
          b = seldat(a,j);
          [v,u] = var(b);
          U(j,:) = +u;
          V(j,:) = +v;
        end
        prwaitbar(0);
        %G = mean(V'*p') * eye(k);
        G = mean(V(:));
        w.mean = +U;
        w.cov = G;
        w.prior =p;
        W = normal_map(w,getlablist(a),k,c);
        W = setname(W,'Nearest Mean');
        W = setcost(W,a);
      else
        u = meancov(a);
        W = knnc(u,1);
        W = setname(W,mapname);
      end

    end

    W = setcost(W,a);  
    
  end
	
	return

