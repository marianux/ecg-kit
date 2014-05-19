%SVC_NU Support Vector Classifier: NU algorithm
% 
%  This routine is outdated, use NUSVC instead
%
% 	[W,J,C] = SVC(A,TYPE,PAR,NU,MC,PD)
%
% INPUT
%   A	    Dataset
%   TYPE  Type of the kernel (optional; default: 'p')
%   PAR   Kernel parameter (optional; default: 1)
%   NU    Regularization parameter (0 < NU < 1): expected fraction of SV 
%         (optional; default: max(leave-one-out 1_NN error,0.05))
%
%   MC    Do or do not data mean-centering (optional; default: 1 (to do))
%   PD    Do or do not the check of the positive definiteness (optional; default: 1 (to do))
%
% OUTPUT
%   W     Mapping: Support Vector Classifier
%   J     Object identifiers of support objects		
%   C     Equivalent C regularization parameter of SVM-C algorithm 
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by 
% quadratic programming. The classifier can be of one of the types 
% as defined by PROXM. Default is linear (TYPE = 'p', PAR = 1). In J 
% the identifiers of the support objects in A are returned.
%
% NU belongs to the interval (0,1). NU close to 1 allows for more class overlap.
% Default NU = 0.25.
% 
% NU is bounded from above by NU_MAX = (1 - ABS(Lp-Lm)/(Lp+Lm)), where
% Lp (Lm) is the number of positive (negative) samples. If NU > NU_MAX is supplied 
% to the routine it will be changed to the NU_MAX.
%
% If NU is less than some NU_MIN which depends on the overlap between the classes 
% the algorithm will typically take a long time to converge (if at all). 
% So, it is advisable to set NU larger than expected overlap.
%
% Output is rescaled in such a manner as if it were returned by SVC with the parameter C.
%
%
% SEE ALSO
% SVO_NU, SVO, SVC, MAPPINGS, DATASETS, PROXM

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Based on SVC.M by D.M.J. Tax, D. de Ridder, R.P.W. Duin
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: svc_nu.m,v 1.6 2009/01/31 18:43:11 duin Exp $

function [W, J, C] = svc_nu(a,type,par,nu,mc,pd)
				warning('SVC_NU is outdated, use NUSVC instead')

  if nargin < 2 | ~isa(type,'prmapping')
    if nargin < 6
      pd = 1;  
    end
    if nargin < 5
      mc = 1;  
    end
    if nargin < 4 | isempty(nu)
      nu = [];
  	  prwarning(3,'Regularization parameter nu set to NN error\n');
    end
    if nargin < 3 | isempty(par)
    	par = 1;
    	prwarning(3,'Kernel parameter par set to 1\n');
    end
    if nargin < 2 | isempty(type)
    	type = 'p';
    	prwarning(3,'Polynomial kernel type is used\n');
    end
    if nargin < 1 | isempty(a)
    	W = prmapping(mfilename,{type,par,nu,mc,pd});
    	W = setname(W,'Support Vector Classifier (nu version)');
    	return;
    end

  	islabtype(a,'crisp');
  	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a,'objects');
  	[m,k,c] = getsize(a);
  	nlab = getnlab(a);  
	
		if isempty(nu), nu = max(testk(a,1),0.01); end

		% The SVC is basically a 2-class classifier. More classes are
  	% handled by mclassc.
  	if c == 2   % two-class classifier
	
  		% Compute the parameters for the optimization:
  		y = 3 - 2*nlab;
      if mc
        u = mean(a);
  	  	a = a -ones(m,1)*u;
      else
        u = [];
			end
      K = a*proxm(a,type,par);
  		% Perform the optimization:
  		[v,J,C] = svo_nu(+K,y,nu,pd);
      % Store the results:
  		W = prmapping(mfilename,'trained',{u,a(J,:),v,type,par},getlablist(a),k,2);
  		%W = cnormc(W,a);
  		W = setname(W,'Support Vector Classifier (nu version)');
  		W = setcost(W,a);
			J = getident(a,J);
  		%J = a.ident(J);
  		
  	else   % multi-class classifier:
	    [W,J,C] = mclassc(a,prmapping(mfilename,{type,par,nu,mc,pd}),'single');
		end

  else % execution
		
		nodatafile(a);
		w = +type;
  	m = size(a,1);
	
	  % The first parameter w{1} stores the mean of the dataset. When it
  	% is supplied, remove it from the dataset to improve the numerical
	  % precision. Then compute the kernel matrix using proxm:
  	if isempty(w{1})
  		d = a*proxm(w{2},w{4},w{5});
  	else
  		d = (a-ones(m,1)*w{1})*proxm(w{2},w{4},w{5});
  	end

  	% Data is mapped by the kernel, now we just have a linear
  	% classifier  w*x+b:

	  d = [d ones(m,1)] * w{3};
    d = sigm([d -d]);
  	W = setdat(a,d,type);
		
	end
	
  return;
