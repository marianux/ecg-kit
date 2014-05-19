%RBSVC Automatic radial basis Support Vector Classifier
%
%   [W,KERNEL,NU] = RBSVC(A)
%
% INPUT
%   A	      Dataset
%
% OUTPUT
%   W       Mapping: Radial Basis Support Vector Classifier
%   KERNEL  Untrained mapping, representing the optimised kernel
%   NU      Resulting value for NU from NUSVC
%
% DESCRIPTION
% This routine computes a classifier by NUSVC using a radial basis kernel
% with an optimised standard deviation by REGOPTC. The resulting classifier
% W is identical to NUSVC(A,KERNEL,NU). As the kernel optimisation is based
% on internal cross-validation the dataset A should be sufficiently large.
% Moreover it is very time-consuming as the kernel optimisation needs
% about 100 calls to SVC.
%
% SEE ALSO
% MAPPINGS, DATASETS, PROXM, SVC, NUSVC, REGOPTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,kernel,nu] = rbsvc(a,sig)

if nargin < 2 | isempty(sig)
	sig = NaN;
end

if nargin < 1 | isempty(a)
	
	w = prmapping(mfilename,{sig});
	
else
	
	islabtype(a,'crisp');
	isvaldfile(a,2,2);              % at least 1 object per class, 2 classes
	a = testdatasize(a,'objects');
  c = getsize(a,3);
  
	if (c > 2)
		
    % Compute c classifiers: each class against all others.	
		w = mclassc(a,prmapping(mfilename,{sig}));	 
    
  else
	
	  if isnan(sig) % optimise sigma
      
		  % find upper bound
		  d = sqrt(+distm(a));
		  sigmax = min(max(d)); % max: smallest furthest neighbor distance
      % find lower bound
		  d = d + 1e100*eye(size(a,1));
		  sigmin = max(min(d)); % min: largest nearest neighbor distance
		  % call optimiser
		  defs = {1};
		  parmin_max = [sigmin,sigmax];
		  [w,kernel,nu] = regoptc(a,mfilename,{sig},defs,[1],parmin_max,testc([],'soft'));
		
	  else % kernel is given
		
		  kernel = proxm([],'r',sig);
		  [w,J,nu] = nusvc(a,kernel);

	  end
    
  end
	
end

w = setname(w,'RB-SVM');
return
