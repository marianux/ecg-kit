%PCLDC Linear classifier using PC expansion on the joint data.
% 
% 	W = PCLDC(A,N)
% 	W = PCLDC(A,ALF)
%
% INPUT
%  A    Dataset
%  N    Number of eigenvectors
%  ALF  Total explained variance (default: ALF = 0.9)
%
% OUTPUT
%  W    Mapping
% 
% DESCRIPTION
% Finds the linear discriminant function W for the dataset A 
% computing the LDC on a projection of the data on the first N  
% eigenvectors of the total dataset (Principle Component Analysis).
% 
% When ALF is supplied the number of eigenvalues is chosen such that at 
% least a part ALF of the total variance is explained. 
% 
% If N (ALF) is NaN it is optimised by REGOPTC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, KLLDC, KLM, FISHERM, REGOPTC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: pcldc.m,v 1.4 2007/06/13 21:59:42 duin Exp $

function W = pcldc(a,n)

			if nargin < 2, n = []; end
	
	if nargin == 0 | isempty(a)
		W = prmapping('pcldc',{n});
		
	elseif isnan(n)    % optimize regularisation parameter
		defs = {1};
		parmin_max = [1,size(a,2)];
		W = regoptc(a,mfilename,{n},defs,[1],parmin_max,testc([],'soft'),0);
		
	else

		islabtype(a,'crisp','soft');
		isvaldfile(a,2,2); % at least 2 object per class, 2 classes
		a = testdatasize(a,'features');
		a = setprior(a,getprior(a));

		% Make a sequential classifier combining PCA and LDC:
		v = pcam(a,n);
		W = v*ldc(a*v);
		W = setcost(W,a);
		
	end
	
	W = setname(W,'PC Bayes-Normal-1');

return
