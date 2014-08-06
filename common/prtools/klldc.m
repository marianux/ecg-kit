%KLLDC Linear classifier built on the KL expansion of the common covariance matrix
% 
%  W = KLLDC(A,N)
%  W = KLLDC(A,ALF)
% 
% INPUT
%  A    Dataset
%  N    Number of significant eigenvectors 
%  ALF  0 < ALF <= 1, percentage of the total variance explained (default: 0.9)
%
% OUTPUT
%  W    Linear classifier 
%
% DESCRIPTION  
% Finds the linear discriminant function W for the dataset A. This is done  
% by computing the LDC on the data projected on the first eigenvectors of
% the averaged covariance matrix of the classes. Either first N eigenvectors
% are used or the number of eigenvectors is determined such that ALF, the 
% percentage of the total variance is explained. (Karhunen Loeve expansion)
%
% If N (ALF) is NaN it is optimised by REGOPTC.
%
%	SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%	MAPPINGS, DATASETS, PCLDC, KLM, FISHERM, REGOPTC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: klldc.m,v 1.4 2007/06/13 21:59:42 duin Exp $

function W = klldc(a,n)
		if nargin < 2
		n = [];
		prwarning(4,'number of significant eigenvectors not supplied, 0.9 variance explained');
	end
	
	if nargin == 0 | isempty(a)
		W = prmapping('klldc',{n});
	
	elseif isnan(n)    % optimize regularisation parameter
		defs = {1};
		parmin_max = [1,size(a,2)];
		W = regoptc(a,mfilename,{n},defs,[1],parmin_max,testc([],'soft'),0);
		
	else
	
		islabtype(a,'crisp','soft');
		isvaldfile(a,2,2); % at least 2 object per class, 2 classes
		a = testdatasize(a,'features');
		a = setprior(a,getprior(a));
		v = klm(a,n);
		W = v*ldc(a*v);
		W = setcost(W,a);
		
	end
	
	W = setname(W,'KL Bayes-Normal-1');
	
return
