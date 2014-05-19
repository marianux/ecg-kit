%FEATSELF Forward feature selection for classification
% 
% [W,R] = FEATSELF(A,CRIT,K,T,FID)
% [W,R] = FEATSELF(A,CRIT,K,N,FID)
%
% INPUT	
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping 
%        (default: 'NN', i.e. the 1-Nearest Neighbor error)
%   K    Number of features to select (default: K = 0, return optimal set)
%   T    Tuning dataset (optional)
%   N    Number of cross-validations (optional)
%   FID  File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W    Output feature selection mapping
%   R    Matrix with step-by-step results
%
% DESCRIPTION
% Forward selection of K features using the dataset A. CRIT sets the
% criterion used by the feature evaluation routine FEATEVAL. If the  
% dataset T is given, it is used as test set for FEATEVAL. Alternatvely a
% a number of cross-validation N may be supplied. For K = 0, the optimal 
% feature set (corresponding to the maximum value of FEATEVAL) is returned.  
% The result W can be used for selecting features using B*W.  
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
% 
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
% 
% SEE ALSO 
% MAPPINGS, DATASETS, FEATEVAL, FEATSELLR, FEATSEL,
% FEATSELO, FEATSELB, FEATSELI, FEATSELP, FEATSELM, PRPROGRESS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: featself.m,v 1.3 2007/04/21 23:06:46 duin Exp $

function [w,r] = featself(a,crit,ksel,t,fid)
		 
  	
  if (nargin < 2) | isempty(crit)
    prwarning(2,'no criterion specified, assuming NN');
    crit = 'NN';
  end
  if (nargin < 3) | isempty(ksel)
    ksel = 0;
  end
  if (nargin < 4)
    prwarning(3,'no tuning set supplied (risk of overfit)');
    t = [];
  end
	if (nargin < 5)
		fid = [];
	end
	
	if nargin == 0 | isempty(a)
		% Create an empty mapping:
		w = prmapping(mfilename,{crit,ksel,t});
	else
		prprogress(fid,'\nfeatself : Forward Feature Selection')
		[w,r] = featsellr(a,crit,ksel,1,0,t,fid);
		prprogress(fid,'featself  finished\n')
	end
	w = setname(w,'Forward FeatSel');

return
