%FEATSELB Backward feature selection for classification
% 
%  [W,R] = FEATSELB(A,CRIT,K,T,FID)
%  [W,R] = FEATSELB(A,CRIT,K,N,FID)
%
% INPUT	
%   A     Dataset
%   CRIT  String name of the criterion or untrained mapping 
%         (optional; default: 'NN', i.e. 1-Nearest Neighbor error)
%   K     Number of features to select 
%         (optional; default: return optimally ordered set of all features)
%   T     Tuning set (optional)
%   N     Number of cross-validations
%   FID   File ID to write progress to (default [], see PRPROGRESS)
%
% OUTPUT
%   W     Output feature selection mapping
%   R     Matrix with step-by-step results of the selection
%
% DESCRIPTION
% Backward selection of K features using the dataset A. CRIT sets the 
% criterion used by the feature evaluation routine FEATEVAL. If the  
% dataset T is given, it is used as test set for FEATEVAL. Alternatvely a
% a number of cross-validation N may be supplied. For K = 0, the optimal 
% feature set (corresponding to the maximum value of FEATEVAL) is returned. 
% The result W can be used for selecting features by B*W. In this case, 
% features are ranked optimally. 
% The selected features are stored in W.DATA and can be found by +W.
% In R, the search is reported step by step as:
% 
% 	R(:,1) : number of features
% 	R(:,2) : criterion value
% 	R(:,3) : added / deleted feature
%
% SEE ALSO
% MAPPINGS, DATASETS, FEATEVAL, FEATSELLR, FEATSEL,
% FEATSELO, FEATSELF, FEATSELI, FEATSELP, FEATSELM, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: featselb.m,v 1.6 2008/07/03 09:08:43 duin Exp $

function [w,r] = featselb(a,crit,ksel,t,fid)

		
	if (nargin < 2) | isempty(crit)
		prwarning(2,'No criterion specified, assuming 1-NN.');
		crit = 'NN';
	end
	if (nargin < 3) | isempty(ksel)
		ksel = 0; 		% Consider all the features and sort them.
	end
	if (nargin < 4)
		prwarning(3,'No tuning set supplied.');
		t = [];
	end
	if (nargin < 5)
		fid = [];
	end
	
	if nargin == 0 | isempty(a)
		% Create an empty mapping:
		w = prmapping(mfilename,{crit,ksel,t});
	else
		prprogress(fid,'\nfeatselb : Backward Feature Selection')
		[w,r] = featsellr(a,crit,ksel,0,1,t,fid);
		%DXD This is a patch: when the number of features has to be
		%optimized, and all features seem useful, when the list of
		%features is not reshuffled to reflect the relative importance of
		%the features:
		% (Obviously, this should be fixed in featsellr, but I don't
		% understand what is happening in there)
		dim = size(a,2);
		if (ksel==0) & (length(getdata(w))==dim)
			rr = -r(:,3); rr(1) = [];
			rr = [setdiff((1:dim)',rr) rr(end:-1:1)'];
			w = setdata(w,rr);
		end
		prprogress(fid,'featselb  finished\n')
	end
	w = setname(w,'Backward FeatSel');

return
