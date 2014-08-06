%FEATSELF Trainable mapping for forward feature selection
% 
%   [W,R] = FEATSELF(A,CRIT,K,T)
%   [W,R] = A*FEATSELF([],CRIT,K,T)
%   [W,R] = A*FEATSELF(CRIT,K,T)
%   [W,R] = FEATSELF(A,CRIT,K,N)
%   [W,R] = A*FEATSELF([],CRIT,K,N)
%   [W,R] = A*FEATSELF(CRIT,K,N)
%
% INPUT	
%   A    Training dataset
%   CRIT Name of the criterion or untrained mapping 
%        (default: 'NN', i.e. the LOO 1-Nearest Neighbor error)
%   K    Number of features to select (default: K = 0, return optimal set)
%   T    Tuning dataset (optional)
%   N    Number of cross-validations (optional)
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% MAPPINGS, DATASETS, FEATEVAL, FEATSELLR, FEATSEL,
% FEATSELO, FEATSELB, FEATSELI, FEATSELP, FEATSELM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: featself.m,v 1.3 2007/04/21 23:06:46 duin Exp $

function [w,r] = featself(varargin)
		 
  varargin = shiftargin(varargin,{'char','prmapping'});
  argin = setdefaults(varargin,[],'NN',0,[],[]);
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained','Forward FeatSel');
    return
  end
    
  [a,crit,ksel,t,fid] = deal(argin{:});
  [w,r] = featsellr(a,crit,ksel,1,0,t);

return
