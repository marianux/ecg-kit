%CDATS Support routine for checking datasets
%
% [B,C,LABLIST,P] = CDATS(A,REDUCE)
%
% INPUT
%   A        Dataset or double
%   REDUCE   0/1, Reduce A to labeled samples (1) or not (0, default), optional.
%
% OUTPUT
%   B        Dataset
%   C        Number of classes
%   LABLIST  Label list of A
%   P        Priors
%
% DESCRIPTION
% This routine supports dataset checking and conversion for mappings
% and density estimators. If A is double it is converted to a one-class
% dataset with all labels set to 1. The same holds if A is an entirely
% unlabeled dataset. The label list of A is returned in LABLIST, but
% for multi-class datasets LABLIST is set to 1. C is the number of classes
% in A. The priors are returned in P (length C), the dataset itself in B.
% If REDUCE is 1, all unlabeled objects are removed.
%
% If A has soft labels, B = A.
% If A does not have labels but targets: C = 1, LABLIST = [], P = 1.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

function [a,c,lablist,p] = cdats(a,red)
	
	if (nargin < 2), red = 0; end

	if (isdataset(a) & islabtype(a,'targets'))
		c = 1; lablist = []; p = 1;
	elseif (isdataset(a) & islabtype(a,'soft'))
		c = getsize(a,3); 
		lablist = getlablist(a);
		if nargout > 3, p = getprior(a); end
	else
		if ~isvaldfile(a)
			a = prdataset(a);
      % we decided to allow doubles instead of datasets where possible.
      % However, if the user calls cdats we assume he needs a dataset
      % definitely.
		end
		c = getsize(a,3); 
		if nargout > 3 % avoid unnecessary warnings
			p = getprior(a);
		end
		if (c == 0)
			a = setlabels(a,1);
			p = 1;
			c = 1;
			lablist = 1;
			prwarning(4,'Dataset unlabeled: all objects will be used')
		elseif (c == 1)
			p = 1;
			if (red == 1)
				a = seldat(a); % remove unlabeled objects
			end
			lablist = getlablist(a);
		else
			prwarning(4,'Multiclass dataset will be combined using priors')
			if (red == 1)
				a = seldat(a); % remove unlabeled objects
			end
			lablist = 1;
		end
		isvaldfile(a,1,1); % at least one object per class
	end;
return