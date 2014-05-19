%KNNM K-Nearest Neighbour based density estimate
%
%   W = KNNM(A,KNN)
%
%   D = B*W
% 
% INPUT
% 	A     Dataset
%   KNN   Number of nearest neighbours	
%
% OUTPUT
%   W     Density estimate 
%
% DESCRIPTION  
% A density estimator is constructed based on the k-Nearest Neighbour rule
% using the labeled objects in A. All objects, however, are used if the
% entire dataset is unlabeled or if A is not a dataset, but a double. 
% In all cases, just single density estimator W is computed.
% Class priors in A are neglected.
% 
% The mapping W may be applied to a new dataset B using DENSITY = B*W.
%
% SEE ALSO
% DATASETS, MAPPINGS, KNNC, QDC, PARZENM, GAUSSM

% $Id: knnm.m,v 1.6 2009/02/02 22:11:27 duin Exp $

function w = knnm(a,arg2)

		if (nargin < 2)
		prwarning(2,'number of neighbours to use not specified, assuming 1');
		arg2 = 1; 
	end

  mapname = 'KNN Density Estimation';
	if (nargin < 1) | (isempty(a)) 		% No arguments given: return empty mapping.

		w = prmapping(mfilename,arg2);
		w = setname(w,mapname);

	elseif (~isa(arg2,'prmapping'))			% Train a mapping.

		knn = arg2;
		
		if isa(a,'prdataset')
			labname = getname(a);
		else
			labname = '';
		end
		if (~isdataset(a) & ~isdatafile(a)) | getsize(a,3) ~= 1
			w = mclassm(a,prmapping(mfilename,knn),'weight');
			w = setlabels(w,labname);
			w = setname(w,mapname);
		return	
			return
		end
		islabtype(a,'crisp');
		isvaldfile(a,1);
		a = testdatasize(a,'objects');
		a = remclass(a);
		[m,k] = size(a);

		w = prmapping(mfilename,'trained',{a,knn},labname,k,1);
		w = setname(w,mapname);

	else															% Execute a trained mapping V.

		v = arg2;
		b   = v.data{1};								% B is the original training data.
		knn = v.data{2};
		%iscomdset(a,b,0);
		[m,k] = size(b); 
		u = scalem(b,'variance');
		a = a*u;
		b = b*u;
		d = sqrt(distm(+a,+b));					% Calculate squared distances
		[s,J] = sort(d,2);							%   and find nearest neighbours.
		ss = s(:,knn);									% Find furthest neighbour
		ss(find(ss==0)) = realmin.^(1/k);	% Avoid zero distances
		f = knn./(m*nsphere(k)*(ss.^k)); 	% Normalize by the volume of sphere
																		%   defined by the furthest neighbour.
		f = f*prod(u.data.rot);         % Correct for scaling																

		w = setdat(a,f,v);							% proper normalisation has still to be done.

	end

return

function v = nsphere(k)

	% Volume k-dimensional sphere
	v = (2*(pi^(k/2))) / (k*gamma(k/2));
	
return
