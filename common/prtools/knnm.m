%KNNM Trainable K-Nearest Neighbour density estimation
%
%   W = KNNM(A,KNN)
%   W = A*KNNM([],KNN)
%   W = A*KNNM(KNN)
%
%   D = B*W
% 
% INPUT
% 	A     Dataset used for training
%   B     Dataset used for evaluation
%   KNN   Number of nearest neighbours
%
% OUTPUT
%   W     Density estimate 
%
% DESCRIPTION  
% A density estimator is constructed based on the k-Nearest Neighbour rule
% using the objects in A. In case A is labeled, density estimates are
% performed classwise and combined by the class priors. The default KNN is 
% the square root of the size of the class. The data is scaled by variance
% normalisation determined by the training set. 
%
% The mapping W may be applied to a new dataset B using DENSITY = B*W.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, KNNC, PARZENM, GAUSSM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = knnm(varargin)

  mapname = 'KNN density estimation';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],[]);
  
  if mapping_task(argin,'definition')
    
    w = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.

		[a,knn] = deal(argin{:});
		if isa(a,'prdataset')
			labname = getname(a);
		else
			labname = '';
		end
		if (~isdataset(a) && ~isdatafile(a)) || getsize(a,3) ~= 1
			w = mclassm(a,prmapping(mfilename,knn),'weight');
			w = setlabels(w,labname);
			w = setname(w,mapname);
			return
		end
		islabtype(a,'crisp');
		isvaldfile(a,1);
		a = testdatasize(a,'objects');
		a = remclass(a);
		[m,k] = size(a);
    knn = setdefaults({knn},round(sqrt(m)));
		w = prmapping(mfilename,'trained',{a,knn},labname,k,1);
		w = setname(w,mapname);

	else															% Execute a trained mapping V.

    [a,v] = deal(argin{:});
    [b,knn] = getdata(v);
		[m,k] = size(b); 
		u = scalem(b,'variance');
		a = a*u;
		b = b*u;
		d = sqrt(distm(+a,+b));					% Calculate squared distances
		s = sort(d,2);						    	% and find nearest neighbours.
		ss = s(:,knn);									% Find furthest neighbour
		ss(ss==0) = realmin.^(1/k);	    % Avoid zero distances
		f = knn./(m*nsphere(k)*(ss.^k));% Normalize by the volume of sphere
																		% defined by the furthest neighbour.
		f = f*prod(u.data.rot);         % Correct for scaling																

		w = setdat(a,f,v);							% proper normalisation has still to be done.

	end

return

function v = nsphere(k)

	% Volume k-dimensional sphere
	v = (2*(pi^(k/2))) / (k*gamma(k/2));
	
return
