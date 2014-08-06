%DISTMAHA Mahalanobis distance
% 
% 	D = DISTMAHA (A,U,G)
% 
% INPUT
%		A		Dataset
%		U		Mean(s) (optional; default: estimate on classes in A)
%		G		Covariance(s) (optional; default: estimate on classes in A)
%
% OUTPUT
%		D		Mahalanobis distance matrix
%
% DESCRIPTION
% Computes the M*N Mahanalobis distance matrix of all vectors in M*K dataset
% A to an N*K dataset of points U, using the covariance matrix or matrices
% G. G should be either be one K*K matrix, or a K*K*N matrix containing a
% covariance matrix for each point in U.
%
% When U and G are not specified, it estimates the C*C Mahalanobis distance
% matrix between all classes in A: the distance between the means,
% relative to the average per-class covariance matrix.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DISTM, PROXM, MEANCOV

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: distmaha.m,v 1.4 2010/02/08 15:31:48 duin Exp $

function d = distmaha (a,u,g)

		[m,k,c] = getsize(a);
	a = testdatasize(a);
	a = testdatasize(a,'features');

	if (nargin == 1)    

		% Calculate Mahalanobis distance matrix between classes in dataset:
		% the distance matrix between the class means, sphered using the
		% average covariance matrix of the per-class centered data.

		u = zeros(c,k);
		for i = 1:c
			J = findnlab(a,i);
			if (~isempty(J))

				% Mean of class I.
				u(i,:) = mean(a(J,:),1);									

				% Remove class mean, for later covariance calculation.
				a(J,:) = a(J,:) - repmat(+u(i,:),length(J),1); 

			end
		end 

		% Sphere means and calculate distance matrix.
		[E,V] = preig(covm(a)); u = u*E*sqrt(prinv(V));
		d = distm(u);															

	elseif (nargin == 2)
		
		% Calculate Euclidean distance between data and means.

		prwarning(4,'covariance matrix not specified, assuming G = I');

		d = distm(a,u);

	elseif (nargin == 3)     

		% Calculate distance between A and distributions specified by U and G.

      % Check arguments.
      [kg1,kg2,cg] = size(g); [cu,ku] = size(u);
      if any([kg1,kg2,ku] ~= k) | (cu ~= cg & cg ~= 1)
         error('sizes of mean(s) and covariance(s) do not match')
      end	

		% Get labels of means if present, or assign labels 1:N. These are
		% only used in the construction of the dataset D.
      if (isa(u,'prdataset')), labels = getlab(u); else, labels = [1:cu]'; end

		% If there is only one covariance matrix, invert it here.
      if (cg == 1), g_inv = prinv(g); end

      d = zeros(m,cu);
      for j = 1:cu

			% If each mean has its own covariance matrix, invert it here.
         if (cg ~= 1), g_inv = prinv(g(:,:,j)); end

			% And calculate the standard Mahalanobis distance.
         d(:,j) = sum((a-repmat(+u(j,:),m,1))'.*(g_inv*(a-repmat(+u(j,:),m,1))'),1)';

      end

		% Construct the distance matrix dataset.
		d = setdata(a,d,labels);
		d = setname(d,'Maha Dist');

  else

      error('incorrect number of arguments')

  end

return
