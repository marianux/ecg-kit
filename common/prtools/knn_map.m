%KNN_MAP Map a dataset on a K-NN classifier
%
%   F = KNN_MAP(A,W)
% 
% INPUT
%   A  Dataset
%   W  K-NN classifier trained by KNNC
%
% OUTPUT
%   F  Posterior probabilities
%
% DESCRIPTION  
% Maps the dataset A by the K-NN classifier W on the [0,1] interval for 
% each of the classes that W is trained on. The posterior probabilities, 
% stored in F, are computed in the following ways:
% soft labeled training set: the normalised average of the soft labels 
%              of the K neighbors.
% crisp labeled training set, K = 1: normalisation of sigm(log(F)) with 
%              F(1:C) = sum(NN_Dist(1:C))./NN_Dist(1:C) - 1
%              in which C is the number of classes and NN_Dist stores
%              the distance to the nearest neighbor of each class.
% crisp labeled training set, K > 1: normalisation of
%              (N(1:C) + 1)/(K+C), in which N stores the number of
%              objects per class within the K first neighbors.
%
% This routine is called automatically to determine A*W if W is trained 
% by KNNC.
%
% Warning: Class prior probabilities in the dataset A are neglected.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, KNNC, TESTK

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: knn_map.m,v 1.3 2007/06/19 11:44:14 duin Exp $

function F = knn_map(T,W)
		% Get the training data and parameters from the mapping:
	data = getdata(W);
	a = data{1};
	knn = data{2};
	[m,k,c] = getsize(a);
	nlab = getnlab(a);

	% If there is no test set, then the leave-one-out is done on the
	% training set (see TESTK).
	if isempty(T) 
		T = a;
		loo = 1;
	else
		loo = 0;
	end
	[mt,kt] = size(T);
	if (kt ~= k), error('Wrong feature size'); end

	r = classsizes(a);
	[num,n] = prmem(mt,m);				% Check the available memory. 
	F = ones(mt,c);
	D = ones(mt,c);
  s = sprintf('Compute distance matrix in %i batches: ',num);
  prwaitbar(num,s);

	% Loop in batches. 
  if isdatafile(a)
    a = prdataset(a);
  end
	for i = 0:num-1
    prwaitbar(num,i+1,[s int2str(i+1)]);									
		if (i == num-1)
			nn = mt - num*n + n;
		else
			nn = n;
		end
		range = [i*n+1:i*n+nn];
		if loo,
			DD = +distm(a,a(range,:));
			dmax = max(DD(:));
			% Set distances to itself at INF to find later the nearest
			% neighbors more easily
			DD(i*n+1:m+1:i*n+nn*m) = inf*ones(1,nn); 
    else
      if isdatafile(T)
        t = +prdataset(T(range));
      else
        t = +T;
      end
			DD = distm(+a,t);
			dmax = max(DD(:));
		end
		J = find(isnan(DD));
		if length(J) > 0
			DD(J) = dmax*10;
		end
		[DD,L] = sort(DD);     	
		
    ss = sprintf('Prepare %i classes: ',c);
    prwaitbar(c,ss);
		switch getlabtype(a)
			
			case 'soft'
				for j=1:c
          prwaitbar(c,j,[ss num2str(j)]);
					F(range,j) = sum(reshape(a.targets(L(1:knn,:),j),knn,length(range)),1)';
				end
				
			case 'crisp'
				L = reshape(nlab(L),size(L));				% Find labels.
				% Find label frequencies.
				for j = 1:c									
          prwaitbar(c,j,[ss num2str(j)]);				
					F(range,j) = sum(L(1:knn,:)==j,1)';
				end
				
			otherwise
				error('Illegal label type')
    end
    prwaitbar(0);

		% Estimate posterior probabilities
		if islabtype(a,'crisp') 
			if (knn >= 2)           % Use Bayes estimators on frequencies.
				F(range,:) = (F(range,:)+1)/(knn+c);
      else			              % Use distances for estimating posteriors
				K = max(F(range,:)');
        ss = sprintf('Posteriors for %i classes: ',c);
        prwaitbar(c,ss);
				for j = 1:c
          prwaitbar(c,j,[ss num2str(j)]);
					K = min(K,r(j));
					J = reshape(find(L==j),r(j),nn);	% Find the distances between 
					J = J(K+[0:nn-1]*r(j));						% that neighbor and other objects.
					D(range,j) = DD(J)';							% Number for all classes.
				end
				F(range,:) = sigm(log(sum(D(range,:),2)*ones(1,c)./ ...
									 (D(range,:)+realmin) - 1 + realmin));
        prwaitbar(0);
			end
		end
		% Normalize the probabilities.
		F(range,:) = F(range,:) ./ (sum(F(range,:),2)*ones(1,c));
	end
  prwaitbar(0);

	if (isdataset(T))
		F = setdata(T,F,getlabels(W));
	end;

return;

