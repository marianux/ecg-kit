%GENDATK K-Nearest neighbor data generation
% 
%   B = GENDATK(A,N,K,S)
%   B = A*GENDATK([],N,K,S)
%   B = A*GENDATK(N,K,S)
%
% INPUT
%   A  Dataset
%   N  Number of points (optional; default: 50)
%   K  Number of nearest neighbors (optional; default: 1)
%   S  Standard deviation (optional; default: 1)
%
% OUTPUT
%   B  Generated dataset
%
% DESCRIPTION 
% Generation of N points using the K-nearest neighbors of objects in the 
% dataset A. First, N points of A are chosen in a random order. Next, to each 
% of these points and for each direction (feature), a Gaussian-distributed 
% offset is added with the zero mean and the standard deviation: S * the mean 
% signed difference between the point of A under consideration and its K 
% nearest neighbors in A. 
%
% The result of this procedure is that the generated points follow the local
% density properties of the point from which they originate.
%
% If A is a multi-class dataset the above procedure is followed class by
% class, neglecting objects of other classes and possibly unlabeled objects.
% 
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I. Default N is 100 objects per class.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, GENDAT, GENDATP

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function B = gendatk(varargin)

  argin = shiftargin(varargin,'vector');
	argin = setdefaults(argin,[],[],1,1);
  if mapping_task(argin,'definition')
    B = define_mapping(argin,'generator','KNN generation');
    return
  end
  
  % execution
  [A,N,k,stdev] = deal(argin{:});

  if isdataset(A)
    Adouble = false;
    A = prdataset(A);
    A = setlablist(A); % remove empty classes first
  else
    Adouble = true;
    A = prdataset(A,1);
  end

	[m,n,c] = getsize(A);
	prior = getprior(A);
	if isempty(N), 
		N = repmat(50,1,c); 				% 50 samples are generated.  		
	end
	N = genclass(N,prior);				% Generate class frequencies according to the priors.			

	lablist = getlablist(A);
	B = [];
	labels = [];
	% Loop over classes.
	for j=1:c
		a = getdata(A,j); 					% The j-th class.
		[D,I] = sort(distm(a)); 
		I = I(2:k+1,:); 						% Indices of the K nearest neighbors.
		alf = randn(k,N(j))*stdev;	% Normally distributed 'noise'.
		nu = ceil(N(j)/size(a,1));	% It is possible that NU > 1 if many objects have to be generated. 
		J = randperm(size(a,1));		
		J = repmat(J,nu,1)';				
		J = J(1:N(j));							% Combine the NU repetitions of J into one column vector.
		b = zeros(N(j),n);

		% Loop over features.
		for f = 1:n
%      Take all objects given by J, consider feature F.
%      Their K nearest neighbors are given by I(:,J)
%      We reshape them as a N(j) by K matrix (N(j) is the length of J)
%      Compute all differences between them and the original objects
%      Multiply these differences by the std dev stored in alf
%      Transpose and sum over the K neighbors, normalize by K
%      Transpose again and add to the original objects 
			 b(:,f) = a(J,f) + sum(( ( a(J,f)*ones(1,k) - ...
								reshape(+a(I(:,J),f),k,N(j))' ) .* alf' )' /k, 1)';
		end
		B = [B;b];
		labels = [labels; repmat(lablist(j,:),N(j),1)];
	end

  if Adouble
    B = +B;
  else
    B = prdataset(B,labels,'prior',A.prior);
    B = set(B,'featlab',getfeatlab(A),'name',getname(A),'featsize',getfeatsize(A));
  end

return;
