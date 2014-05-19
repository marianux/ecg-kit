%GENDATP Parzen density data generation
% 
%   B = GENDATP(A,N,S,G)
% 
% INPUT
%   A  Dataset
%   N  Number(s) of points to be generated (optional; default: 50 per class)
%   S  Smoothing parameter(s) 
%      (optional; default: a maximum likelihood estimate based on A)
%   G  Covariance matrix used for generation of the data 
%      (optional; default: the identity matrix)
%
% OUTPUT
%   B  Dataset of points generated according to Parzen density
%
% DESCRIPTION  
% Generation of a dataset B of N points by using the Parzen estimate of the
% density of A based on a smoothing parameter S. N might be a row/column 
% vector with different numbers for each class. Similarly, S might be 
% a vector with different smoothing parameters for each class. If S = 0, 
% then S is determined by a maximum likelihood estimate using PARZENML. 
% If N is a vector, then exactly N(I) objects are generated for the class I. 
% G is the covariance matrix to be used for generating the data. G may be 
% a 3-dimensional matrix storing separate covariance matrices for each class.
% 
% SEE ALSO
% DATASETS, GENDATK

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gendatp.m,v 1.5 2010/03/25 15:41:39 duin Exp $

function B = gendatp(A,N,s,G)

		if (nargin < 1)
		error('No dataset found');
	end

	A = prdataset(A);
	A = setlablist(A); % remove empty classes first

	[m,k,c] = getsize(A);
	p = getprior(A);
	if (nargin < 2) 
		prwarning(4,'Number of points not specified, 50 per class assumed.');
		N = repmat(50,1,c); 
	end
	if (nargin < 3) 
		prwarning(4,'Smoothing parameter(s) not specified, to be determined be an ML estimate.');
		s = 0; 
	end

	if (length(s) == 1)
		s = repmat(s,1,c);
	end
	if (length(s) ~= c)
		error('Wrong number of smoothing parameters.')
	end

	if (nargin < 4)
		prwarning(4,'Covariance matrices not specified, identity matrix assumed.');
		covmat = 0; 			% covmat indicates whether a covariance matrix should be used
											% 0 takes the identity matrix as the covariance matrix
	else
		covmat = 1;
		if (ndims(G) == 2)
			G = repmat(G,[1 1 c]);
		end
		if any(size(G) ~= [k k c])
			error('Covariance matrix has a wrong size.')
		end
	end
	
	N = genclass(N,p);
	lablist = getlablist(A);

	B = [];
	labels = [];
	% Loop over classes.
	for j=1:c
		a = getdata(A,j);
		a = prdataset(a);
		ma = size(a,1);
		if (s(j) == 0)				% Estimate the smoothing parameter.
			h = parzenml(a);
		else
			h = s(j);
		end
		if (~covmat)
			b = a(ceil(rand(N(j),1) * ma),:) + randn(N(j),k).*repmat(h,N(j),k);
		else 
			b = a(ceil(rand(N(j),1) * ma),:) + ...
			    gendatgauss(N(j),zeros(1,k),G(:,:,j)).*repmat(h,N(j),k);
		end

		B = [B;b];
		labels = [labels; repmat(lablist(j,:),N(j),1)];
	end
	B = prdataset(B,labels);
	B = setprior(B,p);
	B = set(B,'featlab',getfeatlab(A),'name',getname(A),'featsize',getfeatsize(A));	

return
