%GENDATP Parzen density data generation
% 
%   B = GENDATP(A,N,S,G)
%   B = A*GENDATP([],N,S,G)
%   B = A*GENDATP(N,S,G)
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, GENDAT, GENDATK

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function B = gendatp(varargin)

  argin = shiftargin(varargin,'vector');
	argin = setdefaults(argin,[],[],0,[]);
  if mapping_task(argin,'definition')
    B = define_mapping(argin,'generator','Parzen generation');
    return
  end
  
  % execution
  [A,N,s,G] = deal(argin{:});

  if isdataset(A)
    Adouble = false;
    A = prdataset(A);
    A = setlablist(A); % remove empty classes first
  else
    Adouble = true;
    A = prdataset(A,1);
  end

	[m,k,c] = getsize(A);
	p = getprior(A);
	if (isempty(N)) 
		N = repmat(50,1,c); 
  end

	if (length(s) == 1)
		s = repmat(s,1,c);
	end
	if (length(s) ~= c)
		error('Wrong number of smoothing parameters.')
	end

	if (isempty(G))
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
  
  if Adouble
    B = +B;
  else
    B = prdataset(B,labels,'prior',A.prior);
    B = set(B,'featlab',getfeatlab(A),'name',getname(A),'featsize',getfeatsize(A));	
  end

return
