% MDS_CS Trainable mapping for classical scaling
%
% 	W = MDS_CS(D,N)
% 	W = D*MDS_CS([],N)
% 	W = D*MDS_CS(N)
% 
% INPUT
% 	D 	Square dissimilarity matrix of the size M x M
% 	N 	Desired output dimensionality (optional; default: 2)
%
% OUTPUT
% 	W		Classical scaling mapping 
%
% DESCRIPTION  
% A linear mapping of objects given by a symmetric distance matrix D with
% a zero diagonal onto an N-dimensional Euclidean space such that the square 
% distances are preserved as much as possible. 
%
% D is assumed to approximate the Euclidean distances, i.e. 
% D_{ij} = sqrt(sum_k (x_i-x_j)^2). 
%
%
% REFERENCES
% 1. I. Borg and P. Groenen, Modern Multidimensional Scaling, Springer Verlag, Berlin, 
% 	 1997.
% 2. T.F. Cox and M.A.A. Cox, Multidimensional Scaling, Chapman and Hall, London, 1994.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, MDS_INIT, MDS_CS

%
% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function W = mds_cs (varargin)

  mapname = 'Classical MDS';
	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],2);
  
  if mapping_task(argin,'definition')
    
    W = define_mapping(argin,'fixed',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [D,n] = deal(argin{:});
    [m,mm] = size(D);
    % Definition of the projection;
    if ~issym(D,1e-12),
      error('Matrix should be symmetric.')
    end

    D = remove_nan(D);  
    D = +D.^2;
    B = -repmat(1/m,m,m);					
    B(1:m+1:end) = B(1:m+1:end) + 1;		% B = eye(m) - ones(m,m)/m
    B = -0.5 * B * D * B; 							% B is now the matrix of inner products 

    ok = 0;
    p  = 1;
    k  = min (p*n,m);
    options.disp   = 0; 
    options.isreal = 1; 
    options.issym  = 1; 

    while ~ok & k <= m, 								% Find k largest eigenvalues
      k = min (p*n,m);
      [V,S,U] = eigs(B,k,'lm',options);
      s = diag(real(S));
      [ss,I] = sort(-s);
      ok = sum(s>0) >= n;
      p  = p+1;
    end
    if ~ok & k == m,
      error ('The wanted configuration cannot be found.');
    end
    J = I(1:n);
    W = prmapping(mfilename,'trained',{V(:,J), mean(D,2)', s(J)},[],m,length(J));
    W = setname(W,mapname);
  
  else % Evaluation
    
    [D,w] = deal(argin{:});
    [m,n] = size(D);
	  if isa(D,'prdataset'),
      lab = getlab(D); 
			D = +D;
		else
			lab = [];
		end
    D = remove_nan(D);  
    D     = D.^2;
    Q     = w{1}; 
    me    = w{2};
    L     = w{3};

    Z = -repmat(1/n,n,n);                 
    Z(1:n+1:end) = Z(1:n+1:end) + 1;    	 % Z = eye(n) - ones(n,n)/n
    data = -0.5 * (D - me(ones(m,1),:)) * Z * Q * diag(sqrt(abs(L))./L);          
    if ~isempty(lab), 
      W = prdataset(data,lab); 
    else
      W = data;
    end
    
  end
  
return



















































function D = remove_nan(D)
% Remove all appearing NaN's by replacing them by the nearest neighbors.
%
[m,mm] = size(D);

nanindex = find(isnan(D(:)));  
if ~isempty(nanindex),
  for i=1:m
    K = find(isnan(D(i,:)));
    I = 1:mm; 
    [md,kk] = min(D(i,:)); 
    if md < eps, 
      I(kk) = [];
    end
    D(i,K) = min(D(i,I));
  end
  % check whether the diagonal is of all zeros
  if length(intersect(find(D(:) < eps), 1:m+1:(m*m))) >= m,
    D = (D+D')/2;
  end
end
