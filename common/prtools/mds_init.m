%MDS_INIT Initialization for MDS (variants of Sammon) mapping
% 
% 	Y = MDS_INIT (D,N,INIT)
% 
% INPUT
% 	D 		Square dissimilarity matrix of the size M x M
% 	N 		Desired output dimensionality (optional; default: 2)
% 	INIT	Initialization method (optional; default: 'randnp')
%
% OUTPUT
% 	Y 		Initial configuration for the MDS method
%
% DESCRIPTION  
% Finds a configuration of points Y in an N-dimensional space, used as
% a starting configuration for an MDS mapping based on the distance matrix D.
% The parameter INIT is a string standing for the initialization method:
%
%  'randp'  - linear mapping of D on n randomly (uniform distribution) chosen vectors
%  'randnp' - linear mapping of D on n randomly (normal distribution) chosen vectors
%  'randv'  - randomly (uniform distribution) chosen vectors from D
%  'maxv'   - n columns of D with the largest variances
%  'kl'     - Karhunen Loeve projection (linear mapping) of D (first n eigenvectors)
%  'cs'     - Classical Scaling
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, MDS_CS, MDS
%

% Undocumented use: if W is an already trained MDS, then for adding new points Dnew
%     W = mds_init(W,Ystart) 
% would assign Ystart as an initial configuration for he new MDS.
%

%
% Copyright: Elzbieta Pekalska, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function [Y,w] = mds_init(D,n,initm)

w = [];  
	if isa(D,'prmapping') 			% assign an initial configuration for the new MDS
		ww = getdata(D);
		if size(n,2) ~= size(ww{1},2),
			error('The dimensionality of the MDS map and the initial configuration does not match.')
		end
		if nargin < 3,
		 ;
		end	
    ww{6} = n;
    Y     = setdata(D,ww);
    return
    
	else

		if nargin < 2, n = 2; end
		if nargin < 3, initm = 'randnp'; end

    if isa(D,'prdataset'),
			lab = getlab(D);
			D=+D;
		end
    D = remove_nan(D);
		[m,mm] = size(D);
		w = [];
		switch initm		
			case 'randp',
				Y = D * rand(mm,n)/10;

			case 'randnp',
				Y = D * randn(mm,n)/10;

			case 'randv',
				Y = rand(m,n);

			case 'maxv',
				U     = std(D);
				[V,I] = sort(-U);
				Y     = D(:,I(1:n));    

			case 'kl',
				options.disp = 0; 
				[E,L] = eigs(prcov(+D),n,'lm',options);
				 Y    = D * E/100;   

			case 'cs',
				if length(n) == 1 & m == mm,
					w = mds_cs(D,n);   
					Y = D*w;
				else
					error('The CS initialization is not possible');
				end

			otherwise,   
				error('The possible initializations are: randp, randnp, randv, maxv, cs or kl.');  
		end
	end
return;



















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
	% check whether D has a zero diagonal
  if length(intersect(find(D(:) < eps), 1:m+1:(m*m))) >= m,
    D = (D+D')/2;
  end
end
