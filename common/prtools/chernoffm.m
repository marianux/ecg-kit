%CHERNOFFM Suboptimal discrimination linear mapping (Chernoff mapping)
%
%	W = CHERNOFFM(A,N,R)
% 
% INPUT
%	A   Dataset
%	N   Number of dimensions to map to, N < C, where C is the number of classes
%	    (default: min(C,K)-1, where K is the number of features in A)
%	R   Regularization variable, 0 <= r <= 1, default is r = 0, for r = 1 the 
%	    Chernoff mapping is (should be) equal to the Fisher mapping
%
% OUTPUT
% 	W   Chernoff mapping
%
% DESCRIPTION  
% Finds a mapping of the labeled dataset A onto an N-dimensional linear
% subspace such that it maximizes the heteroscedastic Chernoff criterion
% (also called the Chernoff mapping).
%
% REFERENCE
% M. Loog and R.P.W. Duin, Linear Dimensionality Reduction via a Heteroscedastic
% Extension of LDA: The Chernoff Criterion, IEEE Transactions on pattern
% analysis and machine intelligence, vol. PAMI-26, no. 6, 2004, 732-739.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FISHERM, NLFISHERM, KLM, PCA

% Copyright: M. Loog, marco@isi.uu.nl
% Image Sciences Institute, University Medical Center Utrecht 
% P.O. Box 85500, 3508 GA Utrecht, The Netherlands

% $Id: chernoffm.m,v 1.3 2010/02/08 15:31:48 duin Exp $

function W = chernoffm(a,n,r)

		if (nargin < 3)
        r = 0;
	end

	if (nargin < 2)
		prwarning (4, 'number of dimensions to map to not specified, assuming min(C,K)-1');
		n = []; 
	end

	% If no arguments are specified, return an untrained mapping.

	if (nargin < 1) | (isempty(a))
		W = prmapping('chernoffm',{n,r});
		W = setname(W,'Chernoff mapping');
		return
	end

	isvaldset(a,2,2); % at least 2 objects per class, 2 classes

	% If N is not given, set it. 

	[m,k,c] = getsize(a); if (isempty(n)), n = min(k,c)-1; end
	if (n >= m) 
		error('The dataset is too small for the requested dimensionality')
	end

	% To simplify computations, the within scatter W is set to the identity matrix 
    % therefore A is centered and sphered by KLMS.

	v = klms(a); a = a*v;	
	k = size(v,2); % dimensionalities may have changed for small training sets
	
	% Now determine a solution to the Chernoff criterion

	[U,G] = meancov(a); 

	CHERNOFF = zeros(k);
    p = getprior(a);
    
    for i = 1:c-1 % loop over all class pairs
        for j = i+1:c
            p_i = p(i)/(p(i)+p(j)); p_j = p(j)/(p(i)+p(j));
            G_i = (1-r)*G(:,:,i) + r*eye(k); G_j = (1-r)*G(:,:,j) + r*eye(k);
            G_ij = p_i*G_i + p_j*G_j;
            m_ij = sqrtm(real(prinv(G_ij)))*(U(i,:) - U(j,:))';
            CHERNOFF = CHERNOFF + p(i) * p(j) * ...
                ( m_ij * m_ij' ...
                + 1/(p_i*p_j) * ( logm(G_ij) ...
                - p_i * logm(G_i) - p_j * logm(G_j) ) );
        end
    end

	[F,V] = preig(CHERNOFF); 			
	[dummy,I] = sort(-diag(V)); 
	I = I(1:n);
	rot = F(:,I);
	off = -mean(a*F(:,I));
	preW = affine(rot,off,a);
	W = v*preW;					% Construct final mapping.
	W = setname(W,'Chernoff mapping');

return
