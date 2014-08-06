%PCAKLM Principal Component Analysis/Karhunen-Loeve Mapping
%       (PCA or MCA of overall/mean covariance matrix)
% 
% 	[W,FRAC] = PCAKLM(TYPE,A,N)
% 	[W,N]    = PCAKLM(TYPE,A,FRAC)
%
% INPUT
%  A           Dataset
%  TYPE        Type of mapping: 'pca' or 'klm'. Default: 'pca'.
%  N or FRAC   Number of dimensions (>= 1) or fraction of variance (< 1) 
%              to retain; if > 0, perform PCA; otherwise MCA.
%              Default: N = inf.
%
% OUTPUT
% W            Affine Karhunen-Loeve mapping
% FRAC or N    Fraction of variance or number of dimensions retained.
%
% DESCRIPTION
% Performs a principal component analysis (PCA) or minor component analysis
% (MCA) on the overall or mean class covariance matrix (weighted by the
% class prior probabilities). It finds a rotation of the dataset A to an
% N-dimensional linear subspace such that at least (for PCA) or at most (for
% MCA) a fraction FRAC of the total variance is preserved.
%
% PCA is applied when N (or FRAC) >= 0; MCA when N (or FRAC) < 0. If N is 
% given (abs(N) >= 1), FRAC is optimised. If FRAC is given (abs(FRAC) < 1), 
% N is optimised. 
%
% Objects in a new dataset B can be mapped by B*W, W*B or by A*KLM([],N)*B.
% Default (N = inf): the features are decorrelated and ordered, but no 
% feature reduction is performed.
%
% ALTERNATIVE
%
% 	V = PCAKLM(A,TYPE,0)
% 
% Returns the cumulative fraction of the explained variance. V(N) is the 
% cumulative fraction of the explained variance by using N eigenvectors.
%
% This function should not be called directly, only trough PCA or KLM.
% Use FISHERM for optimizing the linear class separability (LDA).
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PCLDC, KLLDC, PCAM, KLM, FISHERM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: pcaklm.m,v 1.15 2010/02/08 15:29:48 duin Exp $

function [w,truefrac] = pcaklm (type,a,frac)

		truefrac = [];

	% Default: preserve all dimensions (identity mapping).
	if (nargin < 3) | (isempty(frac))
		frac = inf; 
		prwarning (3,'no dimensionality given, only decorrelating and ordering dimensions');
	end

	% Default: perform PCAM.
	if (nargin < 1) | (isempty(type))
		type = 'pcam';
		prwarning (3,'no type given, assuming PCA');
	end

	if (strcmp(type,'pcam'))
		mapname = 'PCA';
	elseif (strcmp(type,'klm'))
		mapname = 'Karhunen-Loeve Mapping';
	else
		error('Unknown type specified');
	end

	%DXD Make the name a bit more informative:
	if isfinite(frac)
		if (frac<1)
			mapname = [mapname sprintf(' ret. %4.1f%% var',100*frac)];
		else
			mapname = [mapname sprintf(' to %dD',frac)];
		end
	end

	% Empty mapping: return straightaway.
	if (nargin < 2) | (isempty(a))
		w = prmapping(type,frac);
		w = setname(w,mapname);
		return
	end
	
	%nodatafile(a);
	if ~isdataset(a)
		if isa(a,'double')
			a = prdataset(a,1);   % make sure we have a dataset
		else
			error('nodatafile','Data should be given in a dataset or as doubles')
		end
	end
	
	islabtype(a,'crisp','soft');
	isvaldfile(a,1);   % at least 1 object per class
	a = setfeatdom(a,[]);  % get rid of domain testing

	[m,k,c] = getsize(a);
	p = getprior(a);
	a = setprior(a,p);  % make class frequencies our prior

	% If FRAC < 0, perform minor component analysis (MCA) instead of 
	% principal component analysis.
	mca = (frac < 0); frac = abs(frac);

	% Shift mean of data to origin.
	b = a*scalem(a); 
	
	% If there are less samples M than features K, first perform a lossless
	% projection to the (M-1) dimensional space spanned by the samples.
	if (m <= k)
		testdatasize(b,'objects');
		u = reducm(b); b = b*u;
		korg = k; [m,k] = size(b);
		frac = min(frac,k);
	else
		testdatasize(b,'features');
		u = [];
	end

	% Calculate overall or average class prior-weighted covariance matrix and
	% find eigenvectors F. 

	if (strcmp(type,'pcam'))
		if (c==0) | ~islabtype(a,'crisp')  % we have unlabeled data!
			G = prcov(+b); % use all
		else
			bb = [];
			classsiz = classsizes(b);
			for j = 1:c
				%bb = [bb; seldat(b,j)*filtm([],'double')*p(j)/classsiz(j)];
				bb = [bb; double(+seldat(b,j))*p(j)/classsiz(j)];
			end
			%[U,G] = meancov(remclass(setnlab(bb*m,1)));
			[U,G] = meancov(bb);
		end
	else
		%DXD For high dimensional dataset with many classes, we cannot
		%store all individual cov. matrices in memory (like in the next
		%line), but we have to compute them one by one:
		%[U,GG] = meancov(b,1);
		G = zeros(k,k);
		for i = 1:c
			%G = G + p(i)*GG(:,:,i);
			[U,GG] = meancov(seldat(b,i),1);
			G = G + p(i)*GG;

		end
	end
	[F,V] = preig(G); % overdone if reducm has been called

	% v = V(I) contains the sorted eigenvalues:
	% descending for PCAM, ascending for MCA.
	if (mca)
		[v,I] = sort(diag(V));
	else
		[v,I] = sort(-diag(V)); v = -v;
	end
	
	if (frac == inf)           % Return all dimensions, decorrelated and ordered.
		n = k; truefrac = k;
	elseif (frac == 0)         % Just return cumulative retained variance.
		w = cumsum(v)/sum(v);
		return
	elseif (frac >= 1)         % Return FRAC dimensions.
		n = abs(frac); 
		if (n > k),
			error('illegal dimensionality requested');
		end
		I = I(1:n);
		sv = sum(v); 
		if (sv ~= 0),
			truefrac = cumsum(v(1:n))/sv;
		else,
			truefrac = 0;
		end;
	elseif (frac > 0)         % Return the N dimensions that retain at least (PCA)
                              % or at most (MCA) FRAC variance.
		J = find(cumsum(v)/sum(v) > frac);
		if (mca), n = J(1)-1; else, n = J(1); end;
		truefrac = n; I = I(1:n);
	end

	% If needed, apply pre-calculated projection to (M-1) dimensional subspace.
	if (~isempty(u))
		rot = u.data.rot*F(:,I); 
		off = u.data.offset*F(:,I);
	else
		rot = F(:,I); 
		off = -mean(a*F(:,I));
	end

	% Construct affine mapping.
	
	w = affine(rot,off,a);
  w = setdata(w,v,'eigenvalues');
	w = setname(w,mapname);
		
return
