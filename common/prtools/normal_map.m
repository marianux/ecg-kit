%NORMAL_MAP Map a dataset on normal-density classifiers or mappings
%           Back-end routine
% 
%   F = NORMAL_MAP(A,W)
%
% INPUT
%   A   Dataset
%   W   Mapping
%
% OUTPUT
% 	F   Density estimation for classes in A	
%
% DESCRIPTION
% Maps the dataset A by the normal density based classifier or mapping W.
% For each object in A, F returns the densities for each of the classes or
% distributions stored in W. For classifiers, the densities are weighted
% by the class prior probabilities. This routine is automatically called for
% computing A*W if W is a normal density based classifier or a mapping.
%
% Use W = LOGDENS(W) (or W = W*LOGDENS) if absolute densities are not
% needed and a more accurate posterior probability is desired.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, QDC, UDC, LDC, GAUSSM, LOGDENS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: normal_map.m,v 1.16 2010/02/24 20:28:13 duin Exp $

function F = normal_map(varargin)
		
	if isstruct(varargin{1})  % this is a constructor call like W = normal_map(w,nlab,k,c)
		[w,nlab,k,c] = deal(varargin{:});
		deg = ndims(w.cov)-1; % linear or quadratic
		n = size(w.mean,1);   % number of components
		w.det = zeros(1,n);
		if (deg == 1)
			H = w.cov;
			if min(size(H)) == 1 % cov matrix stored as scalar or diagonal
				E = 1./H;
				w.det = repmat(sum(log(real(H+1e-16))),1,n) ;
			else
				if (prrank(H) < size(H,1))
      		prwarning(2,'Singular case, pseudo-inverse of the covariance matrix is used.');
					E = real(prpinv(H));
				else
					E = real(prinv(H));
				end
				w.det = repmat(sum(log(real(preig(H)+1e-16)+realmin)),1,n);
			end
			w.cov = E;
		elseif deg == 2
			w.det = zeros(1,n);
			for i=1:n
				H = w.cov(:,:,i);
				if (prrank(H) < size(H,1))
					prwarning(1,'Singular case, pseudo-inverse of the covariance matrix is used.');
					E = real(prpinv(H));
				else
					E = real(prinv(H));
				end
				w.cov(:,:,i) = E;
				w.det(i) = sum(log(abs(preig(H)+1e-16)+realmin));
			end
		else
			error('Illegal value for degree')
		end
		F = prmapping(mfilename,'trained',w,nlab,k,c);
		return
	end
	
	% Now we have an execution call like F = normal_map(A,W)
	
	[A,W] = deal(varargin{:});

	if isdatafile(A)
		F = prdataset(A*W); 
		return
	end

	w = +W;          % data field of W (fields: w.mean, w.cov, w.prior, w.nlab)
	                 % each of these data fields has a mean, cov, prior for separate
									 % Gaussian components. The nlab field assigns each component
									 % to a class.

	if ~isfield(w,'det') % if det-field is not available, we are dealing with an old def
		F = normal_map_old(A,W);
		return
	end

	[k,c] = size(W); % c is number of classes
	% DEG = 1 indicates a common cov. matrix and DEG = 2 - separate cov. matrices.
	deg = ndims(w.cov)-1; 		
	U = w.mean; G = w.cov; p = w.prior;
	if (abs(1-sum(p)) > 1e-6)
		error('Class or component probabilities do not sum to one.')
	end
	lablist = getlab(W);

	[m,ka] = size(A);
	if (ka ~= k), 
		error('Feature sizes of the dataset and the mapping do not match.'); 
	end

	n = size(U,1);			% Number of components.
	F = zeros(m,n); 		% Gaussian densities for each component to be computed.

	if (deg == 1)		
		E = G;
	end
	% Loop over components.
  t = sprintf('Processing %i components: ',n);
  prwaitbar(n,t);
	for i=1:n
    prwaitbar(n,i,[t num2str(i)]);
		% Shift A such that the mean lies at the origin.
		X = +A - ones(m,1)*U(i,:);
		if (deg == 2)		
			E = G(:,:,i);
		end
		if min(size(E)) == 1 % diagonal or scalar of cov matrix
			if max(size(E)) == 1, E = repmat(E,1,k); end
			F(:,i) = -0.5*sum(X.*X.*repmat(E(:)',m,1),2) - (w.det(i) + k*log(2*pi))*0.5;
		else
		% Gaussian distribution for the i-th component. Take log of density to preserve tails
			F(:,i) = -0.5*sum(X'.*(E*X'),1)' - (w.det(i) + k*log(2*pi))*0.5;
		end
			if (getout_conv(W) ~= 2) % take log of density to preserve tails
				F(:,i) = exp(F(:,i));
			end
  end
  prwaitbar(0)

	if isfield(w,'nlab') 
		% For mixtures of Gaussians. Relates components to classes
		nlab = w.nlab;
		cc = max(w.nlab);
	else
		nlab = [1:c];
		cc = c;
	end

	if (getout_conv(W) == 2)
		Cmax = max(F(:));  % scale to gain accuracy in the tails
		F = F - Cmax;
	end
	FF = zeros(m,cc);
	% Loop over true classes. Weight the probabilities by the priors.
	for j = 1:cc
		J = find(nlab == j); 
		if (getout_conv(W) == 2) % take log of density to preserve tails
			% difficult to get this right for MOGs, and moreover, probably not of
			% much help. Anyway, we give it a try.
			FF(:,j) = exp(F(:,J))*w.prior(J)';
			% like in parzen_map, use in tails just largest component
			L = find(FF(:,j) <= 1e-300); 
			N = find(FF(:,j) >  1e-300); 
			if ~isempty(L)
				[FM,R] = max(F(L,J),[],2);
				FF(L,j) = FM + log(w.prior(R)');
			end
			if ~isempty(N)
				FF(N,j) = log(FF(N,j));
			end
		else
			FF(:,j) = F(:,J) * w.prior(J)';
		end
	end
	
	if (getout_conv(W) == 2) % scale to gain accuracy in the tails
		Fmax = max(FF,[],2);   % really needed!!!!
		FF = FF - repmat(Fmax,1,cc);
		FF = exp(FF);
	else
		FF = FF + realmin; % avoid devision by 0 in computing posterios later
	end
	
	if isdataset(A)
		F = setdata(A,FF,lablist);
	else
		F = FF;
	end

return;

%NORMAL_MAP_OLD Map a dataset on normal-density classifiers or mappings
%               using an old classifier definition: 
%               cov inverse during execution
% 

function F = normal_map_old(A,W)

		w = +W;          % data field of W (fields: w.mean, w.cov, w.prior, w.nlab)
	                 % each of these data fields has a mean, cov, prior for separate
									 % Gaussian components. The nlab field assigns each component
									 % to a class. 
	[k,c] = size(W); % c is number of classes
	% DEG = 1 indicates a common cov. matrix and DEG = 2 - separate cov. matrices.
	deg = ndims(w.cov)-1; 		
	U = w.mean; G = w.cov; p = w.prior;
	if (abs(1-sum(p)) > 1e-6)
		error('Class or component probabilities do not sum to one.')
	end
	lablist = getlab(W);

	[m,ka] = size(A);
	if (ka ~= k), 
		error('Feature sizes of the dataset and the mapping do not match.'); 
	end

	n = size(U,1);			% Number of components.
	F = zeros(m,n); 		% Gaussian densities for each component to be computed.

	if (deg == 1)		
		H = G;
		if (prrank(H) < size(H,1))
      prwarning(2,'Singular case, pseudo-inverse of the covariance matrix is used.');
			E = real(prpinv(H));
		else
			E = real(prinv(H));
		end
	end
	% Loop over components.
	for i=1:n
		% Shift A such that the mean lies at the origin.
		X = +A - ones(m,1)*U(i,:);
		if (deg == 2)		
			H = G(:,:,i);
			if (prrank(H) < size(H,1))
				prwarning(1,'Singular case, pseudo-inverse of the covariance matrix is used.');
				E = real(prpinv(H));
			else
				E = real(prinv(H));
			end
		end
		% Gaussian distribution for the i-th component. Take log of density to preserve tails
		F(:,i) = -0.5*sum(X'.*(E*X'),1)' - (sum(log(real(preig(H)+1e-16)+realmin)) + k*log(2*pi))*0.5;
		if (getout_conv(W) ~= 2) % take log of density to preserve tails
			F(:,i) = exp(F(:,i));
		end
	end
	

	if isfield(w,'nlab') 
		% For mixtures of Gaussians. Relates components to classes
		nlab = w.nlab;
		cc = max(w.nlab);
	else
		nlab = [1:c];
		cc = c;
	end

	if (getout_conv(W) == 2)
		Cmax = max(F(:));  % scale to gain accuracy in the tails
		F = F - Cmax;
	end
	FF = zeros(m,cc);
	% Loop over true classes. Weight the probabilities by the priors.
	for j = 1:cc
		J = find(nlab == j); 
		if (getout_conv(W) == 2) % take log of density to preserve tails
			% difficult to get this right for MOGs, and moreover, probably not of
			% much help. Anyway, we give it a try.
			FF(:,j) = exp(F(:,J))*w.prior(J)';
			% like in parzen_map, use in tails just largest component
			L = find(FF(:,j) <= 1e-300); 
			N = find(FF(:,j) >  1e-300); 
			if ~isempty(L)
				[FM,R] = max(F(L,J),[],2);
				FF(L,j) = FM + log(w.prior(R)');
			end
			if ~isempty(N)
				FF(N,j) = log(FF(N,j));
			end
		else
			FF(:,j) = F(:,J) * w.prior(J)';
		end
	end
	
	if (getout_conv(W) == 2)
		FF = exp(FF);
	else
		FF = FF + realmin; % avoid devision by 0 in computing posterios later
	end

	if isdataset(A)
		F = setdata(A,FF,lablist);
	else
		F = FF;
	end

return;
