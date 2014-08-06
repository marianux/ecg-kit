%MODECLUST KNN mode-seeking clustering
% 
% 	LAB = MODECLUST(A,K)
% 
% INPUT
%   A       Dataset
%   K       Number of neighbours to search for local mode
%
% OUTPUT
  %   LAB   Indices of modal samples
%
% DESCRIPTION
% A K-NN modeseeking method is used to assign objects to their nearest mode.
% Object densities are defined by one over the distance to the K-th nearest 
% neighbour. Clusters are defined by recursively jumping for every object 
% to object with the highest density in the local neighborhood.
%
% K can also be a vector of neighborhood sizes, which is much faster.
% Default K: a set of values determined by a geometric series.
%
% EXAMPLE
% delfigs
% a = gendatm(1000);      % generate 1000 objects in 8 
% lab = modeclust(a); 
% for j=1:size(lab,2)
%   nclust = numel(unique(lab(:,j)));
%   if nclust < 20 & nclust > 1
%     figure; scatterd(prdataset(a,lab(:,j)));
%     title(['Number of clusters: ' int2str(nclust)]);
%   end
% end
% showfigs
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PRKMEANS, HCLUST, KCENTRES, PROXM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function labels = modeclust(a,K)

	% prepare data
  a = +a;                 % get rid of possible PRTools dataset structure
	m = size(a,1);          % number of objects
	if (nargin < 2  | isempty(K))
    K = intgeos(2,1.2,m); % default set of nb-sizes (neighborhood-sizes)
  end
  nk = numel(K);          % number of neighborhood sizes (clusterings)
	f = zeros(m,nk);        % space for storing densities
  batchsize = 100;        % batch sizes (avoid large distance matrices)
  nbatch = ceil(m/batchsize); % number of batches
  
  % first loop over the data: compute densities for every nb-size
	t = sprintf('Computing %i densities: ',m);
	prwaitbar(m,t);
  for j=1:nbatch
		prwaitbar(m,j*batchsize,[t num2str(j*batchsize)]);
    batch = (j-1)*batchsize+1:j*batchsize;  % take a subset of the objects
    if any(batch>m)              % arrange last batch
      batch = batch(batch <= m);
    end
		prwaitbar(m,(j-1)*batchsize,[t num2str((j-1)*batchsize)]);
    dj = distm(a(batch,:),a);    % distances between batch and all data
    [dd,L] = sort(dj,2);         % sort them per object
    f(batch,:) = 1./dd(:,K);     % density estimates for all nb-sizes
  end
  prwaitbar(0);
  
  % second loop: find indices of maximum density point in every nb
  NN = zeros(m,nk);              % space for storing these indices
	t = sprintf('Checking %i neighborhoods: ',m);
	prwaitbar(m,t);
  for j=1:nbatch
		prwaitbar(m,j*batchsize,[t num2str(j*batchsize)]);
    batch = (j-1)*batchsize+1:j*batchsize;   % take a subset of the objects
    if any(batch>m)              % arrange last batch
      batch = batch(batch <= m)';
    end
    mj = numel(batch);           % number of data points in batch
    dj = distm(a(batch,:),a);    % recomputes distances
    [dd,J] = sort(dj,2);         % sort them per object
    for n=1:nk                   % run over all nb's
      ff = f(:,n);               % densities for this nb-size
      [dummy,I] = max(reshape(ff(J(:,1:K(n))),mj,K(n)),[],2);
      NN(batch,n)  = J([1:mj]'+(I-1)*mj); % for all points the index of
    end                          %  point with maximum local density 
  end
  prwaitbar(0);
  
  % third loop: follow gradient to density modes
	t = sprintf('Computing %i clusterings: ',nk);
	prwaitbar(nk,t);
  for n=1:nk
		prwaitbar(nk,n,[t num2str(n)]);
		N = NN(:,n);                 % initial nb points with maximum density
		% Re-assign samples to the sample their nearest neighbour is assigned to.
		% Iterate until assignments don't change anymore. Samples that then point 
		% to themselves are modes; all other samples point to the closest mode.
		M = N(N);
    while (any(M~=N))
			N = M; M = N(N);
    end
    labels(:,n) = M;  % labels are the indices of the modes to which data 
		                  % points are assigned to. 
    % Data points assigned to the same mode constitute a cluster
	end
	prwaitbar(0);

return
