%HCLUST hierarchical clustering
% 
%  [LABELS, DENDROGRAM] = HCLUST(D,TYPE,K)
%   DENDROGRAM = HCLUST(D,TYPE)
% 
% INPUT
%  D     dissimilarity matrix
%  TYPE  string name of clustering criterion (optional)
%        - 's' or 'single'   : single linkage (default)
%        - 'c' or 'complete' : complete linkage
%        - 'a' or 'average'  : average linkage (weighted over cluster sizes)
%  K     number of clusters (optional)
%
% OUTPUT
%  LABELS       vector with labels
%  DENDROGRAM   matrix with dendrogram
%
% DESCRIPTION 
% Computation of cluster labels and a dendrogram between the clusters for
% objects with a given distance matrix D. K is the desired number of
% clusters. The dendrogram is a 2*K matrix. The first row yields all
% cluster sizes. The second row is the cluster level on which the set of
% clusters starting at that position is merged with the set of clusters
% just above it in the dendrogram. A dendrogram may be plotted by PLOTDG.
%
%   DENDROGRAM = HCLUST(D,TYPE)
%
% As in this case no clustering level is supplied, just the entire
% dendrogram is returned. The first row now contains the object indices.
%
% EXAMPLE
% a = gendats([25 25],20,5);    % 50 points in 20 dimensional feature space
% d = sqrt(distm(a));           % Euclidean distances
% dendg = hclust(d,'complete'); % dendrogram
% plotdg(dendg)
% lab = hclust(d,'complete',2); % labels
% confmat(lab,getlabels(a));    % confusion matrix
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLOTDG, PRKMEANS, KCENTRES, MODESEEK, EMCLUST

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: hclust.m,v 1.3 2008/02/14 11:54:43 duin Exp $

function [labels, dendrogram] = hclust(D,type,k)   
		
	if nargin < 3, k = []; end
	if nargin < 2, type = 's'; end
	[m,m1] = size(D);
	if m ~= m1
		error('Input matrix should be square')
	end
	D = D + diag(inf*ones(1,m));     % set diagonal at infinity.
	W = [1:m+1];               % starting points of clusters in linear object set.
	V = [1:m+2];               % positions of objects in final linear object set.
	F = inf * ones(1,m+1);     % distance of next cluster to previous cluster to be stored at first point of second cluster
	Z = ones(1,m); % number of samples in a cluster (only for average linkage)

  t = sprintf('Analysing %i cluster levels: ',m);
  prwaitbar(m,t);
	for n = 1:m-1
    prwaitbar(m,n,[t int2str(n)]);
		% find minimum distance D(i,j) i<j
		[di,I] = min(D); [dj,j] = min(di); i = I(j);
		if i > j, j1 = j; j = i; i = j1; end
		% combine clusters i,j
		switch type
		case {'s','single'}
			D(i,:) = min(D(i,:),D(j,:));
		case {'c','complete'}
			D(i,:) = max(D(i,:),D(j,:));
		case {'a','average'}
			D(i,:) = (Z(i)*D(i,:) + Z(j)*D(j,:))/(Z(i)+Z(j));
			Z(i:j-1) = [Z(i)+Z(j),Z(i+1:j-1)]; Z(j) = [];
		otherwise
			error('Unknown clustertype desired')
		end
		D(:,i) = D(i,:)';
		D(i,i) = inf;
		D(j,:) = []; D(:,j) = [];
		% store cluster distance
		F(V(j)) = dj;
		% move second cluster in linear ordering right after first cluster
		IV = [1:V(i+1)-1,V(j):V(j+1)-1,V(i+1):V(j)-1,V(j+1):m+1];
		W = W(IV); F = F(IV);   
		% keep track of object positions and cluster distances
		V = [V(1:i),V(i+1:j) + V(j+1) - V(j),V(j+2:m-n+3)];
  end
  prwaitbar(0);
	if ~isempty(k) | nargout == 2
		if isempty(k), k = m; end
		labels = zeros(1,m); 
		[S,J] = sort(-F);      % find cluster level
		I = sort(J(1:k+1));    % find all indices where cluster starts
		for i = 1:k            % find for all objects cluster labels
			labels(W(I(i):I(i+1)-1)) = i * ones(1,I(i+1)-I(i));
		end                    % compute dendrogram
		dendrogram = [I(2:k+1) - I(1:k); F(I(1:k))];
		labels = labels';
	else
		labels = [W(1:m);F(1:m)]; % full dendrogram
	end
return
	
