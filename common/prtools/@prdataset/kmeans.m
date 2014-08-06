%KMEANS PRTools k-means clustering, deprecated
%
%   [LABELS,B] = KMEANS(A,K,MAXIT,INIT)
%
% INPUT
%  A       Dataset
%  K       Number of clusters to be found (optional; default: 2)
%  MAXIT   maximum number of iterations (optional; default: 50)
%  INIT    Labels for initialisation, or
%          'rand'     : take at random K objects as initial means, or
%          'kcentres' : use KCENTRES for initialisation (default)
%
% OUTPUT
%  LABELS  Cluster assignments, 1..K
%  B       Dataset with original data and labels LABELS: 
%          B = PRDATASET(A,LABELS)
% 
% DESCRIPTION
% K-means clustering of data vectors in A. This routine calls PRKMEANS
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRKMEANS, HCLUST, KCENTRES, MODESEEK, EMCLUST, PRPROGRESS
