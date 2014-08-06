%KCENTRES Finds K center objects from a distance matrix
% 
%  [LAB,J,DM] = KCENTRES(D,K,N)
% 
% INPUT
%   D    Distance matrix between, e.g. M objects (may be a dataset)
%   K    Number of center objects to be found (optional; default: 1)
%   N    Number of trials starting from a random initialization
%        (optional; default: 1)
%
% OUTPUT
%   LAB  Integer labels: each object is assigned to its nearest center
%   J    Indices of the center objects
%   DM   A list of distances corresponding to J: for each center in J 
%        the maximum distance of the objects assigned to this center.	
%
% DESCRIPTION  
% Finds K center objects from a symmetric distance matrix D. The center 
% objects are chosen from all M objects such that the maximum of the 
% distances over all objects to the nearest center is minimized. For K > 1,
% the results depend on a random initialization. The procedure is repeated 
% N times and the best result is returned. 
%
% If N = 0, initialisation is not random, but done by a systematic
% selection based on a greedy approach.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% HCLUST, PRKMEANS, EMCLUST, MODESEEK

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: kcentres.m,v 1.5 2008/07/03 09:08:43 duin Exp $

function [labels,Jopt,dm] = kcentres(d,k,n)
	
	if (nargin < 3) | isempty(n)
      n = 1;
      prwarning(4,'Number of trials not supplied, assuming one.');
   end

	if (nargin < 2) | isempty(k)
      k = 1;
      prwarning(4,'Number of centers not supplied, assuming one.');
   end
	
	if(isdataset(d))
		d = +d;
		prwarning(4,'Distance matrix is convert to double.');
	end

	[m,m2] = size(d);
	if ( ~issym(d,1e-12) | any(diag(d) > 1e-14) )
		error('Distance matrix should be symmetric and have zero diagonal')
	end

% checking for a zero diagonal 
	t = eye(m) == 1;	
	if(~all(d(t)==0))
		error('Distance matrix should have a zero diagonal.')
	end

	if (k == 1)
		dmax = max(d);
		[dm,Jopt] = min(dmax);
		labels = repmat(1,m,1);
		return;
	end

	if k > m
		error('Number of centres should not exceed number of objects')
	end
	
	% We are here only if K (> 1) centers are to be found.
	% Loop over number of trials.
	dmax = max(max(d));
	dopt = inf;
	s = sprintf('k-centres, %i attempts: ',n);
	prwaitbar(n,s,n>1);
	if n == 0
		nrep = 1; 
	else 
		nrep = n; 
	end
	for tri = 1:nrep
		prwaitbar(n,tri,[s int2str(tri)]);
		if n == 0
			M = kcentresort(d,k);           % systematic initialisation
		else
			M = randperm(m); M = M(1:k);		% Random initializations
		end
 		J = zeros(1,k); 					    	  % Center objects to be found.

		% Iterate until J == M. See below.
		while 1,
			[dm,I] = min(d(M,:));

			% Find K centers. 
			for i = 1:k                   
				JJ = find(I==i); 
				if (isempty(JJ))
					%JJ can be empty if two or more objects are in the same position of 
					% feature space in dataset
					J(i) = 0;									
				else
						% Find objects closest to the object M(i)
					[dummy,j,dm] = kcentres(d(JJ,JJ),1,1);
					J(i) = JJ(j);
				end
			end

			J(find(J==0)) = [];
			k = length(J);
			if k == 0
				error('kcentres fails as some objects are identical: add some noise')
			end
			if (length(M) == k) & (all(M == J))
				% K centers are found.
				[dmin,labs] = min(d(J,:));
				dmin = max(dmin);
				break;
			end
			M = J;
		end

		% Store the optimal centers in JOPT.
	
		if (dmin <= dopt)
			dopt = dmin;
			labels = labs';
			Jopt = J;
		end
	
	end
	prwaitbar(0)
	
	% Determine the best centers over the N trials.
	dm = zeros(1,k);   
	for i=1:k
		L = find(labels==i);
		dm(i) = max(d(Jopt(i),L));
	end

return;

%KCENTRESORT Sort objects given by dissimilarity matrix
%
%               N = KCENTRESORT(D,P,CRIT)
%
% INPUT
%   D     Square dissimilarity matrix, zeros on diagonal
%   P     Number of prototypes to be selected
%   CRIT  'dist' or 'centre'
%
% OUTPUT
%   N     Indices of selected prototypes
%
% DESCRIPTION
% Sort objects given by square dissim matrix D using a greedy approach
% such that the maximum NN distance from all objects (prototypes)
% to the first K: max(min(D(:,N(1:K),[],2)) is minimized.
%
% This routines tries to sample the objects such that they are evenly
% spaced judged from their dissimilarities. This may be used as
% initialisation in KCENTRES. It works reasonably, but not very good.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% KCENTRES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function N = kcentresort(d,p,crit);

d = +d;
m = size(d,1);
if nargin < 3, crit = 'dist'; end
if nargin < 2 | isempty(p), p = m; end
L = [1:m];
N = zeros(1,p);
[dd,n] = min(max(d,[],2));  % this is the first (central) prototype
e = d(:,n);                 % store here the distances to the nearest prototype (dNNP)
f = min(d,repmat(e,1,m));   % replace distances that are larger than dNNP by dNNP
N(1) = n;                   % ranking of selected prototypes
L(n) = [];                  % candidate prototypes (all not yet selected objects)

for k=2:p                   % extend prototype set
        if strcmp(crit,'centre')
        [dd,n] = min(max(f(L,L),[],1)); % This selects the next prototype out of candidates in L
        e = min([d(:,L(n)) e],[],2);    % update dNNP
        f = min(d,repmat(e,1,m));       % update replacement of distances that are larger
                                              % than dNNP by dNNP
        elseif strcmp(crit,'dist')
                [dd,n] = max(mean(d(N(find(N > 0)),L)));
        else
                error('Illegal crit')
        end
  N(k) = L(n);                    % update list of selected prototypes
  L(n) = [];                      % update list of candidate prototypes
end

