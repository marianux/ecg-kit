%PRKMEANS PRTools k-means clustering
%
%   [LABELS,B] = PRKMEANS(A,K,MAXIT,INIT)
%
% INPUT
%  A       Matrix or dataset
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
% K-means clustering of data vectors in A. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, HCLUST, KCENTRES, MODESEEK, EMCLUST, PRPROGRESS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [assign,a] = prkmeans(a,varargin)

  [kmax,maxit,init] = setdefaults(varargin,2,50,'kcentres');
	n_ini = 100;				% Maximum size of subset to use for initialisation.
  
	% Create dataset with all equal labels and no priors.
	m = size(a,1); 
	a = prdataset(a);
	islabtype(a,'crisp');
  a=set(a,'labels',ones(m,1),'lablist',[1:kmax]','prior',[]); % for speed
	
	n_ini = max(n_ini,kmax*5);  % initialisation needs sufficient samples 
  text = sprintf('k-means clustering, %i iterations: ',maxit);
	prwaitbar(maxit,text);
	% Initialise by performing KCENTRES on...
	if (size(init,1) == 1) & strcmp(init,'kcentres') & (m > n_ini) 
		%prwarning(2,'Initializing by performing KCENTRES on subset of %d samples.', n_ini);
		b = +gendat(a,n_ini); % ... a random subset of A.
		d = +distm(b);
		assign = kcentres(d,kmax,[]);
		bb = setprior(prdataset(b,assign),0);
		w = nmc(bb); % Initial partition W and assignments ASSIGN.
	elseif (size(init,1) == 1) & strcmp(init,'kcentres')
		%prwarning(2,'Initializing by performing KCENTRES on training set.');
		d = +distm(a);    % ... the entire set A.
		assign = kcentres(d,kmax,[]);
		aa = setprior(prdataset(a,assign),0);
		aa = setlablist(aa);
		w = nmc(aa); % mapping trained on the complete dataset
	elseif (size(init,1) == 1) & strcmp(init,'rand')
		%prwarning(2,'Initializing by randomly selected objects');
		R = randperm(m);
		w = nmc(prdataset(a(R(1:kmax),:),[1:kmax]')); % mapping trained on kmax random samples
	elseif (size(init,1) == m)
		assign = renumlab(init);
		kmax = max(assign);
		%prwarning(2,'Initializing by given labels, k = %i',kmax);
		w = nmc(prdataset(a,assign));
	else
		error('Wrong initialisation supplied')
	end
	assign = labeld(a*w);
	a = prdataset(a,assign);
	a = setprior(a,0);
	%tmp_assign = zeros(m,1); % Allocate temporary array.

	% Main loop, while assignments change
	it=1; % number of iterations
	ndif = 1;

	while (it<maxit) & (ndif > 0)
		prwaitbar(maxit,it,[text num2str(it)]);
		tmp_assign = assign;     % Remember previous assignments.
		a = setnlab(a,assign);	
    a = setlablist(a);       % remove empty classes
    % disp([it classsizes(a)])
		w = a * nmc;   % Re-partition the space by assigning samples to nearest mean.
		[dummy,assign] = max(+a*w,[],2);  % Re-calculate assignments.
		it = it+1; % increase the iteration counter
		ndif = sum(tmp_assign ~= assign);
	end
	prwaitbar(0);
	
	if it>=maxit
		prwarning(1,['No convergence reached before the maximum number of %d iterations passed. ' ...
		'The last result was returned.'], maxit);
	end
	
return;
