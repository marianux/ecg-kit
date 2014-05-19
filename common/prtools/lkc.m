%LKC Linear kernel classifier
% 
% 	W = LKC(A,KERNEL)
%
% INPUT
%   A	      Dataset
%   KERNEL  Mapping to compute kernel by A*MAP(A,KERNEL)
%           or string to compute kernel by FEVAL(KERNEL,A,A)
%           or cell array with strings and parameters to compute kernel by
%           FEVAL(KERNEL{1},A,A,KERNEL{2:END})
%           Default: linear kernel (PROXM([],'P',1))
%
% OUTPUT
%   W       Mapping: Support Vector Classifier
%
% DESCRIPTION
% This is a fall-back routine for other kernel procedures like SVC, RBSVC
% and LIBSVC. If they fail due to optimization problems they may fall back
% to this routine which computes a linear classifier in kernelspace using
% the pseudo-inverse of the kernel.
% 
% The kernel may be supplied in KERNEL by
% - an untrained mapping, e.g. a call to PROXM like W = LIBSVC(A,PROXM([],'R',1))
% - a string with the name of the routine to compute the kernel from A
% - a cell-array with this name and additional parameters.
% This will be used for the evaluation of a dataset B by B*W or PRMAP(B,W) as
% well. 
%
% If KERNEL = 0 it is assumed that A is already the kernel matrix (square).
% In this also a kernel matrix should be supplied at evaluation by B*W or 
% PRMAP(B,W). 
%
% SEE ALSO 
% MAPPINGS, DATASETS, SVC, PROXM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
function W = lkc(a,kernel)

		
	if nargin < 2 | isempty(kernel)
		kernel = proxm([],'p',1);
	end

	if nargin < 1 | isempty(a)
		W = prmapping(mfilename,{kernel});
		W = setname(W,'LKC Classifier');
		return;
	end

	if (~ismapping(kernel) | isuntrained(kernel)) % training
	
		islabtype(a,'crisp');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a,'objects');
		[m,k,c] = getsize(a);
		nlab = getnlab(a); 
	
		K = compute_kernel(a,a,kernel);
		K = min(K,K');   % make sure kernel is symmetric
		targets = gettargets(setlabtype(a,'targets'));
		v = prpinv([K ones(m,1); ones(1,m) 0])*[targets; zeros(1,c)];
		
		lablist = getlablist(a);         
		W = prmapping(mfilename,'trained',{v,a,kernel},lablist,size(a,2),c);
		W = setname(W,'LKC Classifier');
  	W = cnormc(W,a);
		W = setcost(W,a);
		
	else % execution
		w = kernel;
		v = getdata(w,1);   % weights
		s = getdata(w,2);   % trainset or empty
		kernel = getdata(w,3); % kernelfunction or 0
		m = size(a,1);
		
		K = compute_kernel(a,s,kernel); % kernel testset
		% Data is mapped by the kernel, now we just have a linear
		% classifier  w*x+b:
		d = [K ones(m,1)]*v;
		if size(d,2) == 1, d = [d  -d]; end
		W = setdat(a,d,w);
	end
	
return;

function K = compute_kernel(a,s,kernel)

	% compute a kernel matrix for the objects a w.r.t. the support objects s
	% given a kernel description

	if  isstr(kernel) % routine supplied to compute kernel
		K = feval(kernel,a,s);
	elseif iscell(kernel)
		K = feval(kernel{1},a,s,kernel{2:end});
	elseif ismapping(kernel)
		K = a*prmap(s,kernel);
	elseif kernel == 0 % we have already a kernel
		K = a;
	else
		error('Do not know how to compute kernel matrix')
	end
		
	K = +K;
		
return
