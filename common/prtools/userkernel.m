%USERKERNEL Construct user defined kernel mapping  
%
%  K = USERKERNEL(B,R,FUNC,P1,P2, ...)
%  K = B*USERKERNEL([],R,FUNC,P1,P2, ...)
%  W = USERKERNEL([],R,FUNC,P1,P2, ...)
%  W = R*USERKERNEL([],[],FUNC,P1,P2, ...)
%  K = B*W
%
% INPUT
%   R     Dataset, representation set, default B
%   B     Dataset
%   FUNC  String with function name to compute kernels (proximities) by
%         K = FEVAL(FUNC,B,R,P1,P2, ...) 
%
% OUTPUT
%   W     Trained kernel mapping
%   K     Kernel (proximity) matrix
%
% DESCRIPTION
% The kernel matrix K is computed according to the definition given by the
% user supplied function FUNC. The size of K is [SIZE(B,1) SIZE(R,1)]

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = userkernel(a,r,kernel,varargin)
  		
	
	if nargin < 3, kernel = []; end
	if nargin < 2, r = []; end
	if nargin < 1, a = []; end
	
	if isempty(a) & isempty(r)
		w = prmapping(mfilename,'untrained',{[],kernel,varargin{:}});
		w = setname(w,'userkernel mapping');
	elseif isempty(r)
		[m,k] = size(a);
		w = prmapping(mfilename,'trained',{a,kernel,varargin},getlab(a),k,m);
	elseif isempty(a)
		[m,k] = size(r);
		w = prmapping(mfilename,'trained',{r,kernel,varargin},getlab(r),k,m);
	elseif ismapping(r) % execution of a*userkernel (trained)
		u = getdata(r);
		[r,kernel,pars] = deal(u{:});
		w = compute_kernel(kernel,a,r,pars);
	elseif isdataset(a) % execution
		w = compute_kernel(kernel,a,r,varargin);
	end
		
	function k = compute_kernel(kernel,a,r,pars)
	
		if isempty(kernel)
			error('No kernel function or kernel mapping supplied')
		elseif ~exist(kernel)
			error('kernel function not found')
		end
		
		if isempty(pars)
			k = feval(kernel,+a,+r);
		else
			k = feval(kernel,+a,+r,pars{:});
		end
		k = setdat(a,k);
		if isdataset(r)
			k = setfeatlab(k,getlabels(r));
		end
	return
	