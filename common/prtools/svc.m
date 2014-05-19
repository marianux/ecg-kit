%SVC Support Vector Classifier
% 
% 	[W,J] = SVC(A,KERNEL,C)
% 	[W,J] = SVC(A,TYPE,PAR,C)
%    W    = A*SVC([],KERNEL,C)
%    W    = A*SVC([],TYPE,PAR,C)
%
% INPUT
%   A	      Dataset
%   KERNEL  - Untrained mapping to compute kernel by A*(A*KERNEL) during
%             training, or B*(A*KERNEL) during testing with dataset B.
%           - String to compute kernel matrices by FEVAL(KERNEL,B,A)
%           Default: linear kernel (PROXM([],'p',1));
%   TYPE    Kernel type (see PROXM)
%   PAR     Kernel parameter (see PROXM)
%   C       Regularization parameter (optional; default: 1)
%
% OUTPUT
%   W       Mapping: Support Vector Classifier
%   J       Object indices of support objects	
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by quadratic
% programming. The non-linearity is determined by the kernel.
% If KERNEL = 0 it is assumed that A is already the kernelmatrix (square).
% In this case also a kernel matrix B should be supplied at evaluation by 
% B*W or PRMAP(B,W).
%
% There are several ways to define KERNEL, e.g. PROXM([],'r',1) for a
% radial basis kernel or by USERKERNEL for a user defined kernel.
%
% If C is NaN this regularisation parameter is optimised by REGOPTC.
%
% See for more possibilties SVCINFO
%
% SEE ALSO
% MAPPINGS, DATASETS, PROXM, USERKERNEL, NUSVC, RBSVC, REGOPTC

% Copyright: D. de Ridder, D. Tax, S. Verzakov, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: svc.m,v 1.11 2009/07/22 19:29:41 duin Exp $

function [W,J,C,nu,alginf] = svc(a,kernel,par1,par2)
		if nargin > 1 & ~isempty(kernel) & isstr(kernel) & exist(kernel) ~= 2 
	% old call: svc(a,type,par,C)
	Options = [];
	if nargin == 4
		C = par2;
	end
	if nargin < 4 | isempty(par2)
		C = 1;
		prwarning(3,'Regularization parameter C set to 1\n');
	end
	if nargin < 3 | isempty(par1)
		par1 = 1;
		prwarning(3,'Kernel parameter par set to 1\n');
	end
	kernel = proxm([],kernel,par1);
else % new call
  if nargin < 4
    Options = [];
	else
		Options = par2;
  end
	if nargin >= 3
		C = par1;
	end
	if nargin < 3 | isempty(par1)
		C = 1;
	end
	if nargin < 2 | isempty(kernel)
		kernel = proxm([],'p',1);
	end
end

DefOptions.mean_centering   = 1;
DefOptions.pd_check         = 1;
DefOptions.bias_in_admreg   = 1;
DefOptions.pf_on_failure    = 1;
DefOptions.multiclass_mode  = 'single';

Options = updstruct(DefOptions,Options,1);

if nargin < 1 | isempty(a)
	W = prmapping(mfilename,{kernel,C});
	W = setname(W,'SVM');
			
elseif (~ismapping(kernel) | isuntrained(kernel)) % training
	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a,'objects');
	[m,k,c] = getsize(a);
	nlab = getnlab(a);
	
	% The SVC is basically a 2-class classifier. More classes are
	% handled by mclassc.
	
	if c == 2   % two-class classifier
	
		if (isnan(C))|length(C)>1     % optimize trade-off parameter
			defs = {proxm([],'p',1),1,[]};
			if isnan(C) % use a default range of C values
				parmin_max = [0,0;1e-2,1e2;0,0];   % kernel and Options can not be optimised
			else % use the range that is supplied by the user
				parmin_max = [0 0; min(C) max(C); 0 0];
				C = NaN;
			end
			[W,J,C,nu,alginf] = regoptc(a,mfilename,{kernel,C,Options},defs,[2],parmin_max,testc([],'soft'));
		else
			% Compute the parameters for the optimization:
			y = 3 - 2*nlab;
			if ~isequal(kernel,0)
      	if Options.mean_centering
        	u = mean(a);         % shift origin for better accuracy
        	b = a - repmat(u,[m,1]);
      	else
        	b = a;
        	u = [];
      	end
				K = b*(b*kernel);
				K = min(K,K'); % make sure kernel matrix is symmetric
				% Perform the optimization:
				[v,J,C,nu] = svo(+K,y,C,Options);
				s = b(J,:);
				insize = k;
			else % kernel is already given!
				K = min(a,a'); % make sure kernel matrix is symmetric
				% Perform the optimization:
				[v,J,C,nu] = svo(+K,y,C,Options);
				s = [];
				u = [];
				insize = size(K,2); % ready for kernel inputs
			end
			% Store the results:
			W = prmapping(mfilename,'trained',{u,s,v,kernel,J},getlablist(a),insize,2);
			% Note: even in case kernel is a mapping, we store the untrained
			% mapping, as training is just administration
  		W = cnormc(W,a);
			W = setcost(W,a);
			
      alginf.svc_type = 'C-SVM';
      alginf.kernel = kernel;
      alginf.C   = C;
      alginf.nu  = nu;
      alginf.nSV = length(J);
      alginf.classsizes = [nnz(y==1), nnz(y==-1)];
      alginf.pf = isnan(nu);
			W = setname(W,'SVM');

		end
		
	else   % multi-class classifier:
		
		[W,J,C,nu,alginf] = mclassc(a,prmapping(mfilename,{kernel,C,Options}),Options.multiclass_mode);
		
	end
	W = setname(W,'SVM');

else % execution
	
	v = kernel;
	w = +v;
	m = size(a,1);
	
	if isempty(w{2})
		K = a; % user supplied testkernel
		J = w{5};
		if size(K,2) > length(J) & size(K,2) >= max(J)
			K = K(:,J); % probably full test kernel
		elseif size(K,2) ~= length(J)
			error('Wrong test kernel supplied')
		end			
	else
		% The first parameter w{1} stores the mean of the dataset. When it
		% is supplied, remove it from the dataset to improve the numerical
		% precision. Then compute the kernel matrix using proxm:
		if isempty(w{1})
			K = a*(w{2}*w{4});
		else
			K = (a-ones(m,1)*w{1})*(w{2}*w{4});
		end
	end
	% Data is mapped by the kernel, now we just have a linear
	% classifier  w*x+b:
	d = [K ones(m,1)] * w{3};
	W = setdat(a,[d -d],v);
	
end

return
