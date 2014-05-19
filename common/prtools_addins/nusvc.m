%NUSVC Support Vector Classifier: NU algorithm
% 
% 	[W,J,NU] = NUSVC(A,KERNEL,NU)
%    W    = A*SVC([],KERNEL,NU)
%
% INPUT
%   A	      Dataset
%   KERNEL  - Untrained mapping to compute kernel by A*(A*KERNEL) during
%             training, or B*(A*KERNEL) during testing with dataset B.
%           - String to compute kernel matrices by FEVAL(KERNEL,B,A)
%           Default: linear kernel (PROXM([],'p',1));
%   NU      Regularization parameter (0 < NU < 1): expected fraction of SV 
%           (optional; default: max(leave-one-out 1_NN error,0.01))
%
% OUTPUT
%   W       Mapping: Support Vector Classifier
%   J       Object indices of support objects		
%   NU      Actual nu_value used
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by quadratic 
% programming. The difference with the standard SVC routine is the use and
% interpretation of the regularisation parameter NU. It is an upperbound
% for the expected classification error. By default NU is estimated by the
% leave-one-error of the 1_NN rule. For NU = NaN an automatic optimisation 
% is performed using REGOPTC. 
%
% If KERNEL = 0 it is assumed that A is already the kernelmatrix (square).
% In this case also a kernel matrix B should be supplied at evaluation by
% B*W or MAP(B,W).
%
% There are several ways to define KERNEL, e.g. PROXM([],'r',1) for a
% radial basis kernel or by USERKERNEL for a user defined kernel.
%
% SEE ALSO
% MAPPINGS, DATASETS, SVC, NUSVO, PROXM, USERKERNEL, REGOPTC

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Based on SVC byby: D. de Ridder, D. Tax, R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: nusvc.m,v 1.6 2010/06/25 09:50:40 duin Exp $

function [W,J,nu,C,alginf] = nusvc(a,kernel,nu,Options)

	prtrace(mfilename);
		
  if nargin < 4
    Options = [];
  end

  DefOptions.mean_centering       = 1;
  DefOptions.pd_check             = 1;
  DefOptions.bias_in_admreg       = 1;
  DefOptions.allow_ub_bias_admreg = 1;
  DefOptions.pf_on_failure        = 1;
  DefOptions.multiclass_mode      = 'single';    

  Options = updstruct(DefOptions,Options,1);
		
  if nargin < 3
    nu = [];
	  prwarning(3,'Regularization parameter nu set to NN error\n');
  end
	if nargin < 2 | isempty(kernel)
		kernel = proxm([],'p',1);
	end
	
	if nargin < 1 | isempty(a)
    W = prmapping(mfilename,{kernel,nu,Options});
    W = setname(W,'nuSVM');

	elseif (~ismapping(kernel) | isuntrained(kernel)) % training
		
		pd = 1; % switches, fixed here.
		mc = 1;
  	islabtype(a,'crisp');
  	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a,'objects');
  	[m,k,c] = getsize(a);
  	nlab = getnlab(a);  
	
	%	if isempty(nu), nu = max(testk(a,1),0.01); end
		if isempty(nu)
			nu = 2*min(max(testk(a,1),0.01),(0.8*min(classsizes(a))/size(a,1)));
		end

		% The SVC is basically a 2-class classifier. More classes are
  	% handled by mclassc.
    if c == 2   % two-class classifier
			if (isnan(nu))                      % optimize trade-off parameter
				defs = {proxm([],'p',1),[],[]};
				parmin_max = [0,0;0.001,0.999;0,0]; % kernel and Options can not be optimised
				[W,J,nu,C,alginf] = regoptc(a,mfilename,{kernel,nu,Options},defs,[2],parmin_max,testc([],'soft'));
			else
      	% Compute the parameters for the optimization:
      	y = 3 - 2*nlab;
			  if ~isequal(kernel,0)
          if mc
            u = nanmean(a);
            b = a -ones(m,1)*u;
          else
            b = a;
            u = [];
					end
          K = b*(b*kernel);
         	% Perform the optimization:
				  [v,J,nu,C] = nusvo(+K,y,nu,Options);
          s = b(J,:);
          insize = k;
				else % kernel is already given!
					K = min(a,a'); % make sure kernel matrix is symmetric
					% Perform the optimization:
					[v,J,nu,C] = nusvo(+K,y,nu,Options);
					u = [];
            s = [];
					insize = size(K,2); % ready for kernel inputs
				end      
        % Store the results, use SVC for execution
      	W = prmapping('svc','trained',{u,s,v,kernel,J},getlablist(a),insize,2);
      	W = cnormc(W,a);
      	W = setcost(W,a);
      	alginf.svc_type = 'nu-SVM';
      	alginf.kernel = kernel;
      	alginf.C   = C;
      	alginf.nu  = nu;
      	alginf.nSV = length(J);
      	alginf.classsizes = [nnz(y==1), nnz(y==-1)];
      	alginf.pf = isnan(C);
			end
  		
  	else   % multi-class classifier:
	    [W,J,nu,C,alginf] = mclassc(a,prmapping(mfilename,{kernel,nu,Options}),'single');
		end
		
 		W = setname(W,'nuSVM');

	end
	
  return;
