%SVC Trainable classifier: Support Vector Machine
% 
% 	[W,J] = SVC(A,KERNEL,C)
%   [W,J] = A*SVC([],KERNEL,C)
%   [W,J] = A*SVC(KERNEL,C)
%
% INPUT
%   A	      Dataset
%   KERNEL  - Untrained mapping to compute kernel by A*(A*KERNEL) during
%             training, or B*(A*KERNEL) during testing with dataset B.
%           - String to compute kernel matrices by FEVAL(KERNEL,B,A)
%           Default: linear kernel (PROXM('p',1))
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
% There are several ways to define KERNEL, e.g. PROXM('r',1) for a
% radial basis kernel or by USERKERNEL for a user defined kernel.
%
% If C is NaN this regularisation parameter is optimised by REGOPTC.
%
% See for more possibilties SVCINFO
%
% EXAMPLE
% a = gendatb;                     % generate banana classes
% [w,J] = a*svc(proxm('p',3));  % compute svm with 3rd order polynomial
% a*w*testc                        % show error on train set
% scatterd(a)                      % show scatterplot
% plotc(w)                         % plot classifier
% hold on; 
% scatterd(a(J,:),'o')             % show support objcts
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PROXM, USERKERNEL, NUSVC, RBSVC, LIBSVC, REGOPTC

% Copyright: D. de Ridder, D. Tax, S. Verzakov, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [W,J,C,nu,alginf] = svc(varargin)

  mapname = 'SVC';
  argin = shiftargin(varargin,{'prmapping','char','cell'});
  argin = setdefaults(argin,[],proxm([],'p',1),1,1);
  
  if mapping_task(argin,'definition')
    
    W = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.

    [a,kernel,C,Options] = check_for_old_call(argin);
    if isequal(Options,1), Options = []; end

    DefOptions.mean_centering   = 1;
    DefOptions.pd_check         = 1;
    DefOptions.bias_in_admreg   = 1;
    DefOptions.pf_on_failure    = 1;
    DefOptions.multiclass_mode  = 'single';

    Options = updstruct(DefOptions,Options,1);

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
        if isequal(kernel,0)
          s = [];
          u = [];
          in_size = 0; % to allow old and new style calls
        else
          if Options.mean_centering
            u = mean(a);         % shift origin for better accuracy
            a = a - repmat(u,[m,1]);
          else
            u = [];
          end
          in_size = k;
        end
        K = compute_kernel(a,a,kernel);
        K = min(K,K');   % make sure kernel is symmetric
        [v,J,C,nu] = svo(+K,y,C,Options);
        
        % Store the results:
        if ~isequal(kernel,0)
          s = a(J,:);
        end   
        W = prmapping(mfilename,'trained',{u,s,v,kernel,J},getlablist(a),in_size,2);
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
      W = W*classc;

    end
    W = setname(W,mapname);

  else % Evaluation

    [a,v] = deal(argin{1:2});
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

function K = compute_kernel(a,s,kernel)

	% compute a kernel matrix for the objects a w.r.t. the support objects s
	% given a kernel description

	if  ischar(kernel) % routine supplied to compute kernel
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

function [a,kernel,C,par] = check_for_old_call(argin)

  [a,kernel,C,par] = deal(argin{:});
  if ischar(kernel) && exist(kernel,'file') ~= 2
    kernel = proxm(kernel,C);
    C = par;
    par = [];
  end
  
return
