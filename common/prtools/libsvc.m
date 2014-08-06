%LIBSVC Trainable classifier: LIBSVM
% 
% 	[W,J] = LIBSVC(A,KERNEL,C)
% 	[W,J] = A*LIBSVC([],KERNEL,C)
% 	[W,J] = A*LIBSVC(KERNEL,C)
%
% INPUT
%   A	    Dataset
%   KERNEL  Mapping to compute kernel by A*MAP(A,KERNEL)
%           or string to compute kernel by FEVAL(KERNEL,A,A)
%           or cell array with strings and parameters to compute kernel by
%           FEVAL(KERNEL{1},A,A,KERNEL{2:END})
%           Default: linear kernel (PROXM([],'P',1))
%   C       Trade_off parameter in the support vector classifier.
%           Default C = 1;
%
% OUTPUT
%   W       Mapping: Support Vector Classifier
%   J       Object idences of support objects. Can be also obtained as W{4}		
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by the libsvm
% package, see http://www.csie.ntu.edu.tw/~cjlin/libsvm/. LIBSVC calls the
% svmtrain routine of libsvm for training. Classifier execution for a
% test dataset B may be done by D = B*W; In D posterior probabilities are
% given as computed by svmpredict using the '-b 1' option. 
% 
% The kernel may be supplied in KERNEL by
% - an untrained mapping, e.g. a call to PROXM like W = LIBSVC(A,PROXM([],'R',1))
% - a string with the name of the routine to compute the kernel from A
% - a cell-array with this name and additional parameters.
% This will be used for the evaluation of a dataset B by B*W or PRMAP(B,W) as
% well. 
%
% If KERNEL = 0 (or not given) it is assumed that A is already the 
% kernelmatrix (square). In this also a kernel matrix should be supplied at 
% evaluation by B*W or PRMAP(B,W). However, the kernel has to be computed with 
% respect to support objects listed in J (the order of objects in J does matter).
%
% EXAMPLE
% a = gendatb;                     % generate banana classes
% [w,J] = a*libsvc(proxm('p',3));  % compute svm with 3rd order polynomial
% a*w*testc                        % show error on train set
% scatterd(a)                      % show scatterplot
% plotc(w)                         % plot classifier
% hold on; 
% scatterd(a(J,:),'o')             % show support objcts
% 
% REFERENCES
% R.-E. Fan, P.-H. Chen, and C.-J. Lin. Working set selection using the second order 
% information for training SVM. Journal of Machine Learning Research 6, 1889-1918, 2005
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% MAPPINGS, DATASETS, SVC, PROXM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
function [W,J,u] = libsvc(varargin)
		
	checktoolbox('libsvm');

  mapname = 'LIBSVM';
  argin = shiftargin(varargin,{'prmapping','char','cell'});
  argin = setdefaults(argin,[],proxm([],'p',1),1,1);
  
  if mapping_task(argin,'definition')
    
    W = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.

    [a,kernel,C] = check_for_old_call(argin);
    opt = ['-s 0 -t 4 -b 1 -e 1e-3 -c ',num2str(C), ' -q'];
    

    islabtype(a,'crisp');
    isvaldfile(a,1,2); % at least 1 object per class, 2 classes
    a = testdatasize(a,'objects');
    [m,k,c] = getsize(a);
    nlab = getnlab(a); 

    K = compute_kernel(a,a,kernel);
    K = min(K,K');   % make sure kernel is symmetric
    K = [[1:m]' K];  % as libsvm wants it
                     % call libsvm
    u = svmtrain(nlab,K,opt);
		if isempty(u)
			prwarning(1,'libsvc: no solution for SVM, pseudo-inverse will be used')
			W = lkc(prdataset(K(:,2:end),getlabels(a)),0);
			J = [1:m]';
			return
		end
                     % Store the results:
    J = full(u.SVs);
    if isempty(J) | J == 0
      % LIBSVM failed, use fisher
      W = fisherc(a);
      prwarning(1,'LIBSVC failed, Fisher used instead')
      return
    end
    if isequal(kernel,0)
      s = [];
      in_size = 0; % to allow old and new style calls
    else
      s = a(J,:);
      in_size = k;
    end

    lablist = getlablist(a);         
    W = prmapping(mfilename,'trained',{u,s,kernel,J,opt},lablist(u.Label,:),in_size,c);

    W = setname(W,'LIBSVM Classifier');
    W = setcost(W,a);

  else % Evaluation

    [a,W] = deal(argin{1:2});
    [u,s,kernel,J,opt] = getdata(W);
    m = size(a,1);

    K = compute_kernel(a,s,kernel);
    k = size(K,2);
    if k ~= length(J)
      if isequal(kernel,0)
        if (k > length(J)) &  (k >= max(J))
          % precomputed kernel; old style call
          prwarning(2,'Old style execution call: The precomputed kernel was calculated on a test set and the whole training set!')  
        else
          error(['Inappropriate precomputed kernel!' newline ...
              'For the execution the kernel matrix should be computed on a test set' ...
              newline 'and the set of support objects']);
        end  
      else
        error('Kernel matrix has the wrong number of columns');
      end
    else  
      % kernel was computed with respect to the support objects
      % we make an approprite correction in the libsvm structure
      u.SVs = sparse((1:length(J))');
    end  
    K = [[1:m]' K];  % as libsvm wants it
    %[lab,acc,d] = svmpredict(getnlab(a),K,u,' -b 1');
    [lab,acc,d] = svmpredict(ones(m,1),K,u,' -b 1');
    W = setdat(a,d,W);
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

function [a,kernel,C] = check_for_old_call(argin)

[a,kernel,C,par] = deal(argin{:});
if ischar(kernel) && exist(kernel,'file') ~= 2
  kernel = proxm(kernel,C);
  C = par;
end