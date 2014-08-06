%GPR Trainable mapping for Gaussian Process regression
%
%     W = GPR(A,KERNEL,S_noise)
%     W = A*GPR([],KERNEL,S_noise)
%     W = A*GPR(KERNEL,S_noise)
%
%INPUT
%  A        Dataset used for training
%  KERNEL   Untrained mapping to compute kernel by A*(A*KERNEL)
%           during training, or B*(A*KERNEL) during evaluation with
%           dataset B
%  S_noise  Standard deviation of the noise 
%
%OUTPUT
%  W        Mapping: Gaussian Process regression
%
%DESCRIPTION
%Fit a Gaussian Process regressor on dataset A. For a nonlinear regressor,
%define kernel mapping KERNEL. For kernel definitions, have a look at
%proxm.m.
%
%SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%DATASETS, MAPPINGS, SVMR, PROXM, LINEARR, TESTR, PLOTR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = gpr(varargin)

	argin = shiftargin(varargin,'prmapping',1);
	argin = shiftargin(argin,'scalar',2);
  argin = setdefaults(argin,[],proxm([],'p',1),1);
  
  if mapping_task(argin,'definition')
    y = define_mapping(argin,'untrained');
    y = setname(y,'Gaussian Proc. regression');
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [x,kernel,s_noise] = deal(argin{:});
    [n,d] = size(x);
    W.X = [+x ones(n,1)];
  % train:
    W.K = W.X*kernel;
    L = chol(+(W.X*W.K) + s_noise*s_noise*eye(n));
    W.w = L\(L'\gettargets(x));
  % store:
    W.kernel = kernel;
    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'Gaussian Proc. regression');
    
  else                                      % Evaluation
    
    [x,v] = deal(argin{1:2});
    W = getdata(v);
    out = [+x ones(size(x,1),1)]*W.K*W.w;
    y = setdat(x,out);
	
  end
