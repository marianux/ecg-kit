%SVMR SVM regression
%
%      W = SVMR(X,NU,KTYPE,KPAR,EP)
%      W = X*SVMR([],NU,KTYPE,KPAR,EP)
%      W = X*SVMR(NU,KTYPE,KPAR,EP)
%
% INPUT
%   X      Regression dataset
%   NU     Fraction of objects outside the 'data tube'
%   KTYPE  Kernel type (default KTYPE='p', for polynomial)
%   KPAR   Extra parameter for the kernel
%   EP     Epsilon, with of the 'data tube'
%
% OUTPUT
%   W      Support vector regression
%
% DESCRIPTION
% Train an nu-Support Vector Regression on dataset X with parameter NU.
% The kernel is defined by kernel type KTYPE and kernel parameter KPAR.
% For the definitions of these kernels, have a look at proxm.m. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%  LINEARR, PROXM, GPR, SVC

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function y = svmr(varargin)

	argin = shiftargin(varargin,{'scalar','char'},1);
	argin = shiftargin(argin,'char',2);
  argin = setdefaults(argin,[],0.01,'p',1,0.2);
  
  if mapping_task(argin,'definition')
    y = define_mapping(argin,'untrained');
    y = setname(y,'SVM regression');
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [x,nu,ktype,kpar,ep] = deal(argin{:});
    [n,d] = size(x);
    y = gettargets(x);
    % kernel mapping:
    wk = proxm(+x,ktype,kpar);
    K = +(x*wk);
    % setup optimization:
    tol = 1e-6;
    C = 1/(n*nu);
    H = [K -K; -K K];
    f = repmat(ep,2*n,1) - [y; -y];
    A = [];
    b = [];
    Aeq = [ones(1,n) -ones(1,n)];
    beq = 0;
    lb = zeros(2*n,1);
    ub = repmat(C,2*n,1);
    %do the optimization:
    if exist('qld')
      alf = qld(H,f,Aeq,beq,lb,ub,[],1);
    else
      alf = quadprog(H,f,A,b,Aeq,beq,lb,ub);
    end
    % find SV's
    w = alf(1:n)-alf((n+1):end);  % the real weights
    Isv = find(abs(w)>tol);
    Ibnd = find(abs(w(Isv))<C-tol);
    % for the 'real' sv's, we want the classifier output to compute the
    % offset:
    I = Isv(Ibnd);
    out1 = sum(repmat(w',length(Ibnd),1).*K(I,:),2);
    sw = sign(w);
    out = y(I) - out1 - sw(I).*ep;

    % store the useful parameters:
    W.b = mean(out);
    W.wk = proxm(+x(Isv,:),ktype,kpar);
    W.w = w(Isv);

    y = prmapping(mfilename,'trained',W,1,d,1);
    y = setname(y,'Linear regression');
  else
    % evaluation
    [x,v] = deal(argin{1:2});
    W = getdata(v);
    [n,d] = size(x);
    Kz = +(x*W.wk);
    out = sum(repmat(W.w',n,1).*Kz,2) + repmat(W.b,n,1);
    y = setdat(x,out);

  end
  
return
