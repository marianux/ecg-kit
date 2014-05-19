%GPR Gaussian Process regression
%
%     W = GPR(A,KERNEL,S_noise)
%
%INPUT
%  A        Dataset
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
%SEE ALSO
% svmr, proxm, linearr, testr, plotr

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function y = gpr(x,kernel,s_noise)

if nargin<3
	s_noise = 1;
end
if nargin<2
	kernel = proxm([],'p',1);
end
if nargin<1 || isempty(x)
	y = prmapping(mfilename,{kernel,s_noise});
	y = setname(y,'Gaussian Proc. regression');
	return
end

if ~ismapping(kernel) || ~istrained(kernel) %training
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
else
	% evaluation
	W = getdata(kernel);
    out = [+x ones(size(x,1),1)]*W.K*W.w;
    y = setdat(x,out);
	
end
