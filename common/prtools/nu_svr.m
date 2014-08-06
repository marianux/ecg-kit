%NU_SVR Support Vector Classifier: NU algorithm
% 
% 	[W,J,C] = NU_SVR(A,TYPE,PAR,C,SVR_TYPE,NU_EPS,MC,PD)
%
% INPUT
%   A	    Dataset
%   TYPE  Type of the kernel (optional; default: 'p')
%   PAR   Kernel parameter (optional; default: 1)
%   C     Regularization parameter (0 < C < 1): expected fraction of SV
%         (optional; default: 0.25)
%   SVR_TYPE  This type can be 'nu' or 'epsilon'
%   NU_EPS The corresponding value for NU or epsilon
%   MC    Do or do not data mean-centering (optional; default: 1 (to do))
%   PD    Do or do not the check of the positive definiteness (optional;
%         default: 1 (to do))
%
% OUTPUT
%   W     Mapping: Support Vector Classifier
%   J     Object identifiers of support objects		
%   C     Equivalent C regularization parameter of SVM-C algorithm
%
% DESCRIPTION
% Optimizes a support vector classifier for the dataset A by 
% quadratic programming. The classifier can be of one of the types 
% as defined by PROXM. Default is linear (TYPE = 'p', PAR = 1). In J 
% the identifiers of the support objects in A are returned.
%
% C belogs to the interval (0,1). C close to 1 allows for more class
% overlap.  Default C = 0.25.
% 
% C is bounded from above by NU_MAX = (1 - ABS(Lp-Lm)/(Lp+Lm)), where
% Lp (Lm) is the number of positive (negative) samples. If NU > NU_MAX
% is supplied to the routine it will be changed to the NU_MAX.
%
% If C is less than some NU_MIN which depends on the overlap between
% classes algorithm will typically take long time to converge (if at
% all).  So, it is advisable to set NU larger than expected overlap.
%
% Output is rescaled in a such manner as if it were returned by SVC with
% the parameter C.
%
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% NU_SVRO, SVO, SVC, MAPPINGS, DATASETS, PROXM

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Based on SVC.M by D.M.J. Tax, D. de Ridder, R.P.W. Duin
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: nu_svr.m,v 1.2 2009/01/31 18:43:11 duin Exp $

function [W, J, epsilon_or_nu] = nu_svcr(a,type,par,C,svr_type,nu_or_epsilon,mc,pd)
if nargin < 2 | ~isa(type,'prmapping')
	if nargin < 8
		pd = 1;  
	end
	if nargin < 7
		mc = 1;  
	end
	if nargin < 6 
		nu_or_epsilon = [];
	end
	if nargin < 5 | isempty(svr_type)
		svr_type = 'epsilon';
	end

	switch svr_type
	case 'nu'
		if isempty(nu_or_epsilon)
			prwarning(3,'nu is not specified, assuming 0.25.');
			nu_or_epsilon = 0.25;
		end
	%nu = nu_or_epsilon;
	case {'eps', 'epsilon'}
		svr_type = 'epsilon';
		if isempty(nu_or_epsilon)
			prwarning(3,'epsilon is not specified, assuming 1e-2.');
			nu_or_epsilon = 1e-2;
		end
	%epsilon = nu_or_epsilon;
	end

	if nargin < 4 | isempty(C)
		prwarning(3,'C set to 1\n');
		C = 1;
	end

	if nargin < 3 | isempty(par)
		par = 1;
		prwarning(3,'Kernel parameter par set to 1\n');
	end
	if nargin < 2 | isempty(type)
		type = 'p';
		prwarning(3,'Polynomial kernel type is used\n');
	end
	if nargin < 1 | isempty(a)
		W = prmapping(mfilename,{type,par,C,svr_type,nu_or_epsilon,mc,pd});
		W = setname(W,['Support Vector Regression (' svr_type ' algorithm)']);
		return;
	end

	islabtype(a,'targets');
	[m,k] = getsize(a);
	y = gettargets(a);	
	% The 1-dim SVR
	if size(y,2) == 1   % 1-dim regression
		uy = mean(y);
		y  = y - uy;
		if mc
			u  = mean(a);
			a  = a - ones(m,1)*u;
		else
			u  = [];
		end

		K = a*proxm(a,type,par);
		% Perform the optimization:
		[v,J,epsilon_or_nu] = nu_svro(+K,y,C,svr_type,nu_or_epsilon,pd);
		% Store the results:
		v(end) = v(end)+uy; 
		W = prmapping(mfilename,'trained',{u,a(J,:),v,type,par},getlablist(a),k,1);
		W = setname(W,['Support Vector Regression (' svr_type ' algorithm)']);
		%W = setcost(W,a);
		J = getident(a,J);
		%J = a.ident(J);

	else   
		error('multivariate SVR is not supported');
	end

else % execution
	w = +type;
	m = size(a,1);

	% The first parameter w{1} stores the mean of the dataset. When it
	% is supplied, remove it from the dataset to improve the numerical
	% precision. Then compute the kernel matrix using proxm.

	if isempty(w{1})
		d = a*proxm(w{2},w{4},w{5});
	else
		d = (a-ones(m,1)*w{1})*proxm(w{2},w{4},w{5});
	end

	% When Data is mapped by the kernel, now we just have a linear
	% regression  w*x+b:
	d = [d ones(m,1)] * w{3};
	W = setdat(a,d,type);
end
	
return;
