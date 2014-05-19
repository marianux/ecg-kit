%SVO_NU Support Vector Optimizer: NU algorithm
%
%   [V,J,C] = SVO(K,NLAB,NU,PD)
%
% INPUT
%   K     Similarity matrix
%   NLAB  Label list consisting of -1/+1
%   NU    Regularization parameter (0 < NU < 1): expected fraction of SV (optional; default: 0.25)
%
%   PD    Do or do not the check of the positive definiteness (optional; default: 1 (to do))
%
% OUTPUT
%   V     Vector of weights for the support vectors
%   J     Index vector pointing to the support vectors
%   C     Equivalent C regularization parameter of SVM-C algorithm
%
% DESCRIPTION
% A low level routine that optimizes the set of support vectors for a 2-class
% classification problem based on the similarity matrix K computed from the
% training set. SVO is called directly from SVC. The labels NLAB should indicate 
% the two classes by +1 and -1. Optimization is done by a quadratic programming. 
% If available, the QLD function is used, otherwise an appropriate Matlab routine.
%
% NU is bounded from above by NU_MAX = (1 - ABS(Lp-Lm)/(Lp+Lm)), where
% Lp (Lm) is the number of positive (negative) samples. If NU > NU_MAX is supplied 
% to the routine it will be changed to the NU_MAX.
%
% If NU is less than some NU_MIN which depends on the overlap between classes 
% algorithm will typically take long time to converge (if at all). 
% So, it is advisable to set NU larger than expected overlap.
%
% Weights V are rescaled in a such manner as if they were returned by SVO with the parameter C.
%
% SEE ALSO
% SVC_NU, SVO, SVC

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Based on SVO.M by D.M.J. Tax, D. de Ridder, R.P.W. Duin
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: svo_nu.m,v 1.3 2010/02/08 15:29:48 duin Exp $

function [v,J,C] = svo_nu(K,y,nu,pd)
	  if (nargin < 4) 
    pd = 1;
	end

	if (nargin < 3)
		prwarning(3,'Third parameter (nu) not specified, assuming 0.25.');
    nu = 0.25;
	end

  nu_max = 1 - abs(nnz(y == 1) - nnz(y == -1))/length(y);

  if nu > nu_max
    prwarning(3,['nu==' num2str(nu)  ' is not feasible; set to ' num2str(nu_max)]);
    nu = nu_max;
  end 

  vmin = 1e-9;		% Accuracy to determine when an object becomes the support object.

	% Set up the variables for the optimization.
	n  = size(K,1);
	D  = (y*y').*K;
	f  = zeros(1,n);
	A  = [y';ones(1,n)];
	b  = [0; nu*n];
	lb = zeros(n,1);
	ub = ones(n,1);
	p  = rand(n,1);

  if pd
  	% Make the kernel matrix K positive definite.
	  i = -30;
  	while (pd_check (D + (10.0^i) * eye(n)) == 0)
	  	i = i + 1;
  	end
  	if (i > -30),
	  	prwarning(2,'K is not positive definite. The diagonal is regularized by 10.0^(%d)*I',i);
  	end
	  i = i+2;
  	D = D + (10.0^(i)) * eye(n);
	end
   
  % Minimization procedure initialization:
	% 'qp' minimizes:   0.5 x' D x + f' x
	% subject to:       Ax <= b
	%
	if (exist('qld') == 3)
		v = qld (D, f, -A, b, lb, ub, p, length(b));
	elseif (exist('quadprog') == 2)
		prwarning(1,'QLD not found, the Matlab routine QUADPROG is used instead.')
		v = quadprog(D, f, [], [], A, b, lb, ub);
	else 
		prwarning(1,'QLD not found, the Matlab routine QP is used instead.')
		verbos = 0;
		negdef = 0;
		normalize = 1;
		v = qp(D, f, A, b, lb, ub, p, length(b), verbos, negdef, normalize);
	end
	
	% Find all the support vectors.
	J = find(v > vmin);

	% First find the SV on the boundary
	I  = J(v(J) < 1-vmin);
	Ip = I(y(I) ==  1);
	Im = I(y(I) == -1);
	if (isempty(v) | isempty(Ip) | isempty(Im))
		%error('Quadratic Optimization failed. Pseudo-Fisher is computed instead.');
		prwarning(1,'Quadratic Optimization failed. Pseudo-Fisher is computed instead.');
		v = prpinv([K ones(n,1)])*y;
		J = [1:n]';
		C = nan;
    return;
	end

  v = y.*v;

  %wxI = K(I,J)*v(J);
  wxIp = mean(K(Ip,J)*v(J),1); %  rho-b
  wxIm = mean(K(Im,J)*v(J),1); % -rho-b

  rho =  0.5*(wxIp-wxIm);
	b   = -0.5*(wxIp+wxIm);
  %b = mean(rho*y(I) - wxI);

	v = [v(J); b]/rho;
  C = 1/rho;   

  return;

