%SVO Support Vector Optimizer, low-level routine
%
%   [V,J,C,NU] = SVO(K,NLAB,C,OPTIONS)
%
% INPUT
%   K     Similarity matrix
%   NLAB  Label list consisting of -1/+1
%   C     Scalar for weighting the errors (optional; default: 1)
%   OPTIONS
%    .PD_CHECK          force positive definiteness of the kernel by adding a small constant 
%                       to a kernel diagonal (default: 1)
%    .BIAS_IN_ADMREG    it may happen that bias of svc (b term) is not defined, then 
%                       if BIAS_IN_ADMREG == 1, b will be taken from the midpoint of its admissible
%                       region, otherwise (BIAS_IN_ADMREG == 0) the situation will be considered 
%                       as an optimization failure and treated accordingly (deafault: 1)
%    .PF_ON_FAILURE     if optimization is failed (or bias is undefined and BIAS_IN_ADMREG is 0)
%                       and PF_ON_FAILURE == 1, then Pseudo Fisher classifier will be computed, 
%                       otherwise (PF_ON_FAILURE == 0) an error will be issued (default: 1)
%
% OUTPUT
%   V     Vector of weights for the support vectors
%   J     Index vector pointing to the support vectors
%   C     C which was actually used for optimization
%   NU    NU parameter of NUSVC algorithm, which gives the same classifier
%
% DESCRIPTION
% A low level routine that optimizes the set of support vectors for a 2-class
% classification problem based on the similarity matrix K computed from the
% training set. SVO is called directly from SVC. The labels NLAB should indicate 
% the two classes by +1 and -1. Optimization is done by a quadratic programming. 
% If available, the QLD function is used, otherwise an appropriate Matlab routine.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% SVC

% Copyright: D.M.J. Tax, D. de Ridder, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: svo.m,v 1.6 2010/02/08 15:29:48 duin Exp $

function [v,J,C,nu] = svo(K,y,C,Options)

	  if nargin < 4
    Options = [];
  end
  
  DefOptions.pd_check       = 1;
  DefOptions.bias_in_admreg = 1;
  DefOptions.pf_on_failure  = 1;
  
  Options = updstruct(DefOptions, Options,1);

	if nargin < 3 | isempty(C)
		prwarning(3,'The regularization parameter C is not specified, assuming 1.');
    C = 1;
	end

	vmin = 1e-9;		% Accuracy to determine when an object becomes the support object.

  vmin1 = min(1,C)*vmin; % controls if an object is the support object
  vmin2 = C*vmin;        % controls if a support object is the boundary support object

	% Set up the variables for the optimization.
	n  = size(K,1);
	D  = (y*y').*K;
	f  = -ones(1,n);
	A  = y';
	b  = 0;
	lb = zeros(n,1);
	ub = repmat(C,n,1);
	p  = rand(n,1);
  D = (D+D')/2; % guarantee symmetry
	% Make the kernel matrix K positive definite.
	if Options.pd_check
    i = -30;
  	while (pd_check (D + (10.0^i) * eye(n)) == 0)
	  	i = i + 1;
	  end
	  if (i > -30),
	  	prwarning(2,'K is not positive definite. The kernel is regularized by adding 10.0^(%d)*I',i);
	  end
	  i = i+2;
	  D = D + (10.0^(i)) * eye(n);
  end

	% Minimization procedure initialization:
	% 'qp' minimizes:   0.5 x' K x + f' x
	% subject to:       Ax <= b
	%
	if (exist('qld') == 3)
		v = qld (D, f, -A, b, lb, ub, p, 1);
	elseif (exist('quadprog') == 2)
		prwarning(1,'QLD not found, the Matlab routine QUADPROG is used instead.')
		opt = optimset; opt.LargeScale='off'; opt.Display='off';
		v = quadprog(D, f, [], [], A, b, lb, ub,[],opt);
	else 
		prwarning(1,'QLD not found, the Matlab routine QP is used instead.')
		verbos = 0;
		negdef = 0;
		normalize = 1;
		v = qp(D, f, A, b, lb, ub, p, 1, verbos, negdef, normalize);
	end

  try
    % check if the optimizer returned anything
    if isempty(v)
      error('Optimization did not converge.');
    end

  	% Find all the support vectors.
	  J  = find(v > vmin1);
    Jp = J(y(J) ==  1);
    Jm = J(y(J) == -1);

    % Sanity check: there are support objects from both classes
    if isempty(J)
      error('There are no support objects.');
    elseif isempty(Jp)
      error('There are no support objects from the positive class.');
    elseif isempty(Jm)
      error('There are no support objects from the negative class.');
    end

    % compute nu parameter
	  nu = sum(v(J),1)/(C*n);	

   % Find the SV on the boundary
	  I = find((v > vmin1) & (v < C-vmin2));

    % Include class information into object weights
    v = y.*v;  

    % There are boundary support objects we can use them to find a bias term
    if ~isempty(I)
      b = mean(y(I)-K(I,J)*v(J));

    elseif Options.bias_in_admreg
      % There are no boundary support objects
      % We try to put the bias into the middle of admissible region

      % non SV
      J0    = (1:n)';
      J0(J) = [];
      J0p = J0(y(J0) ==  1);
      J0m = J0(y(J0) == -1);

      % Jp and Jm are all margin errors 
      lb = max(y([J0p;Jm]) - K([J0p;Jm],J)*v(J));
      ub = min(y([J0m;Jp]) - K([J0m;Jp],J)*v(J));

      if lb > ub
        error('The admissible region of the bias term is empty.');    
      end

      prwarning(2,['The bias term is undefined. The midpoint of its admissible region is used.']);    
      b = (lb+ub)/2;

    else
      error('The bias term is undefined.');
    end

	  v = [v(J); b];

  catch
    err.message = '##';
    lasterror(err);  % avoid problems with prwaitbar
    if Options.pf_on_failure
      prwarning(1,[lasterr ' Pseudo-Fisher is computed instead.']);
      n = size(K,1);
      %v = prpinv([K ones(n,1)])*y;
      v = prpinv([K ones(n,1); ones(1,n) 0])*[y; 0];
      J = [1:n]';
	    nu = nan;		  
    else
      rethrow(lasterror);
    end  
  end

return;
