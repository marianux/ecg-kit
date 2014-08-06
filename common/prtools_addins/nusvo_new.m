%NUSVO Support Vector Optimizer: NU algorithm
%
%   [V,J,NU,C] = NUSVO(K,NLAB,NU,OPTIONS)
%
% INPUT
%   K     Similarity matrix
%   NLAB  Label list consisting of -1/+1
%   NU    Regularization parameter (0 < NU < 1): expected fraction of SV (optional; default: 0.01)
%   OPTIONS
%    .PD_CHECK              force positive definiteness of the kernel by adding a small constant 
%                           to a kernel diagonal (default: 1)
%    .BIAS_IN_ADMREG        it may happen that bias of svc (b term) is not defined, then 
%                           if BIAS_IN_ADMREG == 1, b will be taken from its admissible
%                           region (if the region is bounded, the midpoint will be used);  
%                           if BIAS_IN_ADMREG == 0 the situation will be considered as 
%                           an optimization failure and treated accordingly (deafault: 1)
%    .ALLOW_UB_BIAS_ADMREG  it may happen that bias admissible region is unbounded 
%                           if ALLOW_UB_BIAS_ADMREG == 1, b will be heuristically taken 
%                           from its admissible region, otherwise(ALLOW_UB_BIAS_ADMREG == 0) 
%                           the situation will be considered as an optimization failure and 
%                           treated accordingly (deafault: 1)
%    .PF_ON_FAILURE         if the optimization is failed (optimizer did not converge, or there are 
%                           problems with finding of the bias term and PF_ON_FAILURE == 1, 
%                           then Pseudo Fisher classifier will be computed, 
%                           otherwise (PF_ON_FAILURE == 0) an error will be issued (default: 1)
%
%
% OUTPUT
%   V     Vector of weights for the support vectors
%   J     Index vector pointing to the support vectors
%   NU    NU parameter (useful when NU was automatically selected)
%   C     C regularization parameter of SVC algorithm, which gives the same classifier
%
%
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

% $Id: nusvo.m,v 1.3 2010/02/08 15:29:48 duin Exp $

function [v,J,nu,C] = nusvo_new(K,y,nu,Options,arg)
                 % old call: K,y,nu,pd,abort_on_error
	  prtrace(mfilename);

  % Old call conventions
  if nargin > 4 | (nargin == 4 & ~isstruct(Options) & isempty(Options))
    abort_on_error = [];
    if nargin == 5 
      abort_on_error = arg;  
    end

    pd = [];
    if nargin >= 4 
      pd = Options;
    end

    clear Options
    Options.pd_check             = pd;
    Options.bias_in_admreg       = 1;
    Options.allow_ub_bias_admreg = 1;
    Options.pf_on_failure        = ~abort_on_error;
  end    

  if nargin < 4
    Options = [];
  end


  DefOptions.pd_check             = 1;
  DefOptions.bias_in_admreg       = 1;
  DefOptions.allow_ub_bias_admreg = 1;
  DefOptions.pf_on_failure        = 1;

  Options = updstruct(DefOptions,Options,1);

  if nargin < 3
		prwarning(3,'The regularization parameter NU is not specified, assuming 0.01.');
    nu = 0.01;
	end

  nu_max = 1 - abs(nnz(y == 1) - nnz(y == -1))/length(y);
  nu_max = nu_max*(1-1e-6); % to be on a safe side
  
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
  D = (D+D')/2; % guarantee symmetry
  if Options.pd_check
  	% Make the kernel matrix K positive definite.
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
	% 'qp' minimizes:   0.5 x' D x + f' x
	% subject to:       Ax <= b
	%
% 	if (exist('qld') == 3)
% 		v = qld (D, f, -A, b, lb, ub, p, length(b));
% %         v =[];
% 
%         if (isempty(v) && exist('quadprog') == 2)
%             opts = optimset('Algorithm','active-set','Display','off');
%             v = quadprog(D, f, [], [], A, b, lb, ub, [], opts);
% 
%             if isempty(v)
%                 opts = optimset('Algorithm','trust-region-reflective','Display','off');
%                 v = quadprog(D, f, [], [], A, b, lb, ub, [], opts);
%             end
%         end
% 	elseif (exist('quadprog') == 2)
	if (exist('quadprog') == 2)
% 		prwarning(1,'QLD not found, the Matlab routine QUADPROG is used instead.')
% 		v = quadprog(D, f, [], [], A, b, lb, ub);
        
        opts = optimset('Algorithm','active-set','Display','off');
		v = quadprog(D, f, [], [], A, b, lb, ub, [], opts);
        
        if isempty(v)
            opts = optimset('Algorithm','trust-region-reflective','Display','off');
            v = quadprog(D, f, [], [], A, b, lb, ub, [], opts);
        end
        
	else 
		prwarning(1,'QLD not found, the Matlab routine QP is used instead.')
		verbos = 0;
		negdef = 0;
		normalize = 1;
		v = qp(D, f, A, b, lb, ub, p, length(b), verbos, negdef, normalize);
  end
  
  try
    % check if the optimizer returned anything
    if isempty(v)
      error('prtools:nusvo_new','Optimization did not converge.');
    end

	  % Find all the support vectors.
	  J = find(v > vmin);
    Jp = J(y(J) ==  1);
    Jm = J(y(J) == -1);

    % Sanity check: if there are support objects from both classes
    if isempty(J)
      error('There are no support objects.');
    elseif isempty(Jp)
      error('There are no support objects from the positive class.');
    elseif isempty(Jm)
      error('There are no support objects from the negative class.');
    end

    % Find the SV on the boundary
	  I  = J(v(J) < 1-vmin);
	  Ip = I(y(I) ==  1);
	  Im = I(y(I) == -1);

    % Include class information into object weights
    v = y.*v;  

    [rho, b] = FindRhoAndB(nu,y,K,v,J,Jp,Jm,Ip,Im,Options);
    
    v = [v(J); b]/rho;
    C = 1/rho;   

  catch
    if Options.pf_on_failure
      prwarning(1,[lasterr ' Pseudo-Fisher is computed instead.']);
      n = size(K,1);
      %v = prpinv([K ones(n,1)])*y;
      v = prpinv([K ones(n,1); ones(1,n) 0])*[y; 0];
      J = [1:n]';
	    C = nan;		  
    else
      rethrow(lasterror);
    end
  end  

return;


function [rho, b] = FindRhoAndB(nu,y,K,v,J,Jp,Jm,Ip,Im,Options)
  % There are boundary support objects from both classes, we can use them to find a bias term
  if ~isempty(Ip) & ~isempty(Im)
    % r1 = rho-b
    r1 =  mean(K(Ip,J)*v(J)); 
    
    % r2 = rho+b
    r2 = -mean(K(Im,J)*v(J)); 

  else 
    % Boundary support objects are absent at least in one of classes

    n = length(v);
    
    % non SV, we need them to resctric the addmissible region from above
    J0    = (1:n)';
    J0(J) = [];
    J0p   = J0(y(J0) ==  1);
    J0m   = J0(y(J0) == -1);

    % Only margin errors SV exist
    if isempty(Ip) & isempty(Im)

      % If only margin error SV's are present then it has to be true 
      % nSVp = nSVm, nu*n = nSVp + nSVm, it means that
      % rho and b should be in the region
      % lb1 <= rho-b <= ub1,  lb2 <= rho+b <= ub2
      % rho > 0 (rho = 0 corresponds to the zero margin solution, and cannot be 
      % interpreted as a standard SVC)
      
      if length(Jp) ~= length(Jm)
        error('Inconsistency: only margin errors SV exist and nSVp ~= nSVm.');
      elseif nu*n ~= length(J) 
        error('Inconsistency: only margin errors SV exist and nu*n ~= nSV.');
      elseif ~Options.bias_in_admreg
        error('The bias term is undefined.');
      end

      % suport objects are always exist, the region is alway bounded from below
      % if there are no nonSV, the region is unbounded from above

      % r1 = rho-b
      lb1 = max(K(Jp,J)*v(J));
      ub1 = min([K(J0p,J)*v(J); inf]);

      % r2 = rho+b
      lb2 = -min(K(Jm,J)*v(J));
      ub2 = -max([K(J0m,J)*v(J); -inf]);
      
      if lb1 > ub1 | lb2 > ub2 
        error('prtools:nusvo_new','Admissible region of the bias term is empty.');    
      end
    
      if isinf(ub1) | isinf(ub2)
        Msg = 'Admissible region of the bias term is unbounded';
        if Options.allow_ub_bias_admreg
          prwarning(1,Msg);
        else
          error('prtools:nusvo_new',Msg);
        end      
      end
      
      % admissible region is a a part of 
      % domain [lb1; ub1] x [lb2; ub2] which is above r2 = -r1 line
      %
      % we put (r1,r2) on a rectangle diagonal connecting upper right and bottom left 
      % corners half way from the point of intersection between this diagonal and line r2 = -r1; 
      %

      % first we make a rectangle bounded (put upper right corner above r2 = -r1 line)
      if isinf(ub1) 
        ub1 = max(-lb2,lb1)+1;
      end

      if isinf(ub2) 
        ub2 = max(-lb1,lb2)+1;
      end
      
      r1_intersection = (ub2*lb1 - ub1*lb2)/(ub2-lb2 + ub1-lb1);
      r2_intersection =  -r1_intersection; 
      
      % it may happen that intersecton does not exist: 
      %
      % 1) the whole domain is above line r2 = -r1, 
      % then we use bottom left corner insted of it 
      %
      % 2) the whole domain is below line r2 = - r1, it means that there is no solution 
      % with positive margin and the condition r1+r2 > 0 will be violated (we check it at the end)
      lb1 = max(lb1,r1_intersection);
      lb2 = max(lb2,r2_intersection);
      
      r1 = (lb1 + ub1)/2;
      r2 = (lb2 + ub2)/2;

    % Only margin errors SV are present in the negative class
    elseif ~isempty(Ip) & isempty(Im)
      % r1 = rho-b
      r1 =  mean(K(Ip,J)*v(J)); %  rho-b

      % r2 = rho+b
      lb2 = -min(K(Jm,J)*v(J));
      ub2 = -max([K(J0m,J)*v(J); -inf]);
      
      if lb2 > ub2 
        error('The admissible region of the bias term is empty.');    
      end
      
      % we need to minimize(2*mSVm - nu*n)*r2
      % s.t. lb2 <= r2 <= ub2, r2 > -r1
      coeff_sign = 2*length(Jm) - nu*n; 
      
      if coeff_sign > 0
        r2 = lb2; % put r2 to the lowest value, it has to be large than -r1 

      elseif coeff_sign < 0
        r2 = ub2;
        % If coeff_sign <0 it means, that nu*n > 2*mSVm 
        % Also, nu*n <= nu_max*n = n-abs(np-nm) = 2*min(np,nm)
        % mSVm < min(np,nm) <= nm
        % thus nm-mSVm = nonSVm+bSVm = nonSVm > 0 and ub2 < inf
        
      else % coeff_sing == 0
        if ~Options.bias_in_admreg
          error('The bias term is undefined.');
        end

        % correcting lb2 in such a way that rho >= 0, for (r1,lb2), 
        % it may happen that lb2 > ub2, let's check it later (r1+r2 > 0)
        lb2 = max(lb2,-r1);
        if isinf(ub2)
          Msg = 'Admissible region of the bias term is unbounded';
          if Options.allow_ub_bias_admreg
            prwarning(1,Msg);
            ub2 = lb2 + 1;
          else
            error(Msg);
          end      
        end
        r2 = (lb2+ub2)/2;
      end

    % Only margin errors SV are present in the positive class
    elseif isempty(Ip) & ~isempty(Im)
      % r1 = rho-b
      lb1 = max(K(Jp,J)*v(J));
      ub1 = min([K(J0p,J)*v(J); inf]);

      % r2 = rho+b
      r2 = -mean(K(Im,J)*v(J)); 
      
      if lb1 > ub1
        error('The admissible region of the bias term is empty.');    
      end
      
      % we need to minimize(2*mSVp - nu*n)*r1
      % s.t. lb1 <= r1 <= ub1, r2 > -r1
      coeff_sign = 2*length(Jp) - nu*n; 
      
      if coeff_sign > 0
        r1 = lb1; % put r1 to the lowest value, it has to be large than -r2 

      elseif coeff_sign < 0
        r1 = ub1; 
        % If coeff_sign <0 it means, that nu*n > 2*mSVp 
        % Also, nu*n <= nu_max*n = n-abs(np-nm) = 2*min(np,nm)
        % mSVp < min(np,nm) <= np
        % thus np-mSVp = nonSVp+bSVp = nonSVp > 0 and ub1 < inf
        
      else % coeff_sing == 0
        if ~Options.bias_in_admreg
          error('The bias term is undefined.');
        end

        % correcting lb1 in such a way that rho >= 0, for (lb1,r2), 
        % it may happen that lb1 > ub1, let's check it later (r1+r2 > 0)
        lb1 = max(lb1,-r2);
        if isinf(ub1)
          Msg = 'The admissible region of the bias term is unbounded.';
          if Options.allow_ub_bias_admreg
            prwarning(1,Msg);
            ub2 = lb1 + 1;
          else
            error(Msg);
          end      
        end
        r1 = (lb1+ub1)/2;
      end
    end

    if r1 + r2 <= 0
      error('The solution with the positive margin does not exist.');
    elseif isinf(r1) | isinf(r2)
      % This should never happen, because nu <=nu_max
      error('The contradiction is detected. Although nu <=nu_max, the rho (functional margin size) is infinite');
    end
  end

  rho = (r2+r1)/2;
  b   = (r2-r1)/2;

return
