%NU_SVRO Support Vector Optimizer
%
%   [V,J] = NU_SVRO(K,Y,C)
%
% INPUT
%   K     Similarity matrix
%   NLAB  Label list consisting of -1/+1
%   C     Scalar for weighting the errors (optional; default: 10)
%
% OUTPUT
%   V     Vector of weights for the support vectors
%   J     Index vector pointing to the support vectors
%
% DESCRIPTION
% A low level routine that optimizes the set of support vectors for a 2-class
% classification problem based on the similarity matrix K computed from the
% training set. SVO is called directly from SVC. The labels NLAB should indicate 
% the two classes by +1 and -1. Optimization is done by a quadratic programming. 
% If available, the QLD function is used, otherwise an appropriate Matlab routine.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% NU_SVR, SVO, SVC

% Revisions:
% DR1, 07-05-2003
%      Sign error in calculation of offset

% Copyright: D.M.J. Tax, D. de Ridder, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: nu_svro.m,v 1.2 2010/02/08 15:29:48 duin Exp $

function [v,J,epsilon_or_nu] = nu_svro(K,y,C,svr_type,nu_or_epsilon,pd,abort_on_error)

		if (nargin < 7) | isempty(abort_on_error)
    abort_on_error = 0;
  end

	if (nargin < 6) | isempty(pd)
    pd = 1;
  end

	if (nargin < 5) 
    nu_eps = [];
  end

	if (nargin < 4)  | isempty(svr_type)
    svr_type = 'nu';
  end

  switch svr_type
  case 'nu'
    if isempty(nu_or_epsilon)
      prwarning(3,'nu is not specified, assuming 0.25.');
      nu_or_epsilon = 0.25;
    end
    nu = nu_or_epsilon;
  case {'eps', 'epsilon'}
    svr_type = 'epsilon';
    if isempty(nu_or_epsilon)
      prwarning(3,'epsilon is not specified, assuming 1e-2.');
      nu_or_epsilon = 1e-2;
    end
    epsilon = nu_or_epsilon;
  end

	if (nargin < 3)
		prwarning(3,'C is not specified, assuming 1.');
    C = 1;
	end
	vmin = C*1e-9;		% Accuracy to determine when an object becomes the support object.

	% Set up the variables for the optimization.
	n  = size(K,1);
  D  = K;	
  switch svr_type
  case 'nu'
    f  = [-y', y'];  
    A  = [[ones(1,n), -ones(1,n)]; ones(1,2*n)];
    b  = [0; C*nu*n];
  case 'epsilon'
    f   = epsilon + [-y', y'];
  	A   = [ones(1,n), -ones(1,n)];
    b   = 0;
  end  
  
  lb = zeros(2*n,1);
	ub = repmat(C,2*n,1);
	p  = rand(2*n,1);

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
  D  = [[D, -D]; [-D, D]];	

	% Minimization procedure:
	% minimizes:  0.5 x' D x + f' x
	% subject to: Ax = b
	%
  if (exist('qld') == 3)
		v = qld (D, f, -A, b, lb, ub, p, length(b));
	elseif (exist('quadprog') == 2)
		prwarning(1,'QLD not found, the Matlab routine QUADPROG is used instead.')
		opt = optimset; opt.LargeScale='off'; opt.Display='off';
		v = quadprog(D, f, [], [], A, b, lb, ub,[],opt);
	else 
		prwarning(1,'QLD not found, the Matlab routine QP is used instead.')
		verbos = 0;
		negdef = 0;
		normalize = 1;
		v = qp(D, f, A, b, lb, ub, p, length(b), verbos, negdef, normalize);
	end
  % Find all the support vectors.
  
  if isempty(v) 
    ErrMsg = 'Quadratic Optimization failed.';
    [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
    return;
  end

  v  = v(1:n)-v((n+1):end);
  av = abs(v);
  J  = find(av > vmin);	    % sv's
  I  = J(av(J) < (C-vmin)); % on-tube sv's
  
  %plot(v,y-K(:,J)*v(J),'.')

  switch svr_type
  case 'nu'
    Ip = I(v(I) > 0);
    Im = I(v(I) < 0);  
    if isempty(Ip) | isempty(Im);
      prwarning(2,'epsilon and b are not unique: values from the admissible range will be taken');
      J0 = find(av <= vmin); % non-sv's
      if ~isempty(J0)
        y_minus_wx_0 = y(J0)-K(J0,J)*v(J); 
      else
        prwarning(2,'There are no non-sv''s: admissible (eps,b) range is partially unbounded');
      end
    end
    
    if ~isempty(Ip)
      b_plus_epsilon = mean(y(Ip)-K(Ip,J)*v(J),1); 
    else
      lp_bound = -inf;
      up_bound =  inf;

      if ~isempty(J0)
        lp_bound = max(y_minus_wx_0,[],1); 
      end
    
      Ipo = J((av(J) >= C-vmin) & v(J) > 0); % positive out of tube sv's
      if ~isempty(Ipo)
        up_bound = min(y(Ipo)-K(Ipo,J)*v(J),[],1); 
      end  
    
      if isinf(up_bound)
        Msg = 'Impossible situation: there are no positive sv''s.';    
        [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
        return;       
      elseif isinf(lp_bound)
        b_plus_epsilon = up_bound;
      else
        if lp_bound > up_bound
          ErrMsg = 'Impossible situation: admissible (eps,b) region is empty.';    
          [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
          return;       
        end
        b_plus_epsilon = 0.5*(lp_bound+up_bound);
      end  
    end
        
    if ~isempty(Im)
      b_minus_epsilon = mean(y(Im)-K(Im,J)*v(J),1); 
    else
      lm_bound = -inf;
      um_bound =  inf;
      
      Imo = J((av(J) >= C-vmin) & v(J) < 0); % positive out of tube sv's
      if ~isempty(Imo)
        lm_bound = max(y(Imo)-K(Imo,J)*v(J),[],1); 
      end  
       
      if ~isempty(J0)
        um_bound = min(y_minus_wx_0,[],1); 
      end
      
      if isinf(lm_bound)
        ErrMsg = 'Impossible situation: there are no negative sv''s.';    
        [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
        return;
      elseif isinf(um_bound)
        b_minus_epsilon = lm_bound;
      else
        if lm_bound > um_bound
          ErrMsg = 'Impossible situation: admissible (eps,b) range is empty.';    
          [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
          return;       
        end
        b_minus_epsilon = 0.5*(lm_bound+um_bound);
      end  
    end
    
    % one more paranoic check
    if exist('J0') == 1 & ~isempty(J0)
      ok = 1;
      if isempty(Ip) & ~isempty(Im) 
        ok = b_minus_epsilon <= min(y_minus_wx_0,[],1);
      elseif ~isempty(Ip) & isempty(Im)
        ok = b_plus_epsilon >= max(y_minus_wx_0,[],1);      
      end
      
      if ~ok
        ErrMsg = 'Impossible situation: incosistance in admissible (eps,b) region.';    
        [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
      end  
    end
    
    epsilon = 0.5*(b_plus_epsilon-b_minus_epsilon);
    b = 0.5*(b_plus_epsilon+b_minus_epsilon);

    epsilon_or_nu = epsilon;

  case 'epsilon'
    if ~isempty(I)      
      b = mean(y(I)-K(I,J)*v(J)-epsilon*sign(v(I)));
    else 
      prwarning(2,'b is not unique: value from the admissible range will be taken');
      
      lp_bound = -inf;
      up_bound =  inf;

      lm_bound = -inf;
      um_bound =  inf;

      J0 = find(av <= vmin); % non-sv's
      if ~isempty(J0)
        y_minus_wx_0 = y(J0)-K(J0,J)*v(J); 
        lp_bound = max(y_minus_wx_0,[],1)-epsilon; 
        um_bound = min(y_minus_wx_0,[],1)+epsilon; 
      else
        prwarning(2,'Thers are no non-sv''s');
      end
    
      Ipo = J((av(J) >= C-vmin) & v(J) > 0); % positive out of tube sv's
      if ~isempty(Ipo)
        up_bound = min(y(Ipo)-K(Ipo,J)*v(J),[],1)-epsilon; 
      end  
    
      Imo = J((av(J) >= C-vmin) & v(J) < 0); % negative out of tube sv's
      if ~isempty(Imo)
        lm_bound = max(y(Imo)-K(Imo,J)*v(J),[],1)+epsilon; 
      end  

      l_bound = max(lm_bound,lp_bound);
      u_bound = min(um_bound,up_bound);
      
      ErrMsg = '';
      if isinf(up_bound)
        ErrMsg = 'Impossible situation: there are no positive sv''s.';    
      elseif isinf(lm_bound)
        ErrMsg = 'Impossible situation: there are no negative sv''s.';    
      elseif l_bound > u_bound
        keyboard
        ErrMsg = 'Impossible situation: admissible b region is empty.';    
      end
      
      if ~isempty(ErrMsg)  
        [v,J,epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error);
        return;
      end         
      
      b = 0.5*(l_bound+u_bound); 
    end

    nu = sum(av(J))/(C*n);
    epsilon_or_nu = nu;
    
  end  
	
  v = [v(J); b];

return;

function [v,J, epsilon_or_nu] = ErrHandler(K,y,ErrMsg,abort_on_error)
  if abort_on_error
    error(ErrMsg);
  else
    prwarning(1,[ErrMsg 'Pseudoinverse Regression is computed instead.']);
    n = size(K,1);
    v = prpinv([K ones(n,1)])*y;
    J = [1:n]';
	  epsilon_or_nu = nan;		  
  end  
return

 
