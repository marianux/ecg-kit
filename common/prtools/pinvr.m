%PINVR PSEUDO-INVERSE REGRESSION (PCR)
% 
% 	[W,J,C] = PINVR(A,TYPE,PAR,EPS_TOL,MC,PD)
%
% INPUT
%   A	      Dataset
%   TYPE    Type of the kernel (optional; default: 'p')
%   PAR     Kernel parameter (optional; default: 1)
%   EPS_TOL Tolerance
%   MC      Do or do not data mean-centering (optional; default: 1 (to do))
%   PD      Do or do not the check of the positive definiteness (optional;
%           default: 1 (to do)) (not implemented) 
%
% OUTPUT
%   W       Mapping
%   J       Object identifiers of support objects		
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PROXM

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Based on SVC.M by D.M.J. Tax, D. de Ridder, R.P.W. Duin
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: pinvr.m,v 1.3 2010/02/08 15:29:48 duin Exp $

function [W, J] = pinvr(a,type,par,eps_tol,mc,pd)
    if nargin < 2 | ~isa(type,'prmapping')
    if nargin < 6
      pd = 1;  
    end

    if nargin < 5 |isempty(mc)
      mc = 1;  
    end

    if nargin < 4 
      eps_tol = [];  
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
      W = prmapping(mfilename,{type,par,eps_tol,mc,pd});
    	W = setname(W,['Pseudoinverse Regression']);
    	return;
    end

  	islabtype(a,'targets');
  	[m,k] = getsize(a);
    y = gettargets(a);	
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
      if ~isempty(eps_tol)
        tol = (m+1)*norm([K ones(m,1)])*eps_tol;  
        v = prpinv([K ones(m,1)],tol)*y;
      else
        v = prpinv([K ones(m,1)])*y;
      end
      
      J = [1:m]';

      % Store the results:
      v(end) = v(end)+uy; 
      W = prmapping(mfilename,'trained',{u,a(J,:),v,type,par},getlablist(a),k,1);
  		W = setname(W,['Pseudoinverse Regression']);
  		%W = setcost(W,a);
		J = getident(a,J);
  		%J = a.ident(J);
  		
  	else   
      error('multitarget regeression is not supported');
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
  	% classifier  w*x+b:
    d = [d ones(m,1)] * w{3};
    W = setdat(a,d,type);
	end
	
  return;
