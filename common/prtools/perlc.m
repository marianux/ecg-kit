% PERLC - Train a linear perceptron classifier
% 
%   W = PERLC(A)
%   W = PERLC(A,MAXITER,ETA,W_INI,TYPE)
%
% INPUT
%   A        Training dataset
%   MAXITER  Maximum number of iterations (default 100)
%   ETA      Learning rate (default 0.1)
%   W_INI    Initial weights, as affine mapping, e.g W_INI = NMC(A)
%            (default: random initialisation)
%   TYPE     'batch': update by batch processing (default)
%            'seq'  : update sequentially
%
% OUTPUT
%   W        Linear perceptron classifier mapping
%
% DESCRIPTION
% Outputs a perceptron W trained on dataset A using learning rate ETA for a
% maximum of MAXITER iterations (or until convergence). 
%
% If ETA is NaN it is optimised by REGOPTC.
%
% SEE ALSO
% DATASETS, MAPPINGS, NMC, FISHERC, BPXNC, LMNC, REGOPTC

% Copyright: D. de Ridder, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: perlc.m,v 1.6 2008/07/03 09:11:44 duin Exp $

function w = perlc (a, maxiter, eta, w_ini, type)

	  
  if (nargin  < 5 | isempty(type))
    type = 'batch';
  end
	if (nargin < 4 | isempty(w_ini))
		prwarning(3,'No initial weights W_INI supplied, using random initialization');
		w_ini = []; 
	end
	if (nargin < 3 | isempty(eta))
		prwarning(3,'Learning rate ETA not specified, assuming 0.1');
		eta = 0.1;     
	end
	if (nargin < 2 | isempty(maxiter))
		prwarning(3,'Maximum number of iterations not specified, assuming 100');
		maxiter = 100; 
	end
	if (nargin < 1) | (isempty(a))
		w = prmapping(mfilename,{maxiter,eta,w_ini});
		w = setname(w,'Perceptron');
		return
	end
	
	if isnan(eta)    % optimize regularisation parameter
		defs = {100,0.1,[],'batch'};
		parmin_max = [0,0;1e-6,0.9;0,0;0,0];
		w = regoptc(a,mfilename,{maxiter, eta, w_ini, type},defs,[2],parmin_max,testc([],'soft'),1);
		return
	end
	
	% Unpack the dataset.
	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	[m,k,c] = getsize(a); 
	nlab = getnlab(a);

	% PERLC is basically a 2-class classifier. More classes are
	% handled by mclassc.
	
	if c == 2   % two-class classifier

		ws = scalem(a,'variance');
		a = a*ws;
		% Add a column of 1's for the bias term.
		Y = [+a ones(m,1)]; 

		% Initialise the WEIGHTS with a small random uniform distribution,
		% or with the specified affine mapping.
		if isempty(w_ini)
			weights = 0.02*(rand(k+1,c)-0.5);
		else
			isaffine(w_ini);
			weights = [w_ini.data.rot;w_ini.data.offset];
    end

    converged = 0; iter = 0;
		s = sprintf('perlc, %i iterations: ',maxiter);
		prwaitbar(maxiter,s,m*k>100000);
		while (~converged)

			% Find the maximum output for each sample.
			[maxw,ind] = max((Y*weights)');

      changed = 0;
      if (strcmp(type,'batch'))
        % Update for all incorrectly classified samples simultaneously.
      	changed = 0;
        for i = 1:m
          if (ind(i) ~= nlab(i))
        		weights(:,nlab(i)) = weights(:,nlab(i)) + eta*Y(i,:)';
        		weights(:,ind(i))  = weights(:,ind(i))  - eta*Y(i,:)';
        		changed = 1;
        	end;
        end;
        iter = iter+1;
      else
        % update for the worst classified object only
        J = find(ind' ~= nlab);
        if ~isempty(J)
          [dummy,imax] = min(maxw(J)); i = J(imax);
          weights(:,nlab(i)) = weights(:,nlab(i)) + eta*Y(i,:)';
     			weights(:,ind(i))  = weights(:,ind(i))  - eta*Y(i,:)';
       		iter = iter+1;
          changed = 1;
        end;
      end
			% Continue until things stay the same or until MAXITER iterations.
			converged = (~changed | iter >= maxiter);
			prwaitbar(maxiter,iter,[s int2str(iter)]);

    end
		prwaitbar(0);

		% Build the classifier
		w = ws*affine(weights(1:k,:),weights(k+1,:),a);
		w = cnormc(w,a);
		w = setlabels(w,getlablist(a));
        w = setname(w,'Perceptron');
	
	else   % multi-class classifier:
		
		w = mclassc(a,prmapping(mfilename,{maxiter,eta,w_ini}));
        w = setname(w,'Perceptron');
		
	end	
		
return
