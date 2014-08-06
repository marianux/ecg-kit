% PERLC - Trainable linear perceptron classifier
% 
%   W = PERLC(A,MAXITER,ETA,W_INI,TYPE)
%   W = A*PERLC([],MAXITER,ETA,W_INI,TYPE)
%   W = A*PERLC(MAXITER,ETA,W_INI,TYPE)
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
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, NMC, FISHERC, BPXNC, LMNC, REGOPTC

% Copyright: D. de Ridder, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = perlc (varargin)
  
	mapname = 'Perceptron';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],100,0.1,[],'batch');
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a, maxiter, eta, w_ini, type] = deal(argin{:});
	
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
          w = setname(w,mapname);

    else   % multi-class classifier:

      w = mclassc(a,prmapping(mfilename,{maxiter,eta,w_ini}));
          w = setname(w,mapname);

    end	
    
  end
		
return
