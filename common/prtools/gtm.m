%GTM  Trainable mapping fitting a Generative Topographic Mapping by EM
%
%   [W,L] = GTM (A,K,M,MAPTYPE,REG,EPS,MAXITER)
%
% INPUT
%   A       Dataset or double matrix
%   K       Vector containing number of nodes per dimension (default: [5 5], 2D map)
%   M       Vector containing number of basis functions per dimension (default: [10 10])   
%   MAPTYPE Map onto mean of posterior ('mean', default) or mode ('mode')
%   REG     Regularisation (default: 0)
%   EPS     Change in likelihood to stop training (default: 1e-5)
%   MAXITER Maximum number of iterations (default: inf)
%
% OUTPUT
%   W       GTM mapping
%   L       Likelihood
%
% DESCRIPTION
% Trains a Generative Topographic Mapping of any dimension, using the EM 
% algorithm.
%
% REFERENCES
% Bishop, C.M., Svensen, M. and Williams, C.K.I., "GTM: The Generative 
% Topographic Mapping", Neural Computation 10(1):215-234, 1998.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLOTGTM, SOM, PRPLOTSOM

% (c) Dick de Ridder, 2003
% Information & Communication Theory Group
% Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

function [w,L] = gtm (a, arg2, M, mapto, reg, epsilon, maxiter)

		step = 10;

    if (nargin < 2)
        prwarning(3,'number of nodes K not given, assuming [5 5] (2D map)');
        arg2 = [5 5];
    end;
    
	if (~ismapping(arg2))
		if (nargin < 7)
		    maxiter = inf;
		end;
		if (nargin < 6)
            prwarning(3,'stopping criteria not given, assuming EPSILON = 1e-10');
            epsilon = 1e-10;
		end;
		if (nargin < 5)
            prwarning(3,'regularisation REG not given, assuming 0');
		    reg = 0;
		end;
		if (nargin < 4)
            prwarning(3,'mapping type MAPTO not given, assuming mean');
		    mapto = 'mean';
		end;
        if (nargin < 3)
            prwarning(3,'number of basis functions M not given, assuming 10');
            M = 5;
        end;
    end;
    
	if (nargin == 0) | isempty(a)
		W = prmapping(mfilename); 
		W = setname(W,'GTM');
		return; 
	end

	% Prevent annoying messages: we will tell the user about any problems.
	
	warning off;
	
		a = testdatasize(a);
    t = +a'; [d,N] = size(t); [m,k] = size(a); 

    % If we're to apply a trained mapping, unpack its parameters.
    
    if (ismapping(arg2))
        w = arg2; data = getdata(w);
        K = data{1}; M = data{2};
        W = data{3}; sigma = data{4}; mapto = data{5};
    else
        K = arg2;
    end;
    
    % Create a KK-dimensional grid with dimensions K(1), K(2), ... of 
	% D-dimensional grid points and store it as a D x KK matrix X. Do the
	% same for the basis function centers PHI_MU, with grid dimensions M(1), ...

	K = reshape(K,1,prod(size(K)));		% Turn into vectors.
	M = reshape(M,1,prod(size(M)));

	% Check: either K or M should be a vector, or both should be vectors of
	% the same length.

	if (length(K) > 1) & (length(M) == 1)
		M = M * ones(1,length(K));
	elseif (length(M) > 1) & (length(K) == 1)
		K = K * ones(1,length(M));
	elseif (length(K) ~= length(M))
		error ('Length of vectors K and M should be equal, or K and/or M should be scalar.');
	end;
	
	D = length(K); KK = prod(K); MM = prod(M);

	if (D > d)
		error ('Length of vectors K and M should be <= the data dimensionality.');
    end;
    
    x         = makegrid(K,D);				% Grid points.
	phi_mu    = makegrid(M,D);				% Basis function centers.
	phi_sigma = 2/(mean(M)-1);				% Basis function widths.

	% Pre-calculate Phi.

	for j = 1:KK
  	 	for i = 1:MM
    	    Phi(i,j) = exp(-(x(:,j)-phi_mu(:,i))'*(x(:,j)-phi_mu(:,i))/(phi_sigma^2));
	    end;
    end;
	
    if (~ismapping(arg2))  % Train mapping.

		tries = 0; retry = 1;
		while ((retry == 1) & (tries <= 5))
	
            % Initialisation.
	
		    ptx = zeros(KK,N);
			R   = (1/KK)*ones(KK,N);	
			C = prcov(t'); [U,DD] = preig(C); [dummy,ind] = sort(-diag(DD));
			W = U(:,ind(1:D))*x*prpinv(Phi);
			if (size(C,1) > D)
				sigma = sqrt(DD(D+1,D+1));
			else
				sigma = 1/mean(M);
			end;
	
            done = 0; iter = 0; retry = 0; likelihood = 1.0e20;
	
            while ((~done) & (~retry))
	
      		    iter = iter + 1;
          		done = 1;
	
		  	    factor1 = (1/(2*pi*sigma^2))^(d/2);
		  	    factor2 = (-1/(2*sigma^2));
				WPhi    = W*Phi;
	
				for i = 1:KK
					ptx(i,:) = factor1 * exp(factor2*sum((WPhi(:,i)*ones(1,N)-t).^2));
				end;
	
				if (~retry)
	
					s = sum(ptx); s(find(s==0)) = realmin;		% Prevent divide-by-zero.
					R = ptx ./ (ones(size(ptx,1),1)*s);
					G = diag(sum(R'));
	
                  	% M-step #1 (eqn. 12)

          	        W = (prinv(Phi*G*Phi' + reg*eye(MM))*Phi*R*t')';
	
                  	% M-step #1 (eqn. 13)
          	
					s = 0; WPhi = W*Phi;
					for i = 1:KK
						s = s + sum(R(i,:).*sum((WPhi(:,i)*ones(1,N)-t).^2));
					end;
                  	sigma = sqrt((1/(N*d)) * s);
	
					% Re-calculate log-likelihood.
	
                   	prev_likelihood = likelihood; likelihood = sum(log(mean(ptx)+realmin));
	
					if (rem (iter,step) == 0)
                        prwarning (10, sprintf('[%3d] L: %2.2f (change: %2.2e)', ...
                            iter, likelihood, abs ((likelihood - prev_likelihood)/likelihood)));
					end;
	
					% Continue?
	
                  	done = (abs ((likelihood - prev_likelihood)/likelihood) < epsilon);
      			    done = (done | (iter > maxiter));

      			    if (~isfinite(likelihood))
      			        prwarning(3,'Problem is poorly conditioned, retrying');
      			        retry = 1; tries = tries + 1;
      			    end;
      			    
      	        end;
			end;
        end;
	
		L = likelihood;
	
		if (~retry)
    		w = prmapping(mfilename,'trained',{K,M,W,sigma,mapto},[],k,D);
		    w = setname(w,'GTM');
	    else
            prwarning(3,'Problem is too poorly conditioned, giving up');
            prwarning(3,'Consider lowering K and M or increasing REG');
	        w = [];
	    end;
	    
    else    % Apply mapping.
    
        factor1 = (1/(2*pi*sigma^2))^(d/2);
		factor2 = (-1/(2*sigma^2));
		WPhi    = W*Phi;
	
		for i = 1:KK
		    ptx(i,:) = factor1 * exp(factor2*sum((WPhi(:,i)*ones(1,N)-t).^2));
		end;
	
		s = sum(ptx); s(find(s==0)) = realmin;		% Prevent divide-by-zero.
		R = ptx ./ (ones(size(ptx,1),1)*s);

		switch(mapto)
		    case 'mean',
		        out = (x*R)';
		    case 'mode',
		        [dummy,ind] = max(R); out = x(:,ind)';
            otherwise,
		        error ('unknown mapping type: should be mean or mode')
		 end;

		 w = setdata(a,out,getlabels(w));
    end;
    
    warning on;
    
return

% GRID = MAKEGRID (K,D)
%
% Create a KK = prod(K)-dimensional grid with dimensions K(1), K(2), ... of 
% D-dimensional uniformly spaced grid points on [0,1]^prod(K), and store it 
% as a D x KK matrix X.

function grid = makegrid (K,D)

	KK = prod(K);

	for h = 1:D
		xx{h} = 0:(1/(K(h)-1)):1;	% Support point means
	end;

	% Do that voodoo / that you do / so well...

    if (D==1)
  	    xm = xx;
    else
  	    cmd = '[';
        for h = 1:D-1, cmd = sprintf ('%sxm{%d}, ', cmd, h); end; 
        cmd = sprintf ('%sxm{%d}] = ndgrid(', cmd, D);
        for h = 1:D-1, cmd = sprintf ('%sxx{%d}, ', cmd, h); end; 
        cmd = sprintf ('%sxx{%d});', cmd, D); eval(cmd);
    end;
    
	cmd = 'mm = zeros(D, ';
	for h = 1:D-1, cmd = sprintf ('%s%d, ', cmd, K(h)); end; 
	cmd = sprintf ('%s%d);', cmd, K(D)); eval (cmd);

	for h = 1:D
		cmd = sprintf ('mm(%d,', h);
		for g = 1:D-1, cmd = sprintf ('%s:,', cmd); end; 
        cmd = sprintf ('%s:) = xm{%d};', cmd, h); eval (cmd);
	end;

	grid = reshape(mm,D,KK);

return



