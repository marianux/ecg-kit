%NAIVEBC Naive Bayes classifier
% 
% 	W = NAIVEBC(A,N)
% 	W = A*NAIVEBC([],N)
% 	W = A*NAIVEBC(N)
%
% 	W = NAIVEBC(A,DENSMAP)
% 	W = A*NAIVEBC([],DENSMAP)
% 	W = A*NAIVEBC(DENSMAP)
% 
% INPUT	
%   A       Training dataset
%   N       Scalar number of bins (default: 10)
%   DENSMAP Untrained mapping for density estimation
%
% OUTPUT
%   W       Naive Bayes classifier mapping
%
% DESCRIPTION
% The Naive Bayes Classifier estimates for every class and every feature
% separately. Total class densities are constructed by assuming
% independency and consequently multiplying the separate feature densities.
%
% The default version divides each axis into N bins, counts the number of
% training examples for each of the classes in each of the bins, and 
% classifies the object to the class that gives maximum posterior 
% probability. Missing values will be put into a separate bin.
%
% This routine assumes continuous data. It may be applied to discrete data
% in case all features have the same number of discrete values. For proper
% results the parameter N should be set to this number.
%
% If N is NaN it is optimised by REGOPTC.
%
% Alternatively an untrained mapping DENSMAP may be supplied that will be
% used to estimate the densities per class and per features separately.
% Examples are PARZENM and GAUSSM.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PARZENM, GAUSSM, UDC, QDC, PARZENC, PARZENDC, REGOPTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
  
% $Id: naivebc.m,v 1.5 2007/06/15 09:58:30 duin Exp $

function w = naivebc(varargin)

  mapname = 'NaiveBayes';
	argin = shiftargin(varargin,{'integer','prmapping'});
  argin = setdefaults(argin,[],10);
  
  if mapping_task(argin,'definition')
    
    w = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,arg2] = deal(argin{:});
    
    if ismapping(arg2) && isuntrained(arg2) % train given untrained mapping, e.g. parzen
		
      [m,k,c] = getsize(a);
      v = cell(c,k);
      for i=1:c
        b = seldat(a,i);
        for j=1:k
          v{i,j} = b(:,j)*arg2;
        end
      end
      pars.dens = v;
      pars.p0 = getprior(a);
      w = prmapping(mfilename,'trained',pars,getlablist(a),k,c);
      w = setname(w,'NaiveBayes');
      w = setcost(w,a);
    
    else % construct histograms
    
      N = arg2; % number of bins

      if isnan(N)             % optimize complexity parameter
        defs = {10};
        parmin_max = [2,50];
        w = regoptc(a,mfilename,{N},defs,[1],parmin_max,testc([],'soft'),0);

      else

        islabtype(a,'crisp');
        isvaldfile(a,1,2); % at least 2 object per class, 2 classes
        a = testdatasize(a);

        [m,k,c] = getsize(a); M = classsizes(a);

        % Train the mapping. First, find the scale and offset of the data 
        % and normalise (this is very non-robust, but ok...)

        offset_a = min(a); maxa = max(a); scale_a = maxa - offset_a;

        K = find(scale_a~=0);
  % 			if(any(scale_a==0))
  % 				prwarning (2,'one of the features has the same value for all data; scale change to realmin');
  % 				scale_a(scale_a==0) = realmin;		
  % 			end

        a = a - repmat(offset_a,m,1);
        a = a ./ repmat(scale_a,m,1);

        % P will contain the probability per bin per class, P0 the probability
        % per class. The highest and lowest bounds will not be used; the lowest
        % bound will be used to store the missing values.

        p  = zeros(N+1,k,c);

        % Count the number of objects for each of the classes.

        for i = 1:c

          Ic = findnlab(a,i);											% Extract one class.
          Ia = ceil(N*+(a(Ic,K)));								% Find out in which bin it falls.
          Ia(Ia<1) = 1; Ia(Ia>N) = N;							% Sanity check.

          for j=1:N
            p(j,K,i) = sum(Ia==j);								% Count for all bins.
          end

          p(N+1,K,i) = sum(~isnan(+a(Ic,K))); 	% The missing values.

          % Use Bayes estimators are used, like elsewhere in PRTools.

          p(:,K,i) = (p(:,K,i)+1) / (M(i)+N); 						% Probabilities.
          p(:,K,i) = p(:,K,i) ./ repmat(scale_a(K)/N,N+1,1); 	% Densities.

        end

        % Save all useful data.

        pars.p0 = getprior(a); pars.p = p; pars.N = N;
        pars.offset_a = offset_a; pars.scale_a = scale_a;
        pars.feats = K;

        w = prmapping(mfilename,'trained',pars,getlablist(a),k,c);
        w = setname(w,'Naive Bayes');
        w = setcost(w,a);

      end
		
    end

	else                      % Second argument is a mapping: testing.

    [a,w] = deal(argin{1:2});
		pars = getdata(w);		 	% Unpack the mapping.
		[m,k] = size(a); 
		
		if isfield(pars,'dens')
			
			v = pars.dens;
			[c,k] = size(v);
		
			out = zeros(size(a,1),c);
			for i=1:c
				for j=1:k
					out(:,i) = out(:,i) + +log(a(:,j)*v{i,j});
				end
			end
			
			out = exp(out);
			
		else
			
			c = length(pars.p0);       % Could also use size(w.labels,1)...
			K = pars.feats;            % relevant features

			% Shift and scale the test set.

			a(:,K) = a(:,K) - repmat(pars.offset_a(K),m,1);
			a(:,K) = a(:,K) ./ repmat(pars.scale_a(K),m,1);

			% Classify the test set. First find in which bins the objects fall.

			Ia = ceil(pars.N*+(a(:,K)));
			Ia(Ia<1) = 1; Ia(Ia>pars.N) = pars.N;			% Sanity check.

			% Find the class probability for each object for each feature
			out = zeros(m,length(K),c);
			for i=1:length(K)
				out(:,i,:) = pars.p(Ia(:,i),K(i),:);
			end

			% Multiply the per-feature probs.
			out = squeeze(prod(out,2));
		
			if m == 1
				out = out';
			end

		end
		
		% Weight with class priors
		out = out .* repmat(pars.p0,m,1);
		
		% Done!
		w = setdat(a,out,w);

	end

return



