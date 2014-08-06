%NEURC Automatic neural network classifier
% 
% 	W = NEURC(A,UNITS)
% 	W = A*NEURC([],UNITS)
% 	W = A*NEURC(UNITS)
%
% INPUT
%   A      Dataset
%   UNITS  Number of units
%          Default: 0.2 x size smallest class in A.
%
% OUTPUT
%   W      Trained feed-forward neural network mapping
%
% DESCRIPTION
% Automatically trained feed-forward neural network classifier with UNITS
% units in a single hidden layer. Training, by LMNC, is stopped when the 
% performance on an artificially generated tuning set of 1000 samples per 
% class (based on k-nearest neighbour interpolation) does not improve anymore.
%
% NEURC always tries three random initialisations, with fixed random seeds, 
% and returns the best result according to the tuning set. This is done in 
% order to obtain a reproducable result.
%
% If UNITS is NaN it is optimised by REGOPTC. This may take a long
% computing time and is often not significantly better than the default.
%
% Uses the Mathworks' neural network toolbox.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, LMNC, BPXNC, GENDATK, REGOPTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: neurc.m,v 1.9 2008/07/03 09:11:44 duin Exp $

function argout = neurc (varargin)

  checktoolbox('nnet');
  
	mapname = 'AutoNeuralNet';
  argin = shiftargin(varargin,'integer');
  argin = setdefaults(argin,[],[]);
  
  if mapping_task(argin,'definition')
    argout = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,units] = deal(argin{:});
    n_attempts = 3;			% Try three different random initialisations.
    [m,k] = size(a);
    if isempty(units)
      cs = classsizes(a);
      units = ceil(0.2*min(cs));
    end

  	if isnan(units) % optimize complexity parameter: number of neurons
			defs = {[]};
			parmin_max = [1,30];
			w = regoptc(a,mfilename,{units},defs,[1],parmin_max,testc([],'soft'),0);
			return
  	end
		
		islabtype(a,'crisp');
		isvaldfile(a,1,2); % at least 1 object per class, 2 classes
		a = testdatasize(a);
		a = setprior(a,getprior(a,0));

		% train a network.
		% Reproducability: always use same seeds. 
		randstate = randreset(1); opt_err = inf; opt_mapping = [];

		% Try a number of random initialisations.
		s = sprintf('%i neural network initializations: ',n_attempts);
		prwaitbar(n_attempts,s);
		for attempt = 1:n_attempts
			prwaitbar(n_attempts,attempt,[s int2str(attempt)]);
			prwarning(4,'training with initialisation %d of %d',attempt,n_attempts);
			t = gendatk(a,1000,2,1); 			% Create tuning set based on training set.
			w = lmnc(a,units,inf,[],t);		% Find LMNC mapping.
			e = t*w*testc;								% Calculate classification error.
			if (e < opt_err)							 
				% If this is the best of the three repetitions, store it.
				opt_mapping = w; opt_err = e;	
			end
		end
		prwaitbar(0);
    randreset(randstate); % return original state
		
		% Output is best network found.
		argout = setname(opt_mapping,mapname);

  else % Evaluation
    
    [a,w] = deal(argin{1:2});
		nodatafile(a);
		data = getdata(w); 

		if (length(data) > 1)
			
	    % "Old" neural network - network is second parameter: unpack.
  	  data = getdata(w); weights = data{1};
    	pars = data{2}; numlayers = length(pars);

	    output = a;                       % Output of first layer: dataset.
  	  for j = 1:numlayers-1
    	  % Number of inputs (n_in) and outputs (n_out) of neurons in layer J.
      	n_in = pars(j); n_out = pars(j+1);

	      % Calculate output of layer J+1. Note that WEIGHTS contains both
  	    % weights (multiplied by previous layer's OUTPUT) and biases
    	  % (multiplied by ONES).

      	this_weights = reshape(weights(1:(n_in+1)*n_out),n_in+1,n_out);
	      output = sigm([output,ones(m,1)]*this_weights);

  	    % Remove weights of this layer.
    	  weights(1:(n_in+1)*n_out) = [];
	    end
		else
			% "New" neural network: unpack and simulate using the toolbox.
			net = data{1};
			output = sim(net,+a')';
		end;

		% 2-class case, therefore 1 output: 2nd output is 1-1st output.
		if (size(output,2) == 1)
			output = [output (1-output)]; 
		end

		% Output is mapped dataset.
		argout = setdat(a,output,w);
	
	end

return
	
