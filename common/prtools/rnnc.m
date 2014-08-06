%RNNC Random Neural Net classifier
% 
% 	W = RNNC(A,N,S)
% 
% INPUT
%   A  Input dataset
%   N  Number of neurons in the hidden layer
%   S  Standard deviation of weights in an input layer (default: 1)
%
% OUTPUT
%   W  Trained Random Neural Net classifier
%
% DESCRIPTION
% W is a feed-forward neural net with one hidden layer of N sigmoid neurons.
% The input layer rescales the input features to unit variance; the hidden
% layer has normally distributed weights and biases with zero mean and
% standard deviation S. The output layer is trained by the dataset A.
% Default N is number of objects * 0.2, but not more than 100.
% 
% If N and/or S is NaN they are optimised by REGOPTC.
%
% Uses the Mathworks' Neural Network toolbox.
%
% REFERENCES
% 1. W.F. Schmidt, M.A. Kraaijveld, and R.P.W. Duin, Feed forward neural
% networks with random weights, Proc. ICPR11, Volume II, 1992, 1-4. 
% 2. G.B. Huang, Q.Y. Zhu, C.K. Siew, Extreme learning machine: theory and
% applications, Neurocomputing, 70 (1), 2006, 489-501
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, LMNC, BPXNC, NEURC, RBNC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: rnnc.m,v 1.4 2007/06/19 11:45:08 duin Exp $

function w = rnnc(a,n,s)

		
	checktoolbox('nnet');
	mapname = 'RandNeuralNet';

	if (nargin < 3)
		prwarning(3,'standard deviation of hidden layer weights not given, assuming 1');
		s = 1; 
	end
	if (nargin < 2)
		n = []; 
	end

	% No input arguments: return an untrained mapping.

	if (nargin < 1) | (isempty(a))
		w = prmapping('rnnc',{n,s});
		w = setname(w,mapname);
		return
	end

	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a);
	[m,k,c] = getsize(a);
	
	if isempty(n)
		n = min(ceil(m/5),100);
    prwarning(3,['no number of hidden units specified, assuming ' num2str(n)]);
	end
	
  if isnan(n) | isnan(s) % optimize complexity parameter: number of neurons, st. dev.
		defs = {m/5,1};
		parmin_max = [1,min(m,100);0.01,10];
		w = regoptc(a,mfilename,{n,s},defs,[1,2],parmin_max,testc([],'soft'),[0,1]);
		return
 	end

	% The hidden layer scales the input to unit variance, then applies a
	% random rotation and offset.
	
	w_hidden = scalem(a,'variance');
	w_hidden = w_hidden * cmapm(randn(n,k)*s,'rot');
	w_hidden = w_hidden * cmapm(randn(1,n)*s,'shift');

	% The output layer applies a FISHERC to the nonlinearly transformed output
	% of the hidden layer.

	w_output = w_hidden * sigm;
	w_output = fisherc(a*w_output);
	
	% Construct the network and insert the weights.

  warning('off','NNET:Obsolete');
	transfer_fn = { 'logsig','logsig','logsig' };
	net = newff(ones(k,1)*[0 1],[n size(w_output,2)],transfer_fn,'traingdx','learngdm','mse');
	net.IW{1,1} = w_hidden.data.rot'; net.b{1,1} = w_hidden.data.offset';
  if size(w_output,2) == 2
    net.LW{2,1} = w_output.data.rot'; net.b{2,1} = w_output.data.offset';
  else
    wrot = []; woff = [];
    for j=1:size(w_output,2)
      wrot = [wrot w_output{j}.data.rot];
      woff = [woff w_output{j}.data.offset];
    end
    net.LW{2,1} = wrot'; net.b{2,1} = woff';
  end

	w = prmapping('neurc','trained',{net},getlabels(w_output),k,c);
	w = setname(w,mapname);
	w = setcost(w,a);

return

	