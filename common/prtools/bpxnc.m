%BPXNC Back-propagation trained feed-forward neural net classifier
% 
%   [W,HIST] = BPXNC (A,UNITS,ITER,W_INI,T,FID)
%
% INPUT
%   A      Dataset
%   UNITS  Array indicating number of units in each hidden layer (default: [5])
%   ITER   Number of iterations to train (default: inf)
%   W_INI  Weight initialisation network mapping (default: [], meaning 
%          initialisation by Matlab's neural network toolbox)
%   T      Tuning set (default: [], meaning use A)
%   FID    File descriptor to report progress to (default: 0, no report)
%
% OUTPUT
%   W      Trained feed-forward neural network mapping
%   HIST   Progress report (see below)
%
% DESCRIPTION 
% A feed-forward neural network classifier with length(N) hidden layers with 
% N(I) units in layer I is computed for the dataset A. Training is stopped 
% after ITER epochs (at least 50) or if the iteration number exceeds twice 
% that of the best classification result. This is measured by the labeled 
% tuning set T. If no tuning set is supplied A is used. W_INI is used, if 
% given, as network initialisation. Use [] if the standard Matlab 
% initialisation is desired. Progress is reported in file FID (default 0). 
%
% The entire training sequence is returned in HIST (number of epochs, 
% classification error on A, classification error on T, MSE on A, MSE on T).
% 
% Uses the Mathwork's Neural Network toolbox.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, LMNC, NEURC, RNNC, RBNC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: bpxnc.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function [w,hist] = bpxnc(varargin)

		[w,hist] = ffnc(mfilename,varargin{:});

	return
