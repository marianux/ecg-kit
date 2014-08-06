%GENDATS Generation of a simple classification problem of 2 Gaussian classes
% 
%   A = GENDATS (N,K,D,LABTYPE)
% 
% INPUT
%   N       Dataset size, or 2-element array of class sizes (default: [50 50]).
%   K       Dimensionality of the dataset to be generated (default: 2).
%   D       Distance between class means in the first dimension (default: 1).
%   LABTYPE Label type to generate, 'crisp' or 'soft' (default: 'crisp').
%
% OUTPUT
%   A       Dataset.
%
% DESCRIPTION
% Generation of a K-dimensional 2-class dataset A of N objects. Both classes 
% are Gaussian distributed with identity matrix as covariance matrix. Their 
% means are on a distance D. Class priors are P(1) = P(2) = 0.5.
%
% If N is a vector of sizes, exactly N(I) objects are generated for class I, 
% I = 1,2.
%
% LABTYPE defines the desired label type: 'crisp' or 'soft'. In the latter 
% case true posterior probabilities are set for the labels.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: gendats.m,v 1.3 2009/01/27 13:01:42 duin Exp $

function A = gendats (N,k,d,labtype)

		if (nargin < 1), 
		prwarning(2,'class sizes not specified, assuming [50 50]');
		N = [50 50]; 
	end
	if (nargin < 2), 
		prwarning(2,'dimensionality not specified, assuming 2');
		k = 2; 
	end
	if (nargin < 3),
		prwarning(2,'class mean distance not specified, assuming 1');
		d = 2; 
	end
	if (nargin < 4), 
		prwarning(2,'label type not specified, assuming "crisp"');
		labtype = 'crisp'; 
	end

	% Set equal priors and generate random class sizes according to these.
	p = [0.5 0.5]; N = genclass(N,p);	

	% Unit covariance matrices, zero mean except for distance D in first dim.
	GA = eye(k); GB = eye(k);
	ma = zeros(1,k); mb = zeros(1,k); mb(1) = d;
	U = prdataset([ma;mb],[1 2]');
	U = setprior(U,p);

	% Create dataset.
	A = gendatgauss(N,U,cat(3,GA,GB),labtype);
	A = setname(A,'Simple Problem');

return
