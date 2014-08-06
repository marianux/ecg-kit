%GENDATH Generation of Highleyman classes
% 
%   A = GENDATH(N,LABTYPE)
%
% INPUT
%  N        Number of objects (optional; default: [50,50])
%  LABTYPE  Label type (optional; default: 'crisp')
%
% OUTPUT
%  A        Generated dataset
%
% DESCRIPTION
% Generation of a 2-dimensional 2-class dataset A of N objects
% according to Highleyman. 
%
% The two Highleyman classes are defined by 
% 1: Gauss([1 1],[1 0; 0 0.25]).
% 2: Gauss([2 0],[0.01 0; 0 4]).
% Class priors are P(1) = P(2) = 0.5 
%
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I, I = 1,2.
%
% LABTYPE defines the desired label type: 'crisp' or 'soft'. In the 
% latter case true posterior probabilities are set for the labels.
%
% Defaults: N = [50,50], LABTYPE = 'crisp'.
% 
% EXAMPLES
% PREX_PLOTC, PREX_CLEVAL
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, GAUSS, PRDATASETS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendath.m,v 1.5 2009/01/27 13:01:42 duin Exp $

function A = gendath(N,labtype)

		
	if nargin < 1, N = [50, 50]; end
	if nargin < 2, labtype = 'crisp'; end

	GA = [1 0; 0 0.25];
	GB = [0.01 0; 0 4];
	G = cat(3,GA,GB);
	p = [0.5 0.5];
	N = genclass(N,p);
	U = prdataset([1 1; 2 0],[1 2]','prior',p);
	U = setprior(U,p);
	A = gendatgauss(N,U,G);
	A = setname(A,'Highleyman Dataset');

	switch labtype
	 case 'crisp'
	  ;
	 case 'soft'
	  W = nbayesc(U,cat(3,GA,GB));
	  targets = A*W*classc;
	  A = setlabtype(A,'soft',targets);
	 otherwise
	  error(['Label type ' labtype ' not supported'])
	end

	return
