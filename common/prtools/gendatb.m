%GENDATB Generation of banana shaped classes
% 
%    A = GENDATB(N,S)
% 
% INPUT
%       N    number of generated samples of vector with 
%            number of samples per class
%       S    variance of the normal distribution (opt, def: s=1)
%
% OUTPUT
%       A    generated dataset
%
% DESCRIPTION
% Generation of a 2-dimensional 2-class dataset A of N objects with a
% banana shaped distribution. The data is uniformly distributed along the
% bananas and is superimposed with a normal distribution with standard
% deviation S in all directions. Class priors are P(1) = P(2) = 0.5.
% Defaults: N = [50,50], S = 1.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, PRDATASETS

% Copyright: A. Hoekstra, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: gendatb.m,v 1.3 2008/03/20 07:53:58 duin Exp $

function a = gendatb(N,s)

		
	if nargin < 1 | isempty(N), N = [50,50]; end
	if nargin < 2, s = 1; end

   % Default size of the banana: 
	r = 5;
   % Default class prior probabilities:
	p = [0.5 0.5];
	N = genclass(N,p);

	domaina = 0.125*pi + rand(1,N(1))*1.25*pi;
	a   = [r*sin(domaina') r*cos(domaina')] + randn(N(1),2)*s;

	domainb = 0.375*pi - rand(1,N(2))*1.25*pi;
	a   = [a; [r*sin(domainb') r*cos(domainb')] + randn(N(2),2)*s + ...
         ones(N(2),1)*[-0.75*r -0.75*r]];
	lab = genlab(N);
	a = prdataset(a,lab,'name','Banana Set','lablist',genlab([1;1]),'prior',p);

return
