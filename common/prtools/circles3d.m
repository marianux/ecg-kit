% CIRCLES3D  Create a data set containing 2 circles in 3 dimensions.
%
%   DATA = CIRCLES3D(N) 
%
%	Creates a data set containing N points in 3 dimensions.
%
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I, I = 1,2.Default: N = [50 50].
%
% See also DATASETS, PRDATASETS

% Copyright: E. Pekalska, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: circles3d.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function data = circles3d(N)
			if nargin< 1, N = [50 50]; end
	N = genclass(N,ones(1,2)/2);
	
	n2a = N(1);
	n2b = N(2);
	ha = 0:(2*pi/n2a):2*pi*(n2a/(n2a+1)); ha = ha';
	hb = 0:(2*pi/n2b):2*pi*(n2b/(n2b+1)); hb = hb';

	a = [ sin(ha) cos(ha) zeros(n2a,1) ];
	b = [ sin(hb) cos(hb) ones(n2b,1)  ];

	data = prdataset([a;b],genlab(N));
	data = setname(data,'3D Circles');

return
