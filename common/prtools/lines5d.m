%LINES5D  Generates three 5-dimensional lines
%
%	A = LINES5D(N);
%
% Generates a data set of N points, on 3 non-crossing, non-parallel lines
% in 5 dimensions. 
%
% If N is a vector of sizes, exactly N(I) objects are generated
% for class I, I = 1,2.Default: N = [50 50 50].
%
% See also DATASETS, PRDATASETS

% Copyright: E. Pekalska, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: lines5d.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function data = lines5d(N)
			if nargin< 1, N = [50 50 50]; end

	N = genclass(N,ones(1,3)/3);
	n1 = N(1);
	n2 = N(2);
	n3 = N(3);

  s1 = [0 0 0 1 0];
  s2 = [1 1 1 0 0];
  s3 = [0 1 0 1 0];
  s4 = [1 1 1 1 1];
  s5 = [0 1 1 0 1];
  s6 = [1 0 1 1 1];
  c1 = [0:1/(n1-1):1]';
  c2 = [0:1/(n2-1):1]';
  c3 = [0:1/(n3-1):1]';
  a  = c1*s1 + (1-c1)*s2;
  a  = [a; c2*s3 + (1-c2)*s4];
  a  = [a; c3*s5 + (1-c3)*s6];

	data = prdataset(a,genlab(N));
	data = setname(data,'5D Lines');

return
