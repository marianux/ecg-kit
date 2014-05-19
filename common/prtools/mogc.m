%MOGC Mixture of Gaussian classifier
%
%   W = MOGC(A,N)
%   W = A*MOGC([],N);
%   W = A*MOGC([],N,R,S);
%
%	INPUT
%    A   Dataset
%    N   Number of mixtures (optional; default 2)
%    R,S Regularization parameters, 0 <= R,S <= 1, see QDC
%	OUTPUT
%
% DESCRIPTION
% For each class j in A a density estimate is made by GAUSSM, using N(j)
% mixture components. Using the class prior probabilities they are combined 
% into a single classifier W. If N is a scalar, this number is applied to 
% each class. The relative size of the components is stored in W.DATA.PRIOR.
%
% EXAMPLES
% PREX_DENSITY
%
% SEE ALSO
% DATASETS, MAPPINGS, QDC, PLOTM, TESTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: mogc.m,v 1.6 2009/11/23 09:22:28 davidt Exp $

function w = mogc(a,n,r,s);

		if nargin < 4, s = 0; end
	if nargin < 3, r = 0; end
	if nargin < 2, n = 2; end
	if nargin < 1 | isempty(a)
		w = prmapping(mfilename,{n,r,s});
		w = setname(w,'MoG');
		return
	end
	
	islabtype(a,'crisp','soft');
	isvaldfile(a,n,2); % at least n objects per class, 2 classes
	
	% Initialize all the parameters:
	a = testdatasize(a);
	a = testdatasize(a,'features');
	[m,k,c] = getsize(a);
	p = getprior(a);
	a = setprior(a,p);
	if length(n) == 1
		n = repmat(n,1,c);
	end
    
	if length(n) ~= c
		error('Numbers of components does not match number of classes')
	end
	w = [];
	d.mean = zeros(sum(n),k);
	d.cov = zeros(k,k,sum(n));
	d.prior = zeros(1,sum(n));
	d.nlab = zeros(1,sum(n));
	d.det = zeros(1,sum(n));

	if(any(classsizes(a)<n))
		error('One or more class sizes too small for desired number of components')
	end

	% Estimate a MOG for each of the classes:
	w = [];
	n1 = 1;
	for j=1:c
		n2 = n1 + n(j) - 1;
		b = seldat(a,j);
		%b = setlabtype(b,'soft');
		v = gaussm(b,n(j),r,s);
		d.mean(n1:n2,:) = v.data.mean;
		d.cov(:,:,n1:n2)= v.data.cov;
		d.prior(n1:n2)  = v.data.prior*p(j);
		d.nlab(n1:n2)   = j;
		d.det(n1:n2)    = v.data.det;
		n1 = n2+1;
	end
	
	w = prmapping('normal_map','trained',d,getlablist(a),k,c);
	%w = normal_map(d,getlablist(a),k,c);
	w = setname(w,'MoG');
	w = setcost(w,a);
	
return;
