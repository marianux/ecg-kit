%MOGC Mixture of Gaussian classifier
%
%   W = MOGC(A,N)
%   W = A*MOGC([],N);
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

% Copyright: R.P.W. Duin, r.p.w.duin@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: mogc.m,v 1.5 2008/05/21 11:49:59 duin Exp $

function w = mogc(a,n,r,s);

	prtrace(mfilename);

	if nargin < 4, s = 0; end
	if nargin < 3, r = 0; end
	if nargin < 2, n = 2; end
	if nargin < 1 | isempty(a)
		w = prmapping(mfilename,{n,r,s});
		w = setname(w,'MoG Classifier');
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
% 	d.mean = zeros(sum(n),k);
% 	d.cov = zeros(k,k,sum(n));
% 	d.prior = zeros(1,sum(n));
% 	d.nlab = zeros(1,sum(n));
% 	d.det = zeros(1,sum(n));

	d.mean = [];
	d.cov = [];
	d.prior = [];
	d.nlab = [];
	d.det = [];

	if(any(classsizes(a)<n))
		error('One or more class sizes too small for desired number of components')
	end

	% Estimate a MOG for each of the classes:
	w = [];
	for j=1:c
		b = seldat(a,j);
		%b = setlabtype(b,'soft');
        % para estimar bien cada componente, una regla de pulgar, 10 veces
        % datos la cantidad de variables.
        mb = getsize(b,1);
        ii = n(j) -  max(1, min(n(j), floor(mb/10/k)) );
        bContinue = true;
        while(bContinue)
            try
                v = gaussm(b,n(j)-ii,r,s);
                bContinue = false;
            catch PR_ME
                if(ii < (n(j)-1) && strcmpi(PR_ME.message, 'Not possible to find desired number of components'))
                    ii = ii + 1;
                else
                    rethrow(ME)
                end
            end
        end
        
        cant_comps = n(j)-ii;
        d.mean = [d.mean ; v.data.mean];
		d.cov = cat(3, d.cov, v.data.cov);
		d.prior  = [ d.prior; colvec(v.data.prior*p(j))];
		d.nlab   = [d.nlab; repmat(j, cant_comps, 1)];
		d.det    = [ d.det; colvec(v.data.det)];
        
	end
	
	w = prmapping('normal_map_new','trained',d,getlablist(a),k,c);
	%w = normal_map(d,getlablist(a),k,c);
	w = setname(w,'MoG Classifier');
	w = setcost(w,a);
    
    cFeaturesDomain = getfeatdom(a);
    w = setuser(w,cFeaturesDomain);
	
return;
