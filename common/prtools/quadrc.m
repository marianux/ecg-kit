
%QUADRC Trainable quadratic classifier
% 
%   W = QUADRC(A,R,S)
%   W = A*QUADRC([],R,S)
%   W = A*QUADRC(R,S)
% 
% INPUT
% 	A 		Dataset
%   R,S		0 <= R,S <= 1, regularization parameters (default: R = 0, S = 0)
% 
% OUTPUT
% 	W     Quadratic Discriminant Classifier mapping
%
% DESCRIPTION  
% Computation of the quadratic classifier between the classes of the dataset
% A based on different class covariances. R and S are regularization
% parameters used for finding the covariance matrix as
%
%		G = (1-R-S)*G + R*diag(diag(G)) + S*mean(diag(G))*eye(size(G,1))
%
% NOTE
% This routine differs from QDC; instead of using the densities, it is based
% on the class covariances. The multi-class problem is solved by multiple
% two-class quadratic discriminants. It is, thereby, a quadratic equivalent
% of FISHERC.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FISHERC, NMC, NMSC, LDC, UDC, QDC

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: quadrc.m,v 1.5 2010/02/08 15:29:48 duin Exp $

function w = quadrc(varargin)

  mapname = 'Quadr';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0,0);
  
  if mapping_task(argin,'definition')
    
   w = define_mapping(argin,'untrained',mapname);
    
	elseif mapping_task(argin,'training')			% Train a mapping.
 
    [a,arg2,s] = deal(argin{:});
    [m,k,c] = getsize(a);
		
		islabtype(a,'crisp');
		isvaldfile(a,2,2); % at least 2 objects per class, 2 classes
		a = testdatasize(a,'features'); 
		a = setprior(a,getprior(a));
		r = arg2;

		if (min(classsizes(a)) < 2)
			error('Classes should contain more than one vector.')
		end
			
		if (c == 2)

			% 2-class case: calculate quadratic discriminant parameters.

			%p = getprior(a); pa = p(1); pb = p(2);

			JA = findnlab(a,1); JB = findnlab(a,2);

			ma = mean(a(JA,:)); mb = mean(a(JB,:));

			GA = covm(a(JA,:)); GB = covm(a(JB,:));
			GA = (1-r-s) * GA + r * diag(diag(GA)) + ...
															s * mean(diag(GA))*eye(size(GA,1));
			GB = (1-r-s) * GB + r * diag(diag(GB)) + ...
															s*mean(diag(GB))*eye(size(GB,1));
			DGA = det(GA); DGB = det(GB);

			GA = prinv(GA);
			GB = prinv(GB);

			par1 = 2*ma*GA-2*mb*GB; par2 = GB - GA; 

			% If either covariance matrix is nearly singular, substitute FISHERC.
			% Otherwise construct the mapping.

			if (DGA <= 0) | (DGB <= 0)
				prwarning(1,'Covariance matrix nearly singular, regularization needed; using FISHERC instead')
				w = fisherc(a);
			else
				%par0 = (mb*GB*mb'-ma*GA*ma') + 2*log(pa/pb) + log(DGB) -log(DGA);
				par0 = (mb*GB*mb'-ma*GA*ma') + log(DGB) -log(DGA);
				w = prmapping(mfilename,'trained',{par0,par1',par2},getlablist(a),k,2);
				w = cnormc(w,a);
				w = setname(w,mapname);
				w = setcost(w,a);
			end

		else

			% For C > 2 classes, recursively call this function, using MCLASSC.

			pars = feval(mfilename,[],r,s);
			w = mclassc(a,pars);

		end

	else % Evaluation									

		% Second argument is a trained mapping: test. Note that we can only
		% have a 2-class case here. W's output will be [D, -D], as the distance
		% of a sample to a class is the negative distance to the other class.
    
    [a,v] = deal(argin{1:2});
    [m,k] = size(a);
		pars = getdata(v); 
		d = +sum((a*pars{3}).*a,2) + a*pars{2} + ones(m,1)*pars{1}; 
		w = setdat(a,[d, -d],v);

	end

return
	
function p = logdet(cova,covb)
%
% Numerically feasible computation of log(det(cova)/det(covb)))
% (seems not to help)

k = size(cova,1);

sa = mean(diag(cova));
sb = mean(diag(covb));

deta = det(cova./sa)
detb = det(covb./sb)

p = k*log(sa) - k*log(sb) + log(deta) - log(detb);
