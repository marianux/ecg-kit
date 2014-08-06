%FEATSELM Feature selection map
% 
% [W,R] = FEATSELM(A,CRIT,METHOD,K,T,PAR1,...)
% 
% INPUT
%   A      	Training dataset 
%   CRIT   	Name of criterion: 'in-in', 'maha-s', 'NN' or others 
%           (see FEATEVAL) or an untrained classifier V (default: 'NN')
%   METHOD  - 'forward' : selection by featself (default)
% 	        - 'float'   : selection by featselp
% 	        - 'backward': selection by featselb
% 	        - 'b&b'     : branch and bound selection by featselo
% 	        - 'ind'     : individual
% 	        - 'lr'      : plus-l-takeaway-r selection by featsellr
%           - 'sparse'  : use sparse untrained classifier CRIT
%   K      	Desired number of features (default: K = 0, return optimal set)
%   T      	Tuning set to be used in FEATEVAL (optional)
%   PAR1,.. Optional parameters:
% 	        - L,R       : for 'lr' (default: L = 1, R = 0)
%
% OUTPUT
%   W       Feature selection mapping
%   R       Matrix with step by step results     
%
% DESCRIPTION
% Computation of a mapping W selecting K features. This routines offers a
% central interface to all other feature selection methods. W can be used
% for selecting features in a dataset B using B*W.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, FEATEVAL, FEATSELO, FEATSELB, FEATSELI,
% FEATSELP, FEATSELF, FEATSELLR

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,res] = featselm(a,crit,arg3,ksel,t,par1,par2)

		
	if (nargin < 2 | isempty(crit))
		prwarning(2,'criterion not specified, assuming NN');
		crit = 'NN';        
	end
	if (nargin < 3 | isempty(arg3))
		prwarning(2,'method not specified, assuming forward');
		arg3 = 'forward'; 
	end
	if (nargin < 4)
		ksel = [];
	end
	if (nargin < 5)
		prwarning(3,'no tuning set supplied (risk of overfit)');
		t = [];             
	end
	if (nargin < 6), par1 = []; end;
	if (nargin < 7), par2 = []; end;

	% If no arguments are supplied, return an untrained mapping.

	if (nargin == 0) | (isempty(a))
		w = prmapping('featselm',{crit,arg3,ksel,t,par1,par2});
		w = setname(w,'Feature Selection');
		return
	end

	a = testdatasize(a);
	[m,k] = size(a);

	if (isstr(arg3))
		method = arg3;												% If the third argument is a string,
		switch (method)												%   it specifies the method to use.
		 case {'forward','featself'}											
		  [w,res] = featself(a,crit,ksel,t);
		 case {'float','featselp'}
		  [w,res] = featselp(a,crit,ksel,t);
		 case {'backward','featselb'}
		  [w,res] = featselb(a,crit,ksel,t);
		 case {'b&b','featselo'}
		  [w,res] = featselo(a,crit,ksel,t);
		 case {'ind','featseli'}
		  [w,res] = featseli(a,crit,ksel,t);
		 case {'lr','featsellr'}
		  [w,res] = featsellr(a,crit,ksel,par1,par2,t);
		 case {'sparse'}
		  v = a*crit;
		  if isaffine(v)
		  	v = getdata(v,'rot');
		  	w = featsel(size(a,2),find(v(:,1) == 0));	
			else
				v = getdata(v,'beta');
		  end
		  w = featsel(size(a,2),find(v(:,1) ~= 0));
		 otherwise
		  error('Unknown method specified.')
		end
	elseif (ismapping(arg3))								
		w = arg3;												% If the third argument is a mapping,
		isuntrained(w);									%  assert it is untrained and train
		[w,res] = feval(mfilename,a,crit,w.mapping_file,ksel,t,par1,par2);
	else
		error('Illegal method specified.')
	end

return
