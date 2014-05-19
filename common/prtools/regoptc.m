%REGOPTC Optimise regularisation and complexity parameters by crossvalidation
%
%		[W,PARS] = REGOPTC(A,CLASSF,PARS,DEFS,NPAR,PAR_MIN_MAX,TESTFUN,REALINT)
%
% INPUT
%   A       Dataset, training set
%   CLASSF  Untrained classifiers (mapping)
%   PARS    Cell array with parameters for CLASSF
%   DEFS    Defaults for PARS
%   NPAR    Index in PARS of parameters to be optimised
%   PAR_MIN_MAX Minimum and maximum values of the search interval for
%           the parameters to be optimised
%   TESTFUN Criterion function to be minimized, default TESTC
%   REALINT 0/1 vector, indicating for every parameter in PARS whether
%           it is real (1) or integer (0). Default: all real.
%
% OUTPUT
%   W       Best classifier, trained by A
%   PARS    Resulting parameter vector
%
% DESCRIPTIOM
% This routine is used inside classifiers and mappings to optimise a
% regularisation or complexity parameter. Using cross-validation the
% performance of the classifier is estimated using TESTFUN (e.g. TESTC).
% Matlab's FMINBND is used for the optimisation. Only the parameters in
% PARS that are set to NaN are optimised. For the other ones the given 
% values are used in the internal calls to CLASSF in REGOPTC. In case
% mulitple parameters are set to NaN they are optimised in the order
% supplied by NPAR. 
%
% The final parameters PARS can also be retrieved by GETOPT_PARS. This is
% useful if W is optimised inside training a classifier that does not
% return these parameters in the output.
%
% For examples of usage inside a classifier see LDC and SVC. Consequently
% LDC can be called as in the below example.
%
% EXAMPLE
% A = GENDATD([30 30],50);
% W = LDC(A,0,NaN); % set first reg par to 0 and optimise second.
% GETOPT_PARS       % retrieve optimal paameter set
%
% SEE ALSO
% DATASETS, MAPPINGS, CROSSVAL, TESTC, GETOPT_PARS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [w,varargout] = regoptc(a,classf,parms,defs,regparnum,regparmin_max,testfunc,realint)

	global REGOPT_NFOLDS REGOPT_REPS REGOPT_ITERMAX REGOPT_ITER REGOPT_OPTCRIT REGOPT_PARS

	if isempty(REGOPT_NFOLDS),  REGOPT_NFOLDS = 5; end
	if isempty(REGOPT_REPS),    REGOPT_REPS = 1; end
	if isempty(REGOPT_ITERMAX), REGOPT_ITERMAX = 20; end
	REGOPT_OPTCRIT = inf;
	REGOPT_PARS = [];

	isdataset(a);
	isuntrained(feval(classf,[],parms{:}));
	if nargin < 8, realint = ones(1,length(parms)); end
	if nargin < 7, testfunc = testc([],'crisp'); end
% 	if (length(parms) ~= length(defs)) | (length(parms) < max(regparnum)) | ...
% 			(length(parms) ~= size(regparmin_max,1)) | (length(regparnum) ~= length(realint)) | ...
% 			(size(regparmin_max,2) ~= 2)
% 		error('Some parameters have wrong size')
% 	end
	if (length(parms) ~= length(defs)) | (length(parms) < max(regparnum)) | ...
			(length(parms) ~= size(regparmin_max,1)) | (length(parms) ~= length(realint)) | ...
			(size(regparmin_max,2) ~= 2)
		error('Some parameters have wrong size')
	end

	J = [];
	K = zeros(1,length(parms));
	for j=1:length(parms)
		if ~isempty(parms{j}) & ~ismapping(parms{j}) & ~isstruct(parms{j}) & isnan(parms{j})
			J = [J j];
			K(j) = 1;   % parameters to be optimised
		end
	end
	parms(J) = defs(J);  % store defaults (needed in case of optimal parameters)
	matwarn = warning;
	warning off
	prwarn = prwarning;
	prwarning(0);
	prwaitbar(length(regparnum),'Parameter optimization');
	for j=1:length(regparnum)	
		prwaitbar(length(regparnum),j);
		n = regparnum(j);
		if K(n)
			regparmin = regparmin_max(n,1);
			regparmax = regparmin_max(n,2);
			if regparmin > 0 & regparmax > 0 & realint(n) % if interval positive and real
   		 	setlog = 1;                             % better to use logarithmic scaling
    		regparmin = log(regparmin);
    		regparmax = log(regparmax);
			else
    		setlog = 0;
			end
			REGOPT_ITER = 0;
			%if length(regparnum) == 1
				prprogress([],' par optim: %i steps, %i folds: \n', ...
						REGOPT_ITERMAX,REGOPT_NFOLDS);
					%else
			%	prprogress([],'%i-par optim: %i, %i steps, %i folds: \n', ...
			%			length(regparnum),n,REGOPT_ITERMAX,REGOPT_NFOLDS);
			%end
			prwaitbar(REGOPT_ITERMAX,'Parameter optimization');
			if realint(n) == 1
				regpar = fminbnd(@evalregcrit,regparmin,regparmax, ...
					optimset('Display','off','maxiter',REGOPT_ITERMAX), ...
					classf,a,parms,n,setlog,REGOPT_NFOLDS,REGOPT_REPS,testfunc,1);
			else
				regpar = nfminbnd(@evalregcrit,regparmin,regparmax,REGOPT_ITERMAX, ...
					classf,a,parms,n,setlog,REGOPT_NFOLDS,REGOPT_REPS,testfunc,0);
			end
			prwaitbar(0)
			if setlog
				parms{n} = exp(regpar);
			else
				parms{n} = regpar;
			end
		end
	end
	prwaitbar(0);
	varargout = cell(1,nargout-1);
	[w,varargout{:}] = feval(classf,a,parms{:});
	REGOPT_PARS = parms;
	warning(matwarn);
	prwarning(prwarn);

return

function regcrit = evalregcrit(regpar,classf,a,parms,regparnum, ...
	setlog,nfolds,reps,testfunc,realint);

	global REGOPT_ITER REGOPT_OPTCRIT REGOPT_ITERMAX
	REGOPT_ITER = REGOPT_ITER+1;
	
	prwaitbar(REGOPT_ITERMAX,REGOPT_ITER);
	
	if setlog
		parms{regparnum} = exp(regpar);
	else
		parms{regparnum} =regpar;
	end
	
	if realint
		prprogress([],' %i  %5.3f  %6.2e \n',REGOPT_ITER,REGOPT_OPTCRIT,parms{regparnum});
	else
		prprogress([],' %i  %5.3f  %i \n',REGOPT_ITER,REGOPT_OPTCRIT,parms{regparnum});
	end
		
	w = feval(classf,[],parms{:});
	rand('state',1); randn('state',1);
	regcrit = crossval(a,w,nfolds,reps,testfunc); % use soft error as criterion (more smooth)
	
	REGOPT_OPTCRIT = min(mean(regcrit),REGOPT_OPTCRIT);
	
return

