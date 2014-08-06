%PRROC Receiver-Operator Curve
% 
%   E = PRROC(A,W,C,N)
%   E = PRROC(B,C,N)
%
% INPUT
%   A  Dataset
%   W  Trained classifier, or
%   B  Classification result, B = A*W*CLASSC
%   C  Index of desired class (default: C = 1)
%   N  Number of points on the Receiver-Operator Curve (default: 100)
%
% OUTPUT
%   E  Structure containing the error of the two classes
%
% DESCRIPTION
% Computes N points on the Receiver-Operator Curve (ROC)of the classifier W
% for class C in the labeled dataset B, which is typically the result of
% B = A*W; or for the dataset A labelled by applying the (cell array of)
% trained classifier(s) W.
%
% Note that a ROC is related to a specific class (class C) for which the
% errors are plotted horizontally. The total error on all other classes is
% plotted vertically. The class index C refers to its position in the label
% list of the dataset (A or B). It can be found by GETCLASSI.
%
% The curve is computed for N thresholds of the posteriori probabilities
% stored in B. The resulting error frequencies for the two classes are
% stored in the structure E. E.XVALUES contains the errors in the first
% class, E.ERROR contains the errors in the second class. In multi-class
% problems these are the mean values in a single class, respectively the
% mean values in all other classes. This may not be very useful, but not
% much more can be done as for multi-class cases the ROC is equivalent to a
% multi-dimensional surface.
%
% Use PLOTE(E) for plotting the result. In the plot the two types of error
% are annotated as 'Error I' (error of the first kind) and 'Error II' (error
% of the second kind). All error estimates are weighted according the class
% prior probabilities. Remove the priors in A or B (by setprior(A,[])) to
% produce a vanilla ROC.
%
% EXAMPLES
%	Train set A and test set T:
%	  B = T*NMC(A); E = PRROC(T,50); PLOTE(E); % Plots a single curve
%	  E = PRROC(T,A*{NMC,UDC,QDC});  PLOTE(E); % Plots 3 curves
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
%    John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 
%    2002.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PLOTE, REJECT, TESTC, GETCLASSI

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: roc.m,v 1.8 2009/07/09 14:17:38 davidt Exp $

function e = prroc(a,w,clas,n)

	% Depending on the call, CLAS may the third or second argument.
	% and N the third or the fourth.
	
	if nargin < 4, n = []; end
	if nargin < 3, clas= []; end
	if nargin < 2, w = []; end
	if nargin < 1, a= []; end
	
	if isempty(a)
		e = prmapping('roc','fixed',{w,clas,n});
		return;
	end
	
	name = [];
	num_w = 1;
	compute = 0;
	if isa(w,'double') | isempty(w)
		n = clas; 
		clas = w;
		b = a;
	elseif iscell(w)
		num_w = length(w); 
		compute = 1;
	elseif ismapping(w)
		istrained(w);
		b = a*w*classc;
		name = getname(w);
	else
			error('Classifier or cell array of classifiers expected.')
	end
	
	if isempty(n),
		n = 100;
		prwarning(4,'number of ROC points not given, assuming 100');
	end
	
	if isempty(clas)
		clas = 1;
		prwarning(4,'no class given, class 1 assumed');
	end

	p = getprior(a);
	datname = getname(a);
	lablist = getlablist(a,'string');
	if size(lablist,1)<2
		error('At least two classes should be present for computing the ROC curve.');
	end
	if clas > size(lablist,1)
		error('Wrong class number')
	end
	clasname= lablist(clas,:);
	%DXD: also check the class sizes:
	cs = classsizes(a);
	if cs(clas)==0
		error(['Class ',clas,' does not contain objects.']);
	end
	cs(clas) = 0; 
	if sum(cs)==0
		error(['One of the classes does not contain objects.']);
	end

	% Set up the ROC structure.

	e.error   = zeros(num_w,n+1);
	e.std     = zeros(num_w,n+1);
	e.xvalues = [0:n]/n;;
	e.n       = 1;
	e.xlabel  = 'Error I';
	e.ylabel  = 'Error II';
	e.plot    = 'plot';
	e.names   = name;
    if (~isempty(datname))
		e.title = ['ROC curve for the ' datname ', class ' clasname];
	end

	for j = 1:num_w												

	thr = [];				% Threshold range, will be set below.

		% If a cell array of classifiers was given, apply each one here.

		if (compute)
			v = w{j};
			if (~istrained(v))
				error('The supplied classifier mappings should be trained.')
			end
			b = a*v*classc;
			e.names = char(e.names,getname(v));
		else
			b = b*normm; % make sure we have a normalised classification matrix
		end

		[m,c] = size(b); nlab = getnlab(b); d = sort(b(:));

		% Attempt to compute a good threshold range on the first classified
		% dataset.

		if (isempty(thr))
			thr = [min(d)-eps d(ceil(([1:n-1]/(n-1))*m*c))' 1+eps]; 
        end
		e1 = []; e2 = [];

		% NLAB_OUT will be one where B is larger than THR.
		I = matchlablist(getlablist(b),getfeatlab(b)); % Correct possible changes class orders
		nlab_out = (repmat(+b(:,I(clas)),1,n+1) >= repmat(thr,m,1));

		% ERRS will be 1 where the numeric label is unequal to NLAB_OUT
		% (i.e., where errors occur).
		errs = (repmat((nlab==clas),1,n+1) ~= nlab_out);
			
		% Split the cases where ERRS = 1 into cases for CLAS (E1) and all 
		% other classes (E2).

		e1 = [mean(errs(find(nlab==clas),:),1)];
		e2 = [mean(errs(find(nlab~=clas),:),1)];
		
		% Fill results in the ROC structure.
		e.error(j,:)   = e2;
		e.xvalues(j,:) = e1;

	end

	% First name will be the empty string [], remove it.
	if (num_w > 1)
		e.names(1,:) = [];
	end

return
