%COSTM  Cost mapping, classification using costs
%
%   Y = COSTM(X,C,LABLIST)
%   W = COSTM([],C,LABLIST)
%
% DESCRIPTION
% Maps the classifier output X (assumed to be posterior probability
% estimates) to the cost-outputs, defined by the cost-matrix C:
%
%   C(i,j) = cost of misclassifying an object from class i as class j.
%
% Default C is the cost matrix stored in the dataset X.
% The order of the classes is defined by LABLIST. When no lablist is
% given, the order as given by GETFEATLAB(X) is assumed.
% In order to apply this mapping, it is assumed that the dataset X
% represents posterior class probabilities (i.e. normalized by
% classc).
%
% EXAMPLE
% x = gendatb(100);
% w1 = x*ldc;                  % standard classifier
% C = [0 2; 1 0];              % new cost matrix
% w2 = w1*classc*costm([],C);  % classifier using this cost-matrix
%
% SEE ALSO
% MAPPINGS, CLASSC, TESTD, TESTCOST, SETCOST

% Copyright: D.M.J. Tax, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function w = costm(a,C,lablist)
		
% costm is implemented as a trained mapping to enable proper handling
% of constructs like W = fisherc(train)*costm([],C,Clab); test*W;

if nargin < 3, lablist = []; end
if nargin < 2, C = []; end
[s_in,s_out] = size(C);

if nargin < 1 | isempty(a)
	w = prmapping(mfilename,'combiner',{C,lablist});
	w = setname(w,'Mis-classification-cost');
	return
end

if ismapping(a)  % we do something like   w*costm()

  if isempty(lablist)
    lablist = getlabels(a);
  end
	if ~isempty(C) & ~isempty(lablist)
		% match the labels with that of the previous mapping
		[C,lablist] = matchcost(a.labels,C,lablist);
	end

                  % and now we have a trained mapping:
	w = prmapping(mfilename,'trained',{a,C,lablist},lablist,size(a,1),s_out);

elseif (isdataset(a) | isdouble(a)) & ismapping(C)
                  % we deal now with  d*costm()

	w = feval(mfilename, a*C.data{1}, C.data{2}, C.data{3});
	
elseif (isdataset(a) | isdouble(a)) & ~ismapping(C)
                  % the case that   d*costm([],C,lablist)

	if ~isempty(C) % store in dataset for standard error checking handling
						 % I doubt whether it is really the right thing to do
								 % as it may disturb a different labellingsystem in case
								 % we are dealing with a testset.
		a = setcost(prdataset(a),C,lablist);
  end
  
  if isdataset(a)
    C = a.cost; % retrieve from dataset
  else
    C = [];
  end

	if isempty(C) % still empty, no costs
		w = a;
	else
			% Get the data size
		d = size(a,2);
		if ( (size(C,1) ~= d) | (size(C,2) < d) )
			error('Cost matrix has wrong size');
		end
			% Map the dataset using the costs:
		w = a*C;
			% Now we have a problem: PRTools expects that the estimated class
			% has the highest output, and here we have calculated the cost
			% (which should be low!)
			% Thus we patch:
		w = -w;
		%w = setfeatlab(w,getfeatlab(a));
		if isempty(lablist)
			lablist = getfeatlab(a);
		end
		w = setfeatlab(w,lablist);
	end
else
	error('Input should be mapping or dataset');
end

return
