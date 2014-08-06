%COSTM  Combiner mapping handling costs in classifications
%
%   Y = COSTM(X,C,LABLIST)
%   Y = X*COSTM([],C,LABLIST)
%   Y = X*COSTM(C,LABLIST)
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
% w2 = w1*classc*costm(C);  % classifier using this cost-matrix
% confmat(x*w1)
% confmat(x*w2)
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, CLASSC, TESTD, TESTCOST, SETCOST

% Copyright: D.M.J. Tax, R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function w = costm(varargin)
		
% costm is implemented as a combiner mapping to enable proper handling
% of constructs like W = fisherc(train)*costm([],C,Clab); test*W;

	argin = shiftargin(varargin,'numeric');
  argin = setdefaults(argin,[],[],[]);
  [a,C,lablist] = deal(argin{:});
  [s_in,s_out] = size(C);
  
  if mapping_task(argin,'definition') % definition
    
    w = prmapping(mfilename,'combiner',{C,lablist});
    w = setname(w,'Mis-classification-cost');
    
  elseif ismapping(a)  % we do something like   w*costm()

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
      % Don't forget to remove the costs!!
      w = setcost(w,[]);
      if isempty(lablist)
        lablist = getfeatlab(a);
      end
      w = setfeatlab(w,lablist);
    end
  else
    error('Input should be mapping or dataset');
  end

return
