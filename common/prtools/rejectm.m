%REJECTM Trainable mapping to compute rejection threshold
%
%    W = REJECTM(A,FRAC)
%    W = A*REJECTM([],FRAC)
%    W = A*REJECTM(FRAC)
%
% DESCRIPTION
% Train the threshold of a rejection mapping W such that a fraction FRAC
% of the training data A is rejected. Dataset A is usually the output of
% a classifier. The mapping REJECTM will add one extra reject class.
%
%    W = REJECTM(A,FRAC,REJNAME)
%
% If desired, the rejected objects will be labeled REJNAME. Default is
% REJNAME = 'reject'.
%
% EXAMPLES
%  a = gendatb;              % create trainingset
%  w = ldc(a);               % create supervised classifier
%  wr = rejectm(a*w,0.05);   % reject 5% of the data
%  scatterd(a); plotc(w*wr); % show
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% REJECT, PRROC, PLOTE

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = rejectm(varargin)
  
	mapname = 'Reject mapping';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0.05,'reject');
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,thr,rejname] = deal(argin{:});
    [n,k,c] = getsize(a);
    % add the new outlier class to the lablist
    newll = getlablist(a);
    if isa(newll,'double')
      %newll = [newll; max(newll)+1];
      if nargin>2 & isa(rejname,'char')
        warning('Labels are numeric, user supplied string class label.');
      end
      newll = [newll; rejname];
    else
      newll = char(newll,rejname);
    end

    % find the 'winning' class
     maxa = max(+a,[],2);
    % sort the posteriors for all of the classes:
    sa = sort(maxa);
    % find the thr-percentile and use that as a threshold:
    fracn = max(ceil(thr*n),1);
     thr = sa(fracn);

    % Store the threshold:
    W.thr = thr;
    W.c = c+1;
    w = prmapping(mfilename,'trained',W,newll,k,c+1);
    w = setname(w,'rejection mapping');

  else  % Evaluation
    
    [a,W] = deal(argin{1:2});
    thr = getdata(W,'thr');
    m = size(a,1);

    % just add an extra reject-class, that will have the constant
    % threshold output:
    newout = [a repmat(thr,m,1)];
    w = setdat(a,newout,W);
  end

return
