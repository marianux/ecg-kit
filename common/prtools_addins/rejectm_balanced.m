%REJECTM Rejection mapping
%
%    W = REJECTM(A,FRAC)
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
%  A = GENDATB;              % create trainingset
%  W = LDC(A);               % create supervised classifier
%  WR = REJECTM(A*W,0.05);   % reject 5% of the data
%  SCATTERD(A); PLOTC(W*WR); % show
%
% SEE ALSO
% REJECT, ROC, PLOTE

% Copyright: D.M.J. Tax, D.M.J.Tax@prtools.org
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function w = rejectm_balanced(a,thr,rejname)

if nargin<3
	rejname = 'reject';
end
if nargin<2
	thr = 0.05;
end
if nargin<1 | isempty(a)
	w = mapping(mfilename,{thr,rejname});
	w = setname(w,'rejection mapping');
	return
end

if ~ismapping(thr) %training
    
    %hago el balanceo primero del dataset
    orig_labels = getnlab(a);
    cs = classsizes(a);
    k_2_balance = round(max(cs)./ cs) - 1;
    
    for ii = 1:length(cs)
        aux_idx = findnlab(a,ii);
        for jj = 1:k_2_balance(ii)
            a = [a;a(aux_idx,:)];
        end
    end
    
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
	w = mapping(mfilename,'trained',W,newll,k,c+1);
	w = setname(w,'rejection mapping (balanced)');

else  % evaluation
	W = getdata(thr);
	m = size(a,1);

	% just add an extra reject-class, that will have the constant
	% threshold output:
	newout = [a repmat(W.thr,m,1)];
	w = setdat(a,newout,thr);
end

return
