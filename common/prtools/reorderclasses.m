%REORDERCLASSES Reorder the lablist
%
%     X = REORDERCLASSES(X,LABLIST)
%     X = REORDERCLASSES(X,I)
%
% INPUT
%     X         (labeled) dataset
%     LABLIST   correctly ordered lablist
%     I         permutation vector
%
% OUTPUT
%     X         dataset with reordered classes
%
% DESCRIPTION
% Change the order of the classes in dataset X to the one defined in
% LABLIST. When only a single class label is given in LABLIST, this
% class is put in the first position, and the others are just kept in
% the order they already were.
%
% For example, when I have a dataset A with a labellist {'apple'
% 'banana' 'pear'} and I perform B=REORDERCLASS(A,3) (or identically
% B=REORDERCLASS(A,[3 1 2])), I will get a dataset B with a new
% labellist of {'pear' 'apple' 'banana'}.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%   RENUMLAB

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function x = reorderclasses(x,lablist)

% First deal with the situation that lablist is a cellarray:
if isa(lablist,'cell')
	if isa(lablist{1},'char')
		lablist = strvcat(lablist);
	else % it should be a double...
		if ~isa(lablist,'double')
			error('Only stringlabels or numeric labels are allowed.'); 
		end
		lablist = cell2mat(lablist);
	end
end
% ... and make sure numerical labels are column vectors
if isa(lablist,'double')
	lablist = lablist(:);
end

ll = getlablist(x);
c = size(ll,1);

% Second deal with the situation that we only supplied a single 'target
% class':
if size(lablist,1)==1
	if isa(lablist,'char') % the single label is a string:
		matchnr = strmatch(lablist,ll);
		if isempty(matchnr)
			error('Cannot find %s in the dataset.',lablist);
		end
	else % the single label is a number:
		% is this number a class label, or the index in the lablist??
		if (lablist<1) | (lablist>c) 
			% probably a class label
			matchnr = find(lablist==ll);
			if isempty(matchnr)
				error('Cannot find %d in the dataset.',lablist);
			end
		else
			matchnr = lablist;
		end
	end
	% make a complete permutation matrix: 
	lablist = (1:c)';
   lablist(matchnr) = [];
	lablist = [matchnr; lablist];
end

% test to see if lablist has at least enough elements:
if size(lablist,1)~=c
	error('New lablist or permutation vector does not match with number of classes.');
end

% Do the matching and reordering, depending if the class labels in X are
% string or numerical
if isa(ll,'char')  % X has string labels
	if isa(lablist,'char')  % we have to match 2 string labels
		I = zeros(c,1);
		for i=1:c
			matchnr = strmatch(lablist(i,:),ll);
			if isempty(matchnr)
				error('Label %s cannot be found in the dataset.',lablist(i,:));
			end
			I(i) = matchnr;
		end
	else % we have to apply a simple permutation on strings in X
		if ~isa(lablist,'double')
			error('Expecting a permutation label or a new lablist.');
		end
		I = lablist;
		if (min(I)<1) | (max(I)>c)
			error('Not a valid permutation vector is supplied.');
		end
		lablist = ll(I,:);
	end
else % dataset X has numeric labels

	if (min(lablist)<1) | (max(lablist)>c)
		% we are probably NOT dealing with a permutation vector
		I = zeros(c,1);
		for i=1:c
			matchnr = find(lablist(i)==ll);
			if isempty(matchnr)
				error('Label %d cannot be found in the dataset.',lablist(i));
			end
			I(i) = matchnr;
		end
	else % just a permutation on the numeric labels
		I = lablist;
		lablist = ll(I,:);
	end
end

% So finally perform the transformation:
x = setlablist(x, lablist);
[sortedI,I] = sort(I); %kind of inverse the I
% and be careful with unlabeled data
nlab = getnlab(x);
J = find(nlab>0);
nlab(J) = I(nlab(J));
x = setnlab(x, nlab);

return

