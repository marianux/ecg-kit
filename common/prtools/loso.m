%LOSO  Leave_One_Set_Out crossvalidation
%
%		[E,C,D] = LOSO(A,CLASSF,LABLISTNAME)
%		[E,C,D] = LOSO(A,CLASSF,SET_LABELS)
%		[E,C,D] = LOSO(A,CLASSF,SET_LABELS,SET_LABLIST)
%
% INPUT
%   A            Dataset
%   CLASSF       Untrained classifier
%   LABLISTNAME  Name of label list in case of multiple labeling
%   SET_LABELS   Set of labels for objects in A
%   SET_LABLIST  Order and selection of labels in SET_LABELS
%
% OUTPUT
%   E            Classification error
%   C            Array with numbers of erroneaously clasified objects
%                per label (vertically) and per class (horizontally)
%   D            Classification matrix of classified objects
%                (order may be different from A)
%
% DESCRIPTION
% In crossvalidation it may be desired that sets of corresponding objects
% (e.g. pixels from the same image) are all together in the training set or
% in the test set. This might be enabled by adding an additional labeling to
% the dataset A (see ADDLABELS) corresponding to the sets and running LOSO
% with the corresponding LABLISTNAME.
% Alternatively, the set labels may be supplied in the call. In SET_LABLIST
% a ranking of the used labels can be supplied that will be used in C.
% In case SET_LABLIST does not contain all set labels used in SET_LABELS
% LOSO will only test the set labels given in SET_LABLIST and thereby 
% perform areduced crosvalidation.
%
% The reported error E identical to E = sum(C)./classsizes(D)*getprior(A)';
% By confmat(D) a confusion matrix can be visualised.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, TESTC, CONFMAT

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [e,err,d] = loso(a,classf,set_lab,set_lablist)

if nargin < 4, testfun = testc; end

[curn,curlist] = curlablist(a);  % present lablist number
if nargin == 3 & size(set_lab,1) == 1
	if (isstr(set_lab) & strcmp(set_lab,curlist)) | set_lab == curn
		error('The desired LOO set should differ from the present lablist')
	end
	lablistname = set_lab;
	b = changelablist(a,lablistname);
elseif nargin == 3
	b = addlabels(a,set_lab); 
else
	b = addlablist(a,set_lablist);
	b = setnlab(b,renumlab(set_lab,set_lablist));
end
b = setlablist(b);  % throw out empty sets

nset = getsize(b,3);  
S = [1:nset 0];
c = getsize(a,3);
s = sprintf('Testing %i sets: ',nset);
prwaitbar(nset,s);
err = zeros(nset,c);
d = [];
N = 0;
for j=1:nset
	prwaitbar(nset,j,[s int2str(j)]);
	T = S;
	T(j) = [];
	x = changelablist(seldat(b,T),curlist);
	y = changelablist(seldat(b,j),curlist); 
	N = N+size(y,1);
	if ~isempty(y)
		dd = y*(x*classf);
		dd = prdataset(dd);
		[ee,err(j,:)] = testd(dd);
		d = [d;dd];
	end
end
prwaitbar(0);
e = sum(err)./classsizes(d);
e = e*getprior(a)';
	