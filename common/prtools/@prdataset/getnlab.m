%GETNLAB Get numeric labels
%
%    [NLAB,LABLIST] = GETNLAB(A)
%
% The numeric labels of the dataset A are returned in NLAB.
% These are pointers to the list of labels LABLIST, so LABLIST(NLAB(i))
% is the label of object i. Note, however, that unlabeled objects or
% objects with a special label indication to be set by the user, may have
% numeric labels NLAB <= 0. Note also that for labels of type 'targets'
% the numeric labels NLAB are all set to 0. 
%
% SEE ALSO MULTI_LABELING
