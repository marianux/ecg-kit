%SETPRIOR Reset class prior probabilities of dataset
%
%   A = SETPRIOR(A,PROB,LABLIST)
%
% INPUT
%   A        Dataset
%   PROB     Prior probabilities to be set
%   LABLIST  Label list (optional)
%
% OUTPUT
%   A        Updated dataset
%
% DESCRIPTION
% Resets the class prior probabilities of the dataset A to PROB. PROB should
% be a vector of the length equal to the number of classes in A. In LABLIST,
% the corresponding class labels may be supplied. LABLIST may have only
% class names of the existing classes in A. Reset class names first by 
% SETLABLIST if necessary.
%
% If LABLIST is not given, the order defined by the existing LABLIST for A
% (determined by [NLAB,LABLIST] = renumlab(LABELS)) is used.
%
% PROB = 0 makes all C classes equally probable: 1/C.
% PROB = [] is interpreted as using the existing class frequencies in A as
% prior probabilities. Note that these prior probabilities change, if the
% number of elements in A is changed, or its labeling.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASET, GETPRIOR, ISEMPTY
