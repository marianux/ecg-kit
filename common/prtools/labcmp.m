%LABCMP Compare label sets
%
%	[JNE,JEQ] = LABCMP(LABELS1,LABELS2)
%
% INPUT
%   LABELS1 - list of labels (strings or numeric)
%   LABELS2 - list of labels (strings or numeric)
%
% OUTPUT
%   JNE     - Indices of non-matching labels
%   JEQ     - indices of matching labels
%
% DESCRIPTION
% The comparison of two label sets is particular useful to find
% erroneaously classified objects. For example, if W is a trained
% classifier and A a labeled testset, then the estimated labels are
% given by LAB_EST = A*W*LABELD, while the true labels can be found
% by LAB_TRUE = GETLABELS(A). The erroneously and correctly classfied
% objects can be found by [JNE,JEQ] = LABCMP(LAB_EST,LAB_TRUE).
%
% See also MAPPINGS, DATASETS, LABELD, TESTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [Jne,Jeq] = labcmp(lab1,lab2);

		
	[m1,n1] = size(lab1);
	[m2,n2] = size(lab2);
	if (m1 ~= m2)
		error('Label sets should have equal sizes')
	end
	
	[nn,J] = nlabcmp(lab1,lab2);
	Jne = find(~J);
	Jeq = find(J);

return
