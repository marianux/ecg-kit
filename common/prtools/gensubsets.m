%GENSUBSETS Generate sequence of nested training sets
%
%		[L,R] = GENSUBSETS(NLAB,S)
%		[L,R] = GENSUBSETS(A,S)
%
% INPUT
%   NLAB  Column vector of numeric labels of some dataset A.
%         NLAB = GETNLAB(A)
%   A     Dataset for which subsets are to be created
%   S     Array of growing subset sizes. 
%         S(K,J) should specify the size of training set K for class J with
%         numeric label J.
%
% OUTPUT
%   L     Cell array of length SIZE(S,1)+1 containing a series of growing
%					sets of indices or datasets. Datasets can be reconstructed from
%					indices by A(L{K},:). The last element of L refers to the
%					original dataset A
%   R     Cell array of length SIZE(S,1)+1 containing a series of shrinking
%					sets of indices or datasets. Datasets can be reconstructed from
%					indices by A(R{K},:). The last element of R is empty.
%
% DESCRIPTION
% Learning curves of classifier performances should be based on a
% consistent set of training sets, such that training set K1 is a subset of
% training set K2 if K1 < K2. This routine generates such a set on the
% basis of the numeric labels of the dataset A. L refers to the selected
% objects and R to the deselected ones.
%
% SEE ALSO
% DATASETS, GENDAT

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [L,R] = gensubsets(input,S)

if isdataset(input)
	nlab = getnlab(input);
else
	nlab = input;
end
[n,c] = size(S);  % number of subsets n, number of classes c
if max(nlab) ~= c
	error('Number of columns of size matrix not equal to number of classes')
end

m = length(nlab);
L = cell(1,n+1);
R = cell(1,n+1);
a = prdataset([1:m]',nlab); % make fake dataset with data equal indices
L{n+1} = [1:m]';
R{n+1} = [];
for j=n:-1:1
	[inset,outset] = gendat(a(L{j+1},:),S(j,:));
	L{j} = +inset;
	R{j} = [R{j+1};+outset];
end
if isdataset(input)
	for j=1:n+1
		L{j} = input(L{j},:);
		R{j} = input(R{j},:);
	end
end
		

