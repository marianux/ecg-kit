%TREE_MAP Map a dataset by binary decision tree
% 
%  F = TREE_MAP(A,W)
% 
% INPUT
%  A  Dataset
%  W  Decision tree mapping
%
% OUTPUT
%  F  Posterior probabilities
%
% DESCRIPTION
% Maps the dataset A by the binary decision tree classifier W on the 
% [0,1] interval for each of the classes W is trained on. The 
% posterior probabilities stored in F sum row-wise to one. W should 
% be trained by a classifier like treec. This routine is called 
% automatically to solve A*W if W is trained by treec.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% mappings, datasets, treec

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function [F,lab,N] = tree_map(T,W)

		% N yields for each node and each class in T the fraction of objects
% in T that passes that node.

%[tree,classlist,type,k,c,v,num] = prmapping(W);
% Unpack the classifier:
w = getdata(W);
tree = w{1};
num  = w{2};
[k,c] = size(W);

% Classification of the data vectors in T starts with node 
% num. In F the aposteriori probabilities for the classes
% of the final node are returned.
% N yields for each node and each class in T the fraction 
% of objects in T that passes that node.
% lab returns for each data vector the column for which F is
% maximum.
% tree(n,1) - feature number to be used in node n 
% tree(n,2) - threshold t to be used 
% tree(n,3) - node to be processed if value <= t 
% tree(n,4) - node to be processed if value > t 
% tree(n,5:4+c) - aposteriori probabilities for all classes in
% node n 
% If tree(n,3) == 0, stop, class in tree(n,1)

%[nlabt,lablistt,m,kt,ct,pt] = dataset_old(T);
% Setup the variables, also depending what outputs are requested by
% the user:
% [m,kt,ct] = getsize(T);
[m,kt] = size(T);
if kt ~= k, error('Wrong feature size'); end
[n,d] = size(tree);
lab = zeros(m,1);
if nargout==3
  ct = getsize(T,3);
	b = expandd(getnlab(T),ct);
	N = zeros(n,ct);
end
F = zeros(m,c);
node = num*ones(1,m);
for i = num:n
	S = find(node == i); 
	if nargout==3
		N(i,:) = sum(b(S,:));
	end
	if tree(i,3) > 0
		SL = S(find(+T(S,tree(i,1)) <= tree(i,2)));
		SR = S(find(+T(S,tree(i,1)) > tree(i,2)));
		node(SL) = tree(i,3)*ones(1,length(SL));
		node(SR) = tree(i,4)*ones(1,length(SR));
	elseif tree(i,3) == 0
		node(S) = inf * ones(1,length(S));
		lab(S) = tree(i,1) * ones(1,length(S));
		F(S,:) = tree(i*ones(length(S),1),5:4+c);
	else
      % right, what now?
	end
end
if nargout==3
	N=N./(ones(n,1)*(sum(b,1)./getprior(T)));
end

F = setdat(T,F,W);
%F = dataset_old(F,getlab(T),classlist,pt,lablistt);

return
