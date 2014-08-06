% NAIVEBCC Naive Bayes Combining Classifier
% 
%   W = A*(WU*NAIVEBCC)
%   W = WT*NAIVEBCC(B*WT)
%   D = C*W
%     
% INPUT
%   A   Dataset used for training base classifiers as well as combiner
%   B   Dataset used for training combiner of trained base classifiers
%   C   Dataset used for testing (executing) the combiner
%   WU  Set of untrained base classifiers, see STACKED
%   WT  Set of trained base classifiers, see STACKED
% 
% OUTPUT
%   W   Trained Naive Bayes Combining Classifier
%   D   Dataset with prob. products (over base classifiers) per class
%
% DESCRIPTION
% During training the combiner computes the probabilities
% P(class | classifier outcomes) based on the crisp class assignements 
% made by the base classifiers for the training set. During execution the  
% product of these probabilities are computed, again following the crisp 
% class assignments of the base classifiers. These products are returned 
% as columns in D. Use CLASSC to normalise the outcomes. Use TESTD or 
% LABELD to inspect performances and assigned labels.
%
% NAIVEBCC differs from the classifier NAIVEBC by the fact that the 
% latter uses continuous inputs (no crisp labeling) and does not make a
% distinction between classifiers. Like almost any other classifier
% however, NAIVEBC may be used as a trained combiner as well.
% 
% REFERENCE
% 1. Kuncheva, LI. Combining pattern classifiers, 2004, pp.126-128.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, STACKED, NAIVEBC, CLASSC, TESTD, LABELD

% Copyright: Chunxia Zhang, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function V = naivebcc(A,W)

name = 'Naive Bayes combiner';

if nargin < 1 | isempty(A)
	% If there are no inputs, return an untrained mapping.
	V = prmapping(mfilename);
	V = setname(V,name);

elseif nargin == 1
	% call like V = NBC(A): train a NB classifier on the outcomes of the
	% base classifiers. So A has c x k columns: 
	% c outcomes (classes) x k classifiers

	islabtype(A,'crisp'); % allow crisp labels only
	isvaldfile(A,1,2); % at least one object per class, 2 classes

	A = testdatasize(A,'features'); % test whether they fit
	A = setprior(A,getprior(A,0));  % avoid many warnings
	[m,k,c] = getsize(A); % size of training set (m objects; k features; c classes)
	L = k/c; % compute the number of classifiers

	G = zeros(L*c,c); % confusion matrix for each base classifier
	for i = 1:L
		data = +A(:,(i-1)*c+1:i*c);
		[max_val,max_ind] = max(data,[],2);
		CM = zeros(c,c); % the (k,s)th entry indicates the number of elements 
		                 % whose true labels are omega_k and are assigned to omega_s
		for k = 1:c
			for s = 1:c
				CM(k,s) = sum(A.nlab == k & max_ind == s);
			end
		end
		G((i-1)*c+1:i*c,:) = CM'; % transpose the confusion matrix to facilitate 
		                          % the further computation
	end
	vector = classsizes(A);
	
	% Find all class probs for all classifiers P(class | classifiers)
	% P = (vector/m)*prod((G+1/c)./repmat(vector+1,size(G,1),1)); 
	P = (G+1/c)./repmat(vector+1,size(G,1),1);  % P(x|class)
	P = P.*repmat(getprior(A),size(P,1),1);     % P(x|class) P(class)
	q = sum(P,2);                               % P(x)
	P = P./repmat(q,1,c);                       % P(class|x)
	V = prmapping(mfilename,'trained',P,getlablist(A),L*c,c);
	V = setname(V,name);
	
elseif nargin == 2 & ismapping(W) % execution
	
	% call like V = A*W, in which W is a trained NBC
	
	isdataset(A);
	[m,k] = size(A);
	
	P = getdata(W);       % Prob matrix found by training
	c = size(P,2);        % Number of classes
	L = size(P,1)/c;      % number of classifiers

	% We assume that the k features of A (columns) are the results of
	% L classifiers producing c outcomes. We will construct a matrix 
	% M with for every classifier a 1 for its maximum class (most 
	% confident) class.
	
	M = zeros(m,k);
	for i = 1:L
		NM = zeros(m,c);
		S = (i-1)*c+1:i*c;
		[maxv,maxi] = max(+A(:,S),[],2);
		NM(sub2ind([m,c],[1:m]',maxi)) = ones(m,1);
		M(:,S) = NM;
	end
	
	% Now the confidences of the combiner will be estimated according
	% to the product over the classifiers of P(classifier|class)
	% estimated during training
	
	d = zeros(m,c);
	R = reshape(find(M'==1)-[1:c:m*c*L]',L,m)'+repmat([1:c:c*L],m,1);
	for j=1:c
		q = P(:,j)';
		d(:,j) = prod(q(R),2);
	end
	
	%d = d*normm; % normalise confidences, leave to user
	V = setdat(A,d,W); % store in dataset
	
end
return
