%MDSC Manhatten Dissimilarity Space Classification 
%
%   W = MDSC(A,R,CLASSF)
%   W = A*MDSC([],R,CLASSF)
%   D = X*W
%
% INPUT
%   A       Dateset used for training
%   R       Dataset used for representation
%           or a fraction of A to be used for this.
%           Default: R = A.
%   CLASSF  Classifier used in dissimilarity space
%           Default LIBSVC([],[],100)
%   X       Test set.
%
% OUTPUT
%   W       Resulting, trained feature space classifier
%   D       Classification matrix
%
% DESCRIPTION
% This is a dissimilarity based classifier intended for a feature
% respresentation. The training set A is used to compute for every class
% its own eigenspace. All eigenvectors are used. A dissimilarity space is
% built by the Manhatten (L1, or Minkowsky-1 or city block) distances
% between training objects A or test objects X and the representation
% objects R after transformation (i.e. rotation) to the eigenspace of
% the class of the particular represention object. 
%
% Note that Euclidean distances are not affected by rotation, but Manhatten
% distances are.
%
% New objects in feature space can be classified by D = X*W or by
% D = PRMAP(X,W). Labels can be found by LAB = D*LABELD or LAB = LABELD(D).
%
% SEE ALSO
% DATASETS, MAPPINGS, FDSC

% Copyright: S.W. Kim, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = mdsc(a,r,classf)

	if nargin < 3 | isempty(classf), classf = libsvc([],[],100); end
	if nargin < 2, r = []; end
	if nargin < 1 | isempty(a)
		w = prmapping(mfilename,'untrained',{r,classf}); 
		w = setname(w,'ManhattenDisSpaceC');
		return
	end

	isvaldfile(a,2,2); % at least 2 objects per class, 2 classes

	if isempty(r)
		r = a;         % default representation set is the training set
	elseif is_scalar(r)
		[r,a] = gendat(a,r); % or we take a random fraction 
	end
	
	% Training set and representation set should belong to the same set of
	% classes. Let us check that.
	laba = getlablist(a);
	labr = getlablist(r);
	[nlab1,nlab2,lablist] = renumlab(laba,labr);
	c = size(lablist,1);
	if any([length(unique(nlab1)) length(unique(nlab2))] ~= c)
		error('Training set and representation set should contain the same classes')
	end
	
	% We are now ready to compute the classifier. The set of class dependent
	% rotations will be stored in a stacked set of mappings v.
	v = cell(1,c);
	for j=1:c  % compute the mapping for every class
		b = seldat(a,getclassi(a,lablist(j,:))); % training set of class j
		[e,d] = preig(covm(b)); % its eigenvalue decomposition
		u = affine(e);          % the rotation, apply it to the represention
		s = seldat(r,getclassi(r,lablist(j,:)))*u; % objects of the same class
		v{j} = u*proxm(s,'m',1); % store the rotation and the proximity mapping
		                         % to the rotated representation objects 
	end
	v = stacked(v);     % combine all mappings as a stacked mapping
	d = a*v;            % compute dissimilarity matrix for the training set
	n = disnorm(d);     % find a proper normalisation and ...
	w = a*(v*n*classf); % apply it to the training set, compute the classifier
	                    % and include the mapping for use by the test objects 
	w = setname(w,'ManhattenDisSpaceC');
	
return
	

	

