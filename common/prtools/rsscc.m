%RSSCC Trainable random subspace combining classifier
%
%    W = RSSCC(A,CLASSF,NFEAT,NCLASSF)
%    W = A*RSSCC([],CLASSF,NFEAT,NCLASSF)
%    W = A*RSSCC(CLASSF,NFEAT,NCLASSF)
%
% INPUT
%   A       Dataset
%   CLASSF  Untrained base classifier
%   NFEAT   Number of features for training CLASSF
%   NCLASSF Number of base classifiers
% 
% OUTPUT
%   W       Combined classifer
%
% DESCRIPTION
% This procedure computes a combined classifier consisting out of NCLASSF
% base classifiers, each trained by a random set of NFEAT features of A.
% W is just the set of base classifiers and still needs a combiner, e.g.
% use W*MAXC or W*VOTEC.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PARALLEL

function w = rsscc(varargin)
  
	mapname = 'rsscc';
  argin = shiftargin(varargin,'prmapping');
  argin = setdefaults(argin,[],nmc,[],[]);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'untrained',mapname);
    
  elseif mapping_task(argin,'training')			% Train a mapping.
  
    [a,classf,nfeat,nclassf] = deal(argin{:});
    isvaldset(a,1,1);
    [m,k] = size(a);
    if isempty(nfeat)
      nfeat = max(round(m/10),2); % use at least 2D feature spaces
    end
    if isempty(nclassf)
      nclassf = max(ceil(k/nfeat),10); % use at least 10 classifiers
    end
    if nfeat >= k  % allow for small feature sizes (k < nfeat)
      nfeat = k;   % use all features
      nclassf = 1; % compute a single classifier
    end
    nsets = ceil(nfeat*nclassf/k);
    featset = zeros(k,nsets);
    for j=1:nsets
      featset(:,j) = randperm(k)';
    end
    featset = featset(1:nfeat*nclassf);
    featset = reshape(featset,nclassf,nfeat);

    w = [];
    s = sprintf('Compute %i classifiers: ',nclassf);
    prwaitbar(nclassf,s);
    for j=1:nclassf
      prwaitbar(nclassf,j,[s num2str(j)]);
      w = [w; a(:,featset(j,:))*classf];
    end
    prwaitbar(0)
    w = prmapping(mfilename,'trained',{w,featset},getlablist(a),k,getsize(a,3));
    
  else % Evaluation
    
    [a,w] = deal(argin{1:2});
    wdata = getdata(w);
    w = wdata{1};
    featset = wdata{2}';
    w = a(:,featset(:))*w;
    
  end
  
return
	
	
		