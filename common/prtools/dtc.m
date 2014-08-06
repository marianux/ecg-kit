%DTC Verzakov Tree - Trainable Decision Tree Classifier
% 
%   W = DTC(A, CRIT, CHISQSTOPVAL, PRUNE, T, NFEAT)
% 
% INPUT
%   A       Training dataset. 
%           Object weights and class priors:
%             It is possible to assign individual weights to dataset
%             objects by using 'weights' identifier:
%
%             A = SETIDENT(A, WEIGHTS, 'weights')
%
%             If weights are not defined they assumed to be equal to 1 for
%             all objects. The actual weights used by the training routine
%             are computed by using supplied weights along with the class 
%             priors (if priors are not set then the apparent priors are used):
%
%             ACTUAL_WEIGHTS(I) = M*(WEIGHTS(I)/CW(C))*(PRIOR(C)/SUM(PRIOR))
%
%             where M is the total amount of (labelled) objects in A, 
%             CW is the vector weights sums for each class, and C is
%             the class of the object I. The sum of all actual weights is M.
%
%           Feature types:
%             Features are treated diffrently based on their domain information.
%             If the feature domain is empty or is the interval/collection of intervals  
%             then this feature is considered to be the continuous one and
%             branches are created by splitting at the threshold value.
%             Otherwise (feature domain specifies a set of values or set of
%             names), the feature is considered to be the nominal one and
%             branches corresponding to all present feature values are
%             created.
%
%           Unknown and 'non-applicable' values:
%           If the feature value is NaN then it is assumed that its value
%           is unknown. In such a situation the object with unknown value
%           is split into fractions and sent down all branches.
%           The concept of the 'non-applicable' feature value is different
%           from the concept of the unknown (missing) value. 
%           ...
%           The 'non-applicable' values for the continues features are
%           encoded as INF. If feature domain is the (set of) interval(s), 
%           then INF value has to be explicitly added to the domain
%           defintion. The 'non-applicable' value of nominal features
%           does not have predefined encoding. If it is necessary user
%           have to include such value into domain definition.
%
%   CRIT   Splitting citerion name.
%          'igr' Information Gain Ratio (default):
%          As defined by Quinlan. The penalty on the number of the distinct
%          values of the continues feature is used. If the gain is zero or
%          negative due to such penalization, the split is not performed.
%          This leads to smaller trees and may give non-zero training error.
%          This criterion does not use costs. (Costs are used only at the classification step).
%          
%          'gini' Gini impurity index. More precisely, the change in this
%          index. GINI index can be interpreted as a misclassification rate
%          for the stochastic prior based classifier, so costs are
%          naturally embedded. If the change in the (absolute) error less 
%          or equal to 0.1 (change in the cost less or equal to 0.1 of minimal 
%          absolute value of non-zero costs) the split is not performed.
%          This leads to smaller trees and may give non-zero training error.
%
%          'miscls' Classification error criterion.
%          To be used only for educational puposes because  
%          it gives rather inferior results. Costs are naturally embedded. 
% 
%   CHISQSTOPVAL Early stopping crtitical value for the chi-squared test 
%          on the difference between the original node class distribution and
%          branches class distributions. CHISQSTOPVAL is 0 by default. 
%          Which means that branches will be discarded only if they bring no
%          change in class distribution.
%
%   PRUNE  Pruning type name
%          'prunep' - pessimistic (top-down) pruning as defined by Quinlan. 
%          Pessimistic pruning be perfromed if cost matrix is defined.
%          'prunet' - test set (bottom-up) pruning using the (required)
%          dataset T.
%          These implementations of both pruning algorithms may be not
%          exactly correct if there are unknown values in datastes.
%
%   T      Test test for the test set pruning
%
% OUTPUT
%   W      Classifier mapping
%
% DESCRIPTION
%    If (full grown) branches of (sub)tree do not improve classification error 
%    (misclassification cost) they are immediatley discarded. 
%    This may happen because we use regularized posteriors. As a result
%    the algorithm is more stable, trees are smaller, but split on
%    the training set may be not perfect.
%
% REFERENCES
% [1] J.R. Quinlan, Simplifying Decision Trees, International Journal of 
%     Man-Machine Studies, 27(3), pp. 221-234, 1987.
% [2] J.R. Quinlan, Improved use of continuous attributes in C4.5. Journal
%     of AI Research, 4(96), pp. 77-90, 1996.
%
% see also DATASETS, MAPPINGS, TREEC

% Copyright: S. Verzakov, s.verzakov@gmail.com
% Based on prtools' treec, tree_map, and their subroutines
% by Guido te Brake and R.P.W Duin

function w = dtc(a, crit, chisqstopval, prune, t, nfeat)

		% When no input data is given, an empty tree is defined:
	if (nargin == 0) || isempty(a)
    
    if nargin < 2 || isempty(crit)
			crit = 'igr';
    end
    
    if nargin < 3 || isempty(chisqstopval)
      chisqstopval = 0;
    end
    
    if nargin < 4 || isempty(prune)
      prune = '';
    end
    
    if nargin < 5
      t = [];
    end
    
    if nargin < 6
      nfeat = [];
    end
    
    w = prmapping('dtc', {crit, chisqstopval, prune, t, nfeat});
    w = setname(w, ['DecTree' upper(crit)]);
	
  elseif (nargin == 2) && ismapping(crit)
    % Execution
    w = crit;    
    tree = +w;
    
    % B contains posteriors.
    % We do not need to convert posteriors to costs.
    % PRTools will do it automatically

    b = mapdt(tree, +a);
    b = setdat(a,b,w);
%    b = setfeatlab(b, getlabels(w));
    
    w = b;
  
  else
    % Training
	
    if nargin < 2 || isempty(crit)
      crit = 'igr'; 
    end

    if nargin < 3 || isempty(chisqstopval)
      chisqstopval = 0; 
    end
    
    if nargin < 4 || isempty(prune)
      prune = ''; 
    end

    if nargin < 5
      t = []; 
    end
    
    if nargin < 6
      nfeat = []; 
    end
    
    
    if ~any(strcmpi(crit, {'igr', 'gini', 'miscls'}))
      error('Unknown splitting criterion');
    end

    islabtype(a,'crisp');
    isvaldfile(a, 1, 2); % at least 1 object per class, 2 classes
    a = seldat(prdataset(a)); % get rid of all unlabelled objects
    
    m = size(a,1);
    % Get weights (if defined)
    weights = getident(a, 'weights');
    if isempty(weights)
      weights = ones([m 1]);
    else
     idx = (weights > 0);
     if nnz(idx) < m
       a = a(idx, :);  
       weights = weights(idx);
     end
     isvaldset(a, 1, 2); 
    end
      
    % First get some useful parameters:
    
    % Sizes
    [m, k, c] = getsize(a);
    cs = classsizes(a);

    % Determine if features are categorical
    featdom = getfeatdom(a);
    featdom = featdom(:)';
    if isempty(featdom)
      feattype = zeros(1, k);
    else
      feattype = cellfun(@(x) ischar(x) || (size(x,1)==1), featdom);
    end
    
    % Features to use
    usefeat = true([1, k]);
    if ~isempty(nfeat)
      nfeat = min(nfeat, k);
    end

    % Labelling
    nlab = getnlab(a);
    
    % Get priors
    prior = getprior(a);
    prior = prior/sum(prior);
    
    % Define weights:
    % the total sum of weights is equal to m
    % Compute class weights
    cw = zeros(size(cs));
    for i=1:c
      cw(i) = sum(weights(nlab == i));
    end
    % Rescale objects weights
    % by class weight factors derived from priors
    cwf = m*prior./cw;
    weights = weights.*(cwf(nlab)).';

    % Miscalssification cost
    cost = a.cost;

    % Now the training can really start:
    tree = makedt(+a, feattype, usefeat, nfeat, c, nlab, weights, cost, crit, chisqstopval);
    
    if ~isempty(prune)
      if strcmpi(prune, 'prunep')
        if ~isempty(cost)
          error('Pessimistic prunning based on misclassification costs is not implemented');
        end
        tree = prunep(tree);
      elseif strcmpi(prune, 'prunet')
        if isempty(t)
          error('Test set is not specified for the test set prunning');
        end
         tree = prunet(tree, t, cost);
      else
        error('Unkown prunning method');        
      end
      
      tree = cleandt(tree);
    end
    

    % Store the results:
    w = prmapping('dtc', 'trained', {tree}, getlablist(a), k, c);
    w = setname(w, ['DecTree' upper(crit)]);
    w = setcost(w, cost);
  end
  
  return

%MAKEDT General tree building algorithm
% 
%   TREE = MAKEDT(A, FEATTYPE, USEFEAT, NFEAT, C, NLAB, WEIGHTS, COST, CRIT, CHISQSTOPVAL)
% 
% INPUT
%   A      Data matrix
%
%   FEATTYPE Row defining the feature types (0 - continues, 1 - nominal)    
%
%   USEFEAT Row defining the features to be used for splitting
%
%   NFEAT  The maximum number of features for which splitting has to be
%          attempted. If the number of the available features is geater
%          than NFEAT the random subset of NFEAT features will be used for
%          splitting.
%
%  C       Total number of classes in the initial dataset
%
%  NLAB    Numeric class labels of objects in A
%
%  WEIGHTS Weights of objects in A
%  
%  COST    Misclassification costs
%
%  CRIT    Name of splitting criterion to be used
%
%  CHISQSTOPVAL Crtitical value for the chi-squared test
%
% OUTPUT
%   TREE   Structure containing arrays defining decision tree.
%          .NSMP The number of samples. NSMP(J) is the sum of weights of objects 
%          which reached the node J.
%
%          .POST Posterior probabilities. POST(J, :) is the class
%          distribution at the node J. Object weights are taken into
%          account. 0 and 1 probablities are avoided by perfroming Bayes
%          uniform priors regularization.
%          
%          .CIDX Class index. CIDX(J) is the class corresponding to the
%          node J if it considered as a leaf.
%
%          .ERRL Leaf error. ERRL(J) is the misclassification error (cost) 
%          on training samples which reached the node J if this node is
%          considered to be a leaf.
%
%          .ERRT Leaf error. ERRT(J) is the misclassification error (cost) 
%          on training samples which reached the node J if this node is
%          considered to be a root of the (sub)tree.
%
%          .SIZE Tree size. SIZE(J) is the (sub)tree size with root at the
%          node J.
%
%          .FIDX Feature index. FIDX(J) is the index of the feature on
%          which node J is split. FIDX(J) == 0 means that J is a leaf.
%        
%          .FVAL Feature value(s). For nominal features FVAL{J} is the set
%          of feature values observed at the node J (LENGTH(FVAL{J} == 0 
%          for the leaf, otherwise it is > 1).
%          For continues features, if LENGTH(FVAL{J}) == 1 then it contains
%          threshold THR for splitting into the left (<= THR) and the
%          right (> THR) branches. If FVAL{J} == [] (and J is not a leaf) 
%          then it means that split is perfromed between applicable and 
%          non-applicable values. 
%         
%          .NIDX Branch node indices. For nominal features NIDX{J,K} is
%          the index of branch node of node J with value FVAL{J,K} 
%          of feature FIDX(J). For continues features (if LENGTH(FVAL{J}) == 1) 
%          NIDX{J,1} is the index of the left branch node, NIDX{J,2} is the 
%          index of the right branch node. If FVAL{J} == [] 
%          (and J is not a leaf) then NIDX{J,1} is the index of branch node
%          containing objects with applicable values of feature FIDX(J) and
%          NIDX{J,2} is the index of branch node containing objects with 
%          non-applicable values of the same feature.
% 
% This is a low-level routine called by DTC.
% 
% See also IGR, GINI, MISCLS

function tree = makedt(a, feattype, usefeat, nfeat, c, nlab, weights, cost, crit, chisqstopval) 
	
  	
  % Construct the tree:    
  
  % Find (absolute) class frequencies
  C = zeros(1, c);
  for j=1:c
    C(j) = sum(weights(nlab == j)); 
  end
  
  NC = nnz(C);
  tree.nsmp = sum(C);
  
  % regularization by 'uniform' Bayesian priors;
  C0 = C;
  C = C + 1;
  p = C/sum(C);
  
  if isempty(cost)
    [maxpost, cidx] = max(p);
    errc = tree.nsmp*(1-maxpost);
    %errc = tree.nsmp - C0(cidx);
  else
    costp = p*cost;
    [mincost, cidx] = min(costp);
    errc = mincost;    
  end
  
  tree.post = p;
  tree.cidx = cidx;
  tree.errl = errc;
  tree.errt = errc;
  tree.size = 1;

  if NC ~= 1 % not a pure class dataset
    % now the tree is recursively constructed further:
		% use desired split criterion
    [fidx, fval, nb, bidx, chisqval] = findsplit(+a, feattype, usefeat, nfeat, c, nlab, weights, cost, crit);
    
    % When the stop criterion is not reached yet, we recursively split
		% further:
    if ~isempty(fidx) && (chisqval > chisqstopval)
      tree.fidx = fidx;
      tree.fval = {fval};

      if feattype(fidx) > 0 
        usefeat(fidx) = 0;
      end
      
      uidx = bidx == 0;
      nu = nnz(uidx);
      if nu > 0
        knsmp = tree.nsmp - sum(weights(uidx));
      end
      
      tree.nidx = {zeros(1, nb)};
      tree.errt(1) = 0; 
      for j=1:nb
        tree.nidx{1}(j) = tree.size(1) + 1;
        
        J = bidx == j;
        if nu == 0
          bweights = weights(J);
        else
          J = J | uidx;
          bweights = weights(J);
          buidx = uidx(J);
          bknsmp = sum(bweights(~buidx)); 
          bweights(buidx) = (bknsmp/knsmp)*bweights(buidx);
        end
        
        branch = makedt(+a(J, :), feattype, usefeat, nfeat, c, nlab(J), bweights, cost, crit, chisqstopval);
        
        branch.nidx = cellfun(@(x) x + tree.size(1)*(x>0), branch.nidx, 'UniformOutput', false);
        tree.errt(1) = tree.errt(1) + branch.errt(1);
        tree.size(1) = tree.size(1) + size(branch.nidx, 1);

        tree.nsmp = [tree.nsmp; branch.nsmp];
        tree.post = [tree.post; branch.post];
        tree.cidx = [tree.cidx; branch.cidx];
        tree.errl = [tree.errl; branch.errl];
        tree.errt = [tree.errt; branch.errt];
        tree.size = [tree.size; branch.size];          
        tree.fidx = [tree.fidx; branch.fidx];
        tree.fval = [tree.fval; branch.fval];
        tree.nidx = [tree.nidx; branch.nidx];
      end
    end
  end
  
  % no improvement in error (cost), rollback
  if (tree.size(1) > 1) && (tree.errt(1) >= tree.errl(1))
    tree.nsmp = tree.nsmp(1);
    tree.post = tree.post(1,:);    
    tree.cidx = tree.cidx(1);        
    tree.errl = tree.errl(1);    
    tree.errt = tree.errl(1); % sic!
    tree.size = 1;
  end  

  if tree.size(1) == 1
    % We reached the stop criterion or no further split is possible
    % so we make a leaf node:
    tree.fidx = 0;
    tree.fval = {[]};
    tree.nidx = {[]};
  end
	
	return

%MAPDT Tree mapping and node statistic calculation
% 
% 	[P, S] = MAPDT(TREE, A, NLAB, WEIGHTS)
% 

function [p, s] = mapdt(tree, a, nlab, weights)
    
  persistent dt st
  
  if isstruct(tree)
    dt = tree;
    
    if nargin < 4
      weights = [];
    end
  
    if nargin < 3
      nlab = [];
    end
  
    m = size(a, 1);
    [n, c] = size(dt.post);
    p = zeros([m, c]);
  
    if (nargout < 2) || isempty(nlab)
      st = [];    
      
      for i=1:m
        p(i, :) = mapdt(1, +a(i, :));
      end
      
    else
      st = zeros([n, c]);
      nlab (nlab > c) = 0;
      if isempty(weights)
        weights = ones(m, 1);
      end
      
      for i=1:m
        p(i, :) = mapdt(1, +a(i, :), nlab(i), weights(i));
      end
    end
    
    
    s = st;
    clear dt st
    
  else
    j = tree;

    while j > 0  
      if ~isempty(st) && (nlab > 0)
        st(j, nlab) = st(j, nlab) + weights;
      end
      
      k = 0;
      fidx = dt.fidx(j);

      if fidx ~= 0
        nidx = dt.nidx{j};
        aval = a(fidx);

        if ~isnan(aval)
          fval = dt.fval{j};

          if isempty(fval)
            k = nidx(2 - (aval ~= inf));            
          elseif length(fval) == 1
            if aval ~= inf
              k = nidx(2 - (aval <= fval));            
            elseif length(nidx) == 3
              k = nidx(3);
            end;
          else
            k = nidx(aval == fval);
            if isempty(k)
              k = 0;
            end
          end
          
        else
          p = zeros([1, size(dt.post, 2)]);
          if isempty(st) || (nlab <= 0)
            for b=1:length(nidx)
              k = nidx(b);
              f = (dt.nsmp(k)/dt.nsmp(j));
              p = p + f*mapdt(nidx(b), a);
            end
          else
            for b=1:length(nidx)
              k = nidx(b);
              f = (dt.nsmp(k)/dt.nsmp(j));
              p = p + f*mapdt(nidx(b), a, nlab, f*weights);
            end
          end
          k = -1;
        end
      end

      if k == 0
        p = dt.post(j, :);
      end

      j = k;
    end
  end
    
  return

  
%FINDSPLIT General routine for finding the best split
% 
% 	[FIDX, FVAL, NB, BIDX, CHISQVAL, CRITVAL] = FINDSPLI(A, FEATTYPE, USEFEAT, NFEAT, C, NLAB, WEIGHTS, COST, CRIT)
% 

function [fidx, fval, nb, bidx, chisqval, critval] = findsplit(a, feattype, usefeat, nfeat, c, nlab, weights, cost, crit)
	
    
  selfeatidx = find(usefeat);
  nf = length(selfeatidx);
  if ~isempty(nfeat) && (nfeat < nf)
    permidx = randperm(length(selfeatidx));
    selfeatidx = selfeatidx(permidx(1:nfeat));
    nf = nfeat;
  end

  fval = cell([1, nf]);
  chisqval = zeros([1, nf]);
  critval = nan([1, nf]);
  
  % repeat for all selected features
  for f=1:nf
    fidx = selfeatidx(f);
    af = a(:, fidx);
    
    kidx = ~isnan(af); % known values index
    if nnz(kidx) == 0
      continue
    end
    
    MU = sum(weights(~kidx)) + realmin;

    af = af(kidx);
    wk = weights(kidx);
    nlabk = nlab(kidx);
    
    switch feattype(fidx)
      case 0 % continous/ordered feature
        naidx = af == inf;
        NA = repmat(realmin, [1 c]);
        MNA = c*realmin;
        nlabna = nlabk(naidx);
        nna = length(nlabna);
        if nna > 0
          for j = 1:c
            NA(1,j) = sum(wk(nlabna == j)) + realmin;
          end
          MNA = sum(NA);  
        end
        
        apidx = find(~naidx);
        af = af(apidx);
        wap = wk(apidx);
        nlabap = nlabk(apidx);
        [af, sortidx] = sort(af);
        wap = wap(sortidx);
        nlabap = nlabap(sortidx);
        
        if length(af) == 1
          uv = af;
          ns = 0;
          sli = [];
        else
          labchngcount = cumsum([1; double(diff(nlabap) ~= 0)]);
          uniquevallowidx = find([true; ((diff(af)./(0.5*abs(af(1:end-1) + af(2:end)) + realmin)) > 1e-8)]);
          uniquevalhighidx = [(uniquevallowidx(2:end) - 1); length(af)];
          % unique values 
          uv = af(uniquevalhighidx);

          % split low indices in unique values
          sli = find(labchngcount(uniquevalhighidx(2:end)) - labchngcount(uniquevallowidx(1:end-1)) > 0);
          % split low indices in af
          splitlowidx = uniquevalhighidx(sli);
          ns = length(splitlowidx);
        end
        
        if (ns == 0) && (nna == 0)
          continue
        end
        
        % applicable, left, and right branch class counts
        AP = zeros(1, c);
        L = zeros(ns, c); 
        R = zeros(ns, c);        
        
        for j = 1:c
          J = find(nlabap == j);
          mj = length(J);
          AP(j) = sum(wap(J));
          if (ns > 0) && (mj > 0)
            L(:, j) = (repmat(splitlowidx, [1, mj]) >= repmat(J.', [ns, 1]))*wap(J) + realmin;
            R(:, j) = AP(j) - L(:, j) + realmin;
          end
        end
        AP = AP + 2*realmin;
        
        % total count of applicable
        MAP = sum(AP);
        
        % known
        K = AP + NA;
        MK = MNA + MAP;

        % object counts for branches
        ML = sum(L, 2);
        MR = sum(R, 2); 
        
        [cv, i] = feval(crit, 0, MU, MK, K, MNA, NA, MAP, AP, ML, L, MR, R, cost, weights, uv, sli);
        critval(f) = cv;

        if ~isnan(cv) && (cv > -inf) 
          if (ns > 0) && ~isempty(i)
            t = splitlowidx(i);
            fval{f} = 0.5*(af(t) + af(t+1));

            if nna == 0 
              APL = AP*(ML(i)/MAP);
              APR = AP*(MR(i)/MAP);
              chisqval(f) = sum(((L(i, :) - APL).^2)./APL + ((R(i, :) - APR).^2)./APR);
            else
              KNA = K*(MNA/MK);
              KL = K*(ML(i)/MK);
              KR = K*(MR(i)/MK);
              chisqval(f) = sum(((NA - KNA).^2)./KNA + ((L(i, :) - KL).^2)./KL + ((R(i, :) - KR).^2)./KR);
            end
          else
            fval{f} = [];
            KNA = K*(MNA/MK);
            KAP = K*(MAP/MK);
            chisqval(f) = sum(((NA - KNA).^2)./KNA + ((AP(i, :) - KAP).^2)./KAP);
          end
        end
        
      case 1 % nominal feature
        vf = unique(af);
        v = length(vf);

        B = zeros(v, c);
        for j=1:c
          J = find(j == nlab);
          mj = length(J);
          if mj > 0
            B(:, j) = (repmat(vf, [1, mj]) ==  repmat(af(J).', [v 1]))*weights(J);
          end
        end
        
        B = B + realmin;
        
        % class counts for known
        K = sum(B, 1);

        % total known count
        MK = sum(K);

        % object counts for branches
        MB = sum(B, 2);
        
        cv = feval(crit, 1, MU, MK, K, MB, B, cost, weights);        
        critval(f) = cv;
        
        if ~isnan(cv) && (cv > -inf)
          KB = MB.*K/MK;
          chisqval(f) = sum(sum(((B - KB).^2)./KB));
          fval{f} = bf;
        end
    end
  end
    
  % best criterion over all features
  testfeatidx = find(~isnan(critval) & (critval > -inf));
  if isempty(testfeatidx)
    fidx = [];
    fval = [];
    nb = [];
    bidx = [];
    chisqval = [];
    critval = [];    

  else
    [critval, fidx] = max(critval(testfeatidx));
    fidx = testfeatidx(fidx);
    fval = fval{fidx};
    chisqval = chisqval(fidx);    
    fidx = selfeatidx(fidx);
    
    af = a(:, fidx);
    m = size(af, 1);
    bidx = zeros(m, 1);

    switch feattype(fidx)
      case 0
        apidx = af < inf;

        if ~isempty(fval)
          lidx = af <= fval;
          ridx = apidx & ~lidx;
          bidx(lidx) = 1;
          bidx(ridx) = 2;
          nb = 2;
        else  
          bidx(apidx) = 1;
          nb = 1;
        end
        
        naidx = ~isnan(af) & ~apidx;
        if nnz(naidx) > 0
          nb = nb + 1;
          bidx(naidx) = nb;            
        end
      
      case 1
        nb = length(fval);
        for i=1:nb
          bidx(af == fval(i)) = i;
        end
    end
  end
  
  return

  
%IGR The information gain ratio
% 
% 	[CRITVAL, IDX] = IGR(FEATTYPE, MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, COST, WEIGHTS, UV, SLI)
%   [CRITVAL, IDX] = IGR(FEATTYPE,MU, MK, K, MB, B, COST, WEIGHTS

function [critval, idx] = igr(feattype, varargin) 
	
  	
  switch feattype
    case 0
      [MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, cost, weights, uv, sli] = deal(varargin{:});

      M = MK + MU;
      
      infoold = - K * log2(K.'/MK);
      infonew = - NA * log2(NA.'/MNA);
      infosplit = - MU * log2(MU/M) - MNA * log2(MNA/M);      

      if ~isempty(R) 
        infolr = - ( ...
          sum(L .* log2(L./(repmat(ML, [1 size(L, 2)]))), 2) + ...
          sum(R .* log2(R./(repmat(MR, [1 size(R, 2)]))), 2) ...
        );
      
        % best criterion value over all thresholds
        [infolr, idx] = min(infolr);
        
        infonew = infonew + infolr;
        infosplit = infosplit - [ML(idx), MR(idx)] * log2([ML(idx); MR(idx)]/M);

      else
        idx = [];
        infonew = infonew - AP * log2(AP.'/MAP);
        infosplit = infosplit - MAP * log2(MAP/M);
      end
      
      infothresh = log2(max(length(uv) - 1, 1));
      infogain = (infoold - infonew - infothresh);

    case 1
      idx = [];
       [MU, MK, K, MB, B, cost, weights] = deal(varargin{:});
      
      M = MK + MU;
      
      infoold = - K * log2(K.'/MK);
      infosplit = - MU.' * log2(MU/M);
      
      infonew = - sum(sum(B .* log2(B./(repmat(MB, [1 size(B, 2)]))), 2), 1);
      infosplit = infosplit - MB.' * log2(MB/M);
      
      infogain = (infoold - infonew);
  end
  
  % infogain = infogain / M;
  % infosplit = infospit / M;
  
  if infogain > 0
    critval = infogain / infosplit;    
  else
    critval = -inf;
  end
  
  return

%GINI 
% 
% 	[CRITVAL, IDX] = GINI(FEATTYPE, MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, COST, WEIGHTS, UV, SLI)
%   [CRITVAL, IDX] = GINI(FEATTYPE,MU, MK, K, MB, B, COST, WEIGHTS% 

function [critval, idx] = gini(feattype, varargin) 
	
  	
  switch feattype
    case 0
      [MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, cost, weights, uv, sli] = deal(varargin{:});
      
      if isempty(cost)
        purityold = K * (K.'/MK);
        %impurityold = MK - purityold;

        if ~isempty(R)
          puritylr = ( ...
            sum(L .* (L./(repmat(ML, [1 size(L, 2)]))), 2) + ...
            sum(R .* (R./(repmat(MR, [1 size(R, 2)]))), 2) ...
          );  

          % best criterion value over all thresholds
          [puritylr, idx] = max(puritylr);
          puritynew = puritylr + NA * (NA.'/MNA);
        
        else
          idx = [];
          puritynew = AP * (AP.'/ MAP) + NA * (NA.'/MNA);
        end
        
        %impuritynew = MK - puritynew;
        deltaimpurity = puritynew - purityold;

      else
        impurityold = K * cost * (K.'/MK);
        
        if ~isempty(R)
          impuritylr = ...
            sum((L*cost) .* (L./(repmat(ML, [1 size(L, 2)]))), 2) + ...
            sum((R*cost) .* (R./(repmat(MR, [1 size(R, 2)]))), 2);
 
          [impuritylr, idx] = min(impuritylr);
          impuritynew = impuritylr + NA * cost * (NA.'/MNA);
        
        else
          idx = [];
          impuritynew = AP * cost* (AP.'/ MAP) + NA * cost *(NA.'/MNA);
        end
        
        deltaimpurity = impurityold - impuritynew;
      end
      
    case 1
      idx = [];
      [MU, MK, K, MB, B, cost, weights] = deal(varargin{:});
      
      if isempty(cost)
        purityold = K * (K.'/MK);
        %impurityold = MK - purityold;        
        puritynew = sum(sum(B .* (B./(repmat(MB, [1 size(B, 2)]))), 2));
        %impuritynew = MK - puritynew;        
        deltaimpurity = puritynew - purityold;
      
      else
        impurityold = K * cost * (K.'/MK);        
        impuritynew = ( ...
          sum(sum((B*cost) .* (B./(repmat(MB, [1 size(B, 2)]))), 2), 1) ...
        );
      
        deltaimpurity = impurityold - impuritynew;
      end
  end
  
  
  if isempty(cost)
    deltamin = 0.1;
  else
    deltamin = 0.1*min(abs(cost(abs(cost) > 0)));
  end
  
  
  if deltaimpurity > deltamin
    critval = deltaimpurity/(MU+MK);
  else
    critval = -inf;
  end
  
  return

%MISCLS Miscalssification error
% 
% 	[CRITVAL, IDX] = MISCLS(FEATTYPE, MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, COST, WEIGHTS, UV, SLI)
%   [CRITVAL, IDX] = MISCLS(FEATTYPE,MU, MK, K, MB, B, COST, WEIGHTS
% 

function [critval, idx] = miscls(feattype, varargin) 
	
  	
  switch feattype
    case 0
      [MU, MK, K, MNA, NA, PRMAP, AP, ML, L, MR, R, cost, weights, uv, sli] = deal(varargin{:});      

      if isempty(cost)
        ccold = max(K, [], 2);
        if ~isempty(R)
          ccnew = max(L, [], 2) + max(R, [], 2);
          [ccnew, idx] = max(ccnew);
          ccnew = ccnew + max(NA, [], 2);
        else
          idx = [];
          ccnew = max(AP, [], 2) + max(NA, [], 2);
        end
        
        deltamc = ccnew - ccold;
      
      else
        mcold = min(K*cost, [], 2);
        if ~isempty(R)
          mcnew = min(L*cost, [], 2) + min(R*cost, [], 2);
          [mcnew, idx] = min(mcnew);     
           mcnew = mcnew + min(NA*cost, [], 2);
        else
          idx = [];
          mcnew = min(AP*cost, [], 2) + min(NA*cost, [], 2);          
        end
        
        deltamc = mcold - mcnew; 
      end

    case 1
      idx = [];
      [MU, MK, K, MB, B, cost, weights] = deal(varargin{:});      

      if isempty(cost)
        ccold = max(K, [], 2);        
        ccnew = sum(max(B, [], 2), 1);
        deltamc = ccnew - ccold;        
      else
        mcold = min(K*cost, [], 2);        
        mcnew  = sum(min(B*cost, [], 2), 1);
        deltamc = mcold - mcnew; 
      end
  end
  

  %if isempty(cost)
  %  deltamin = min(weights(weights > 0));
  %else
  %  deltamin = min(weights(weights > 0))*min(cost(cost > 0));
  %end
  %
  %if ~isempty(deltamin)
  %  deltamin = 0.5*deltamin;
  %else
  %  deltamin = 0;
  %end
  %
  %if deltamc <= deltamin
  %  critval = -inf;
  %else
  %  critval = deltamc / (MU+MK);
  %end
  
  critval = deltamc / (MU+MK);
  
  return

%PRUNEP Pessimistic pruning of a decision tree
% 
% 	TREE = PRUNEP(TREE,NODE)
% 
% Pessimistic pruning defined by Quinlan.

function tree = prunep(tree, node)
	
    
  persistent pt;
  
  if (nargout ~= 0) || (nargin ~= 1) || ~isscalar(tree) || ~isnumeric(tree) || ~isint(tree)
    if nargin < 2 || isempty(node)
      node = 1;
    end
    
    pt = tree;
    prunep(node);
    tree = pt;
    clear pt;
  
  else
    node = tree;
    if pt.fidx(node) > 0
      tnidx = (node+1):(node + pt.size(node) - 1);

      nleaves = nnz(pt.fidx(tnidx) == 0);
      errt = pt.errt(node) + 0.5*nleaves;
      errl = pt.errl(node) + 0.5;
      nsmp = pt.nsmp(node);
      sd = sqrt(errt*(1-errt/nsmp));

      if errl < (errt + sd)
        pt.fidx(tnidx) = -1;

        pt.errt(node) = pt.errl(node);
        pt.size(node) = 1;
        pt.fidx(node) = 0;
        pt.fval(node) = {[]};
        pt.nidx(node) = {[]};

      else
        errt = 0;
        for i=1:length(pt.nidx{node})
          idx = pt.nidx{node}(i);
          prunep(idx);
          errt = errt + pt.errt(idx); 
        end
        pt.errt(idx) = errt;
      end
    end
  end

  return

%PRUNET Prune tree by testset
% 
% 	TREE = PRUNET(TREE,T,COST)
% 
% The test set a is used to prune a decision tree. 

function tree = prunet(tree, t, cost)
	
  	
  persistent pt;
  
  if (nargout ~= 0) || (nargin ~= 1) || ~isscalar(tree) || ~isnumeric(tree) || ~isint(tree)
    if nargin < 3
      cost = [];
    end
    
    [m, k, c] = getsize(t);
	
    cs = classsizes(t);
    prior = getprior(t);
    prior = prior/sum(prior);
    weights = m*prior./cs;

    pt = tree;
    mt = size(pt.post,1);
    [p, s] = mapdt(pt, +t, getnlab(t), weights);
	
    % error (cost) in each node as if there were leafs
    if isempty(cost)
      idx = sub2ind([mt, c], (1:mt).', pt.cidx);  
      pt.test_errl = sum(s,2) - s(idx);
    else
      pt.test_errl = sum(s .* cost(:, pt.cidx).', 2);  
    end
    
    pt.test_errt = pt.test_errl;
    
    prunet(1);
    tree = pt;
    clear pt;

  else
    node = tree;
    if pt.fidx(node) > 0
      test_errt = 0;
      errt = 0;

      for i=1:length(pt.nidx{node})
        idx = pt.nidx{node}(i);
        prunet(idx);
        test_errt = test_errt + pt.test_errt(idx);
        errt = errt + pt.errt(idx);
      end

      if pt.test_errl(node) <= test_errt
        pt.fidx(pt.nidx{node}) = -1;

        pt.errt(node) = pt.errl(node);    
        pt.size(node) = 1;
        pt.fidx(node) = 0;
        pt.fval(node) = {[]};
        pt.nidx(node) = {[]};

      else
        pt.test_errt(node) = test_errt;
        pt.errt(node) = errt;
      end
    end
  end;
  
  return
  
function tree = cleandt(tree)
	
    
  rnidx = tree.fidx == -1;
  if nnz(rnidx) == 0
    return
  end
  
  rncount = cumsum(rnidx);
  fn = fieldnames(tree);
  for i=1:length(fn)
    tree.(fn{i})(rnidx, :) = [];
  end
  
  tree.nidx = cellfun(@(x) x - rncount(x).', tree.nidx, 'UniformOutput', false);

  idx = (1:length(rnidx)).';
  idx(rnidx) = [];
  tree.size = tree.size - (rncount(idx) - rncount(idx + tree.size - 1));

 	return

