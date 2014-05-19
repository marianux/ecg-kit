%SUBSC Subspace Classifier
%
%   W = SUBSC(A,N)
%   W = SUBSC(A,FRAC)
%
% INPUT
%   A          Dataset
%   N or FRAC  Desired model dimensionality or fraction of retained 
%              variance per class
%
% OUTPUT
%   W          Subspace classifier  
%
% DESCRIPTION
% Each class in the trainingset A is described by linear subspace of
% dimensionality N, or such that at least a fraction FRAC of its variance
% is retained. This is realised by calling PCAM(AI,N) or PCAM(AI,FRAC) for
% each subset AI of A (objects of class I). For each class a model is
% built that assumes that the distances of the objects to the class
% subspaces follow a one-dimensional distribution. 
%
% New objects are assigned to the class of the nearest subspace.
% Classification by D = B*W, in which W is a trained subspace classifier
% and B is a testset, returns a dataset D with one-dimensional densities
% for each of the classes in its columns.
%
% If N (ALF) is NaN it is optimised by REGOPTC.
%
% REFERENCE
% E. Oja, The Subspace Methods of Pattern Recognition, Wiley, New York, 1984.
%
% SEE ALSO
% DATASETS, MAPPINGS, PCAM, FISHERC, FISHERM, GAUSSM, REGOPTC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function W = subsc(A,N)

  name = 'SubspaceC';
  
  % handle default
  if nargin < 2, N = 1; end
  
  % handle untrained calls like subsc([],3);
  if nargin < 1 | isempty(A)
    W = prmapping(mfilename,{N});
    W = setname(W,name);
    return
  end
  
	
	if isa(N,'double') & isnan(N)    % optimize regularisation parameter
		defs = {1};
		parmin_max = [1,size(A,2)];
		W = regoptc(A,mfilename,{N},defs,[1],parmin_max,testc([],'soft'),0);
		
	elseif isa(N,'double')
    
  % handle training like A*subsc, A*subsc([],3), subsc(A)
  % PRTools takes care that they are all converted to subsc(A,N)
      
    islabtype(A,'crisp');     % allow crisp labels only
    isvaldfile(A,1,2);         % at least one object per class, two objects
		A = testdatasize(A,'features');
		A = setprior(A,getprior(A));
    [m,k,c] = getsize(A);     % size of the training set
    for j = 1:c               % run over all classes
      B = seldat(A,j);        % get the objects of a single class only
      u = mean(B);            % compute its mean
      B = B - repmat(u,size(B,1),1); % subtract mean
      v = pcam(B,N);           % compute PCA for this class
      v = v*v';               % trick: affine mappings in the original space
      B = B - B*v;            % differences of objects and their mappings
      s = mean(sum(B.*B,2));  % mean square error w.r.t. the subspace
      data(j).u = u;          % store mean
      data(j).w = v;          % store mapping
      data(j).s = s;          % store mean square distance
    end
                              % define trained mapping,
                              % store class labels and size
    W = prmapping(mfilename,'trained',data,getlablist(A),k,c);
    W = setname(W,name);
    
  elseif isa(N,'prmapping')
    
  % handle evaluation of a trained subspace classifier W for a dataset A.
  % The command D = A*W is by PRTools translated into D = subsc(A,W)
  % Such a call is detected here by the fact that N appears to be a mapping.
    
    W = N;                    % avoid confusion: call the mapping W
    m = size(A,1);            % number of test objects
    [k,c] = size(W);          % mapping size: from K features to C classes
    d = zeros(m,c);           % output: C class densities for M objects
    
    for j=1:c                 % run over all classes
      u = W.data(j).u;        % class mean in training set
      v = W.data(j).w;        % mapping to subspace in original space
      s = W.data(j).s;        % mean square distance
      B = A - repmat(u,m,1);  % substract mean from test set
      B = B - B*v;            % differences objects and their mappings
      d(:,j) = sum(B.*B,2)/s; % convert to distance and normalise 
    end
    d = exp(-d/2)/sqrt(2*pi); % convert to normal density
    
    A = prdataset(A);           % make sure A is a dataset 
    d = setdata(A,d,getlabels(W)); % take data from D and use 
                              % class labels as given in W
                              % other information in A is preserved   
    W = d;                    % return result in output variable W
    
  else     
    
    error('Illegal call')     % this should not happen
    
  end
  
return
