%DISNORM Trainable mapping for dissimilarity matrix normalization
%
% 	V = DISNORM(D,OPT)
% 	V = D*DISNORM([],OPT)
% 	V = D*DISNORM(OPT)
%   F = E*V
%
% INPUT
%   D 	 NxN dissimilarity matrix or dataset, which defines the norm
%   E    Matrix to be normalized, e.g. D itself
% 	OPT  'max' : maximum dissimilarity is set to 1 by global rescaling
% 		   'mean': average dissimilarity is set to 1 by global rescaling (default)
%
% OUTPUT
%   V 	Trained mapping
%   F  	Normalized dissimilarity data
%
% DEFAULT
%   OPT = 'mean'
%
% DESCRIPTION
% Operation on dissimilarity matrices, like the computation of classifiers
% in dissimilarity space, may depend on the scaling of the dissimilarities
% (a single scalar for the entire matrix). This routine computes a scaling
% for a giving matrix, e.g. a training set and applies it to other
% matrices, e.g. the same training set or based on a test set.
%
% Note that DNORM = (DISNORM*MAPEX); F = D*DNORM is equivalent to 
% F = D*DISNORM(D). So the user defined DNORM can be treated as a fixed
% mapping for normalizing dissimilarity matrices by themselves.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, MAPEX

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function out = disnorm(varargin)

	argin = shiftargin(varargin,'char');
  argin = setdefaults(argin,[],'mean');
  if mapping_task(argin,'definition')
    % call like U=disnorm or U=disnorm('mean') or U=disnorm([],'mean')
    out = define_mapping(argin,'untrained');
    out = setname(out,'Disnorm');
  elseif mapping_task(argin,'training')
    % call like W=disnorm(D,'mean') or W=D*U
    [D,opt] = deal(argin{:});
    if ~isdataset(D)
      D = prdataset(D,1);
      D = setfeatlab(D,getlabels(D));
    end
    if strcmpi(opt,'mean')
      n = size(D,1);
      m = sum(sum(+D))/(n*(n-1));
    elseif strcmpi(opt,'max')
      m = max(D(:));
    else
      error('Wrong OPT')
    end
    out = prmapping(mfilename,'trained',{m},[],size(D,2),size(D,2));
  elseif mapping_task(argin,'trained execution')
    % call like E=D*W
    [D,W] = deal(argin{:});
    m = getdata(W,1);
    out = D./m;
  else
    error('Illegal input')
  end			
	
