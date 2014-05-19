%DISNORM Normalization of a dissimilarity matrix
%
% 	V = DISNORM(D,OPT)
%   F = E*V
%
% INPUT
%   D 	 NxN dissimilarity matrix or dataset, which sets the norm
%   E    Matrix to be normalized, e.g. D itself
% 	OPT  'max' : maximum dissimilarity is set to 1 by global rescaling
% 		   'mean': average dissimilarity is set to 1 by global rescaling (default)
%
% OUTPUT
%   V 	Fixed mapping
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

% Copyright: Elzbieta Pekalska, ela.pekalska@googlemail.com
% Faculty EWI, Delft University of Technology and
% School of Computer Science, University of Manchester

function V = disnorm(D,opt)
if nargin < 2,
	opt = 'mean';
end

if nargin == 0 | isempty(D)
  V = prmapping(mfilename,{opt});
  V = setname(V,'Disnorm');
  return
end

if ~isdataset(D)
	D = prdataset(D,1);
	D = setfeatlab(D,getlabels(D));
end

%DEFINE mapping
if isstr(opt)
%	discheck(D);
	opt = lower(opt);
	if strcmp(opt,'mean')
		n = size(D,1);
		m = sum(sum(+D))/(n*(n-1));
	elseif strcmp(opt,'max')
		m = max(D(:));
	else
		error('Wrong OPT.')
	end
	if nargout > 1
		D = D./m;
	end	
    V = prmapping(mfilename,'trained',{m},[],size(D,2),size(D,2));
	return;
end

% APPLY mapping
if ismapping(opt)
	opt = getdata(opt,1);
	V = D./opt;
	return;
end			
	
