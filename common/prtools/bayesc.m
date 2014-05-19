%BAYESC Bayes classifier
%
%		W = BAYESC(WA,WB, ... ,P,LABLIST)
%
% INPUT
%   WA, WB, ... Trained mappings for supplying class density estimates
%   P           Vector with class prior probabilities 
%               Default: equal priors
%   LABLIST    	List of class names (labels)
%
% OUTPUT
%    W          Bayes classifier.
%
% DESCRIPTION
% The trained mappings WA,WB, ... should supply proper densities estimates
% D for a dataset X by D = X*WA, etcetera. E.g. they should be trained by
% commands like GAUSSM(A), PARZENM(A), KNNM(A). Consequently, they should
% have a size of K x 1 (assuming that X and A are K-dimensional). Also
% sizes of K x N are supported, assuming a combined density estimate for N
% classes simultaneously. BAYESC weighs the class densitites by the class
% priors in P and names the classes by LABLIST. If LABLIST is not supplied,
% the labels stored in the mappings are used.
%
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd edition, 
% John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 2002.
%
% SEE ALSO
% DATASETS, MAPPINGS, GAUSSM, PARZENM, KNNM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function w = bayesc(varargin)

	
	n = nargin;
	p = [];
	lablist = [];
	
	if (nargin > 1)
		if nargin > 2 & (~ismapping(varargin{end-1}))
			p = varargin{end-1};
			lablist = varargin{end};
			n = n-2;
		elseif (~ismapping(varargin{end}))
			p = varargin{end};
			n = n-1;
		end
	end
	
	if (nargin < 1 | isempty(varargin{1}))
							% Definition
		w = prmapping(mfilename,'combiner',{p,lablist});
		w = setname(w,'Bayes Classifier');
	
	elseif ismapping(varargin{1})
							% Construction of the trained Bayes Classifier						
		w = [];
		k = size(varargin{1},1);
        wsize = 0;
		for j=1:n
			v =  varargin{j};
			if ~ismapping(v) | getout_conv(v) ~= 0
				error('Density estimating mapping expected and not found')
			end
			if size(v,1) ~= k
				error('Mappings / density estimators should be defined for the same dimensionality')
			end
			w = [w v];
            wsize = wsize + size(v,2);
        end
        w = setsize(w,[k,wsize]);

		c = size(w,2);
		if isempty(p), p = ones(1,c)/c; end
		if length(p) ~= c
			error('Vector with prior probabilities has wrong length')
		end

		if isempty(lablist)
			lablist = getlabels(w);
		end
		if size(lablist,1) ~= c
			error('Label list has wrong size')
		end
		w = w*affine(p(:)');
		w = setlabels(w,lablist);
		w = setname(w,'Bayes Classifier');
		
	else
		
		error('Wrong input')
		
	end
		

