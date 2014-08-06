%SPATM Augment image dataset with spatial label information
%
%   E = SPATM(D,S)
%   E = D*SPATM([],S)
% 
% INPUT
%        D     image dataset classified by a classifier
%        S     smoothing parameter (optional, default: sigma = 1.0)
%
% OUTPUT
%        E     augmented dataset with additional spatial information
%
% DESCRIPTION
% If D = A*W*CLASSC, the output of a classification of a dataset A
% containing feature images, then E is an augmented version of D:
% E = [D T]. T contains the spatial information in D, such that
% it adds for each class of which the objects in D are assigned to,
% a Gaussian convoluted (std. dev s) 0/1 image with '1'-s on the
% pixel positions (objects) of that class. T is normalized such that
% its row sums are 1. It thereby effectively contains Parzen estimates
% of the posterior class probabilities if the image is considered as a
% feature space. Default: S = 1.
%
% Spatial and feature information can be combined by feeding E into
% a class combiner, e.g: A*W*CLASSC*SPATM([],2)*MAXC
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, PREX_SPATM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: spatm.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function b = spatm(a,s)
	
		
	if nargin < 2, s = 1; end
	if nargin < 1 | isempty(a)
		b = prmapping('spatm','fixed',s);
		return
	end
	
	% assertion: the image with pixels being objects is required
	isfeatim(a);
	
	% initialize the label-image y:
	[m,k,c] = getsize(a);
	[n1,n2] = getobjsize(a);
	%DXD Avoid that the feature labels might be reordered...
	%[labt,x,newlablist] = renumlab(getfeatlab(a),labeld(a));
	[dummy,x] = max(a,[],2);
	y = zeros(n1,n2,max(x));
	y((x(:)-1)*n1*n2 + [1:n1*n2]') = ones(n1,n2);
	% make the label-image a dataset:
	z = im2feat(y);
	% Store also all the useful things like prior or lablist
	z = setdat(a,z);
	% smooth the label-image and add it to the dataset a:
	b = [a datgauss(z,s)];
	%b = datgauss(z,s);

return
