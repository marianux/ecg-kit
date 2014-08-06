%CLASSIM Classify image and return resulting label image
%
%   LABELS = CLASSIM(Z)
%   LABELS = CLASSIM(A,W)
%   LABELS = A*W*CLASSIM
%
% INPUT
%   Z      Classified dataset, or
%   A,W    Dataset and classifier mapping
%
% OUTPUT
%   LABELS Label image
%          When no output is requested, the label image is displayed.
%
% DESCRIPTION
% Returns an image with the labels of the classified dataset image Z
% (typically, the result of a mapping or classification A*W in which A is 
% a set of images stored as features using IM2FEAT). For each object in
% Z (a pixel), a numeric class label is returned. The colormap is loaded
% automatically.
%
% Note that if the number of classes is small, e.g. 2, an appropriate 
% colormap has to be loaded for displaying the result by IMAGE(LABELS), 
% or more appropriately, LABELS should be multiplied such that the minimum 
% and maximum of LABELS are well spread in the [1,64] interval of the 
% standard colormaps.
%
% EXAMPLES
% PREX_SPATM
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, IM2FEAT, LABELD

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: classim.m,v 1.3 2007/09/09 21:21:20 duin Exp $

function labels = classim(a,w)

		
	% Untrained mapping
	if nargin == 0
		labels = prmapping('classim','fixed');
		return
	end
	
	if nargin == 2
		ismapping(w);
		a = a*w;
	end

	% Assertion: generate an error, if a is not a dataset with objects as pixels
	isfeatim(a);
	
	if size(a,2) == 1 
		% Assuming the 2-class case
		J = 2 - (double(a) >= 0);
	else
		% Multi-class problem
		[mx,J] = max(double(a),[],2);
	end

	%fl = renumlab(getfeatlab(a));
	%labels = reshape(fl(J),getobjsize(a));
	labels = reshape(J,getobjsize(a));

	% Display the label image	
	if nargout == 0
		n = 61/(size(a,2) +0.5);
		imagesc(labels*n)
		colormap colorcube
		clear labels
	end

return;
