%LABELIM Construct image of object (pixel) labels 
%
%		IM = LABELIM(A)
%		IM = A*LABELIM
%
% INPUT
% 	A 	Dataset containing images stored as features
%
% OUTPUT
%		IM	Image containing the labels of the objects
%
% DESCRIPTION
% For a dataset A containing images stored as features, where each pixel
% corresponds to an object, the routine creates an image presenting the
% labels of the objects. Note that if the number of classes is small, e.g.
% 2, an appropriate colormap will have to be loaded for displaying the
% result using IMAGE(LABELS). More appropriate, LABELS should be multiplied
% such that the minimum and maximum of LABELS are well spread in the [1,64]
% interval of the standard colormaps. bb
% The resulting image may also directly	be displayed by:
%		LABELIM(A) or 
%		A*LABELIM
% for which a suitable colormap is loaded automatically.
% 
%	SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%	DATASETS, CLASSIM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: labelim.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function labels = labelim(a)

		% No arguments given: return an untrained mapping.

	if (nargin == 0)
		labels = prmapping('labelim','fixed');
		return
	end

	isfeatim(a);									% Assert that A is a feature image dataset.
	[n,m] = getobjsize(a);				% Get image size and reshape labels to image.
	J = getnlab(a); labels = reshape(J,n,m);

	if (nargout == 0)
		n = 61/(size(a,2)+0.5);			% If no output is requested, display the
		imagesc(labels*n);						%    image with a suitably scaled colormap.
		colormap colorcube;
		clear labels;
	end

return
