%NLABELD Return numeric labels of classified dataset (c
% 
% 	NLABELS = NLABELD(Z)
% 	NLABELS = Z*NLABELD
%	  NLABELS = NLABELD(A,W)
% 	NLABELS = A*W*NLABELD
%
% INPUT
%		Z        Classified dataset, or
%		A,W      Dataset and classifier mapping
%
% OUTPUT
%		NLABELS	 Column vector of numeric labels)
%
% DESCRIPTION 
% Returns the numberic labels of the classified dataset Z (typically the
% result of a mapping or classification A*W). For each object in Z (i.e. 
% each row) the feature label or class label (i.e. the column label) of the
% maximum column value is returned. This corresponds with the classes
% stored in W, which can be found by GETLABELS(W).
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, TESTC, PLOTC, GETLABELS

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function labels = nlabeld(a,w)

		if (nargin == 0)

		% Untrained mapping.
		labels = prmapping(mfilename,'fixed');

	elseif (nargin == 1)

		% In a classified dataset, the feature labels contain the output
		% of the classifier.
		[m,k] = size(a); featlist = getfeatlab(a);

		if (k == 1)
			% If there is one output, assume it's a 2-class discriminant: 
			% decision boundary = 0. 
			J = 2 - (double(a) >= 0); 
		else
			% Otherwise, pick the column containing the maximum output.
			[dummy,J] = max(+a,[],2);
		end
		labels = J;
	elseif (nargin == 2)

		% Just construct classified dataset and call again.
		labels = feval(mfilename,a*w);

	else
		error ('too many arguments');
	end

return

