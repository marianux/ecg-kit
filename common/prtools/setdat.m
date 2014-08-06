%SETDAT Reset data and feature labels of dataset for classification output
%
%   A = SETDAT(A,DATA,W)
%
% INPUT
%    A      Dataset
%    DATA   Dataset or double  
%    W      Mapping  (optional)
%
% OUTPUT
%    A      Dataset 
%
% DESCRIPTION
% The data in the dataset A is replaced by DATA (dataset or double). The
% number of objects in A and DATA should be equal. Optionally, A is given
% the feature labels and the output size as defined by the the mapping W.
% This call is identical to:
%    A = SETDATA(A,DATA);
%    A = SET(A,'featlab',GETLABELS(W),'featsize',GETSIZE_OUT(W));
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASET, SETDATA

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: setdat.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function a = setdat(a,b,w)
		if ~isdataset(a)
		a = prdataset(a);
	end
		
	a = setdata(a,b);

	if (nargin > 2)
		if (~isa(w,'prmapping'))
			error('Third parameter should be mapping')
		end
		% Special case: 2 classes, p(class2) = 1-p(class1).
	%No! this should not be done here!
		%if (size(a,2) == 1) & (size(w,2) == 2)
		%	a = [a 1-a];
		%end
		% Add attributes of W to A.
    a = setfeatsize(a,getsize_out(w));
    featlab = getlabels(w);
    if ~isempty(featlab);
      a = setfeatlab(a,featlab);
    end
		%a = set(a,'featlab',getlabels(w),'featsize',getsize_out(w));
	end
		
	return;
