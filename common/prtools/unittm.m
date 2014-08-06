%UNITTM Trainable unit mapping
% 
%   W = UNITTM(A)
%   W = A*UNITTM
%
% INPUT
%   A   Array or dataset
%
% OUTPUT
%   W   Unit mapping, if applied to dataset, it is returned unchanged
%
% DESCRIPTION
% This is a trainable unit mapping that maps any dataset on itself. 
% The difference with the fixed mapping is UNITM is that UNITM(A) returns
% A and that W = UNITTM(A) returns a mapping with input and output 
% dimensionality the feature size of A. Consequently B*W returns B, but
% checks the feature size against A.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, UNITM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = unittm (a,v)

		if (nargin == 0) | (isempty(a))     % untrained definition
		w = prmapping(mfilename,'untrained');
		w = setname(w,'Unit Mapping');
	elseif (nargin == 1) | isempty(v)   % training
		isdataset(a);
		[m,k] = size(a);
		w = prmapping(mfilename,'trained',[],getfeatlab(a),k,k);
	elseif ismapping(v)                 % execution
		w = a;                            % return input  
	else
		error('Unparsable call')
	end

	return
