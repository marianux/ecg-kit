%REJECTC Construction of a rejecting classifier
%
%   WR = REJECTC(A,W,FRAC,TYPE)
%
% INPUT
%   A      Dataset
%   W      Trained or untrained classifier
%   FRAC   Fraction to be rejected. Default: 0.05
%   TYPE   String with reject type: 'ambiguity' or 'outlier'.
%          'a' and 'o' are supported as well. Default is 'a'.
%
% OUTPUT
%   WR     Rejecting classifier
%
% EXAMPLE
% a = gendatb
% w = ldc(a);
% v = rejectc(a,w,0.2);
% scatterd(a);
% plotc(w);
% plotc(v,'r')
%
% DESCRIPTION
% This command extends an arbitrary classifier with a reject option. If WR
% is used for classifying a dataset B, then D = B*WR has C+1 columns
% ('features'), one for every class in A and an additional one that takes
% care of the rejection: a NaN for numeric labels (classnames in A) or en
% empty string for string labels. 
% NOTE: Objects that are rejected are not counted as an error in TESTC. The
% classification error estimated by TESTC just considers the total number
% of objects for wich B*WR*LABELD has a correct classname and neglects all
% others. So by rejection the error estimate by TESTC may increase,
% decrease or stay equal.
%
% SEE ALSO
% DATASETS, MAPPINGS, LABELD, TESTC, REJECTM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w_out = rejectc(a,w,frac,type)

if nargin < 4, type = 'a'; end
if nargin < 3, frac = 0.05; end

ismapping(w);
isdataset(a);

if isuntrained(w), w = a*w; end
istrained(w);

conv = getout_conv(w);
if isstr(getlablist(a))
	rejname = '';
else
	rejname = NaN;
end

switch(type)
	case{'a','ambiguity'}
		w_out = w*classc*rejectm(a*w*classc,frac,rejname);
	case{'o','outlier'}
		conv = getout_conv(w);
		if conv > 1
			w = setout_conv(w,conv-2);
		elseif conv == 1
			error('Outlier reject not (yet) supported for this classifier')
		end
		w_out = w*rejectm(a*w,frac,rejname);
	otherwise
		error('Unknown reject type')
end

	
