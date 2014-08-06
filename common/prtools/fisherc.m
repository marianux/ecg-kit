%FISHERC Trainable classifier: Fisher's Least Square Linear Discriminant
% 
%   W = FISHERC(A)
%   W = A*FISHERC
% 
% INPUT
%   A  Dataset
%
% OUTPUT
%   W  Fisher's linear classifier 
%
% DESCRIPTION  
% Finds the linear discriminant function between the classes in the 
% dataset A by minimizing the errors in the least square sense. This 
% is a multi-class implementation using the one-against-all strategy.
%
% This classifier also works for soft and  target labels.
%
% For high dimensional datasets or small sample size situations, the 
% Pseudo-Fisher procedure is used, which is based on a pseudo-inverse.
%
% This classifier, like all other non-density based classifiers, does not
% use the prior probabilities stored in the dataset A. Consequently, it
% is just for two-class problems and equal class prior probabilities 
% equivalent to LDC, which assumes normal densities with equal covariance
% matrices.
%
% Note that A*(KLMS([],N)*NMC) performs a very similar operation, but uses
% the prior probabilities to estimate the mean class covariance matrix used
% in the pre-whitening operation performed by KLMS. The reduced
% dimensionality N controls some regularization.
% 
% REFERENCES
% 1. R.O. Duda, P.E. Hart, and D.G. Stork, Pattern classification, 2nd ed.
% John Wiley and Sons, New York, 2001.
% 2. A. Webb, Statistical Pattern Recognition, Wiley, New York, 2002.
% 3. S. Raudys and R.P.W. Duin, On expected classification error of the
% Fisher linear classifier with pseudo-inverse covariance matrix, Pattern
% Recognition Letters, vol. 19, no. 5-6, 1998, 385-392.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, TESTC, LDC, NMC, FISHERM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: fisherc.m,v 1.6 2010/02/08 15:31:48 duin Exp $

function W = fisherc(a)
		
	% No input arguments, return an untrained mapping.
	if (nargin < 1) | (isempty(a))
		W = prmapping(mfilename); 
		W = setname(W,'Fisher');
		return;
	end

	isvaldfile(a,1,2); % at least 1 object per class, 2 classes
	a = testdatasize(a); 
		
	[m,k,c] = getsize(a);
	if islabtype(a,'crisp') & c > 2
		W = mclassc(a,fisherc);
		return
	end
	y = gettargets(a);
	if islabtype(a,'soft')
		y = invsigm(y); % better to give [0,1] targets full range
	end
	u = mean(a);    
	% Shift A to the origin. This is not significant, just increases accuracy.
	% A is extended by ones(m,1), a trick to incorporate a free weight in the 
	% hyperplane definition. 
	b = [+a-repmat(u,m,1), ones(m,1)]; 

	if (rank(b) <= k)
		% This causes Fisherc to be the Pseudo-Fisher Classifier
		prwarning(2,'The dimensionality is too large. Pseudo-Fisher is trained instead.');  
		v = prpinv(b)*y;                 
	else
		% Fisher is identical to the Min-Square-Error solution.		
		v = b\y;               
	end

	offset = v(k+1,:) - u*v(1:k,:); 	% Free weight. 
	W = affine(v(1:k,:),offset,a,getlablist(a),k,c);
	% Normalize the weights for good posterior probabilities.
 	W = cnormc(W,a);									
	W = setname(W,'Fisher');

return;
