%LOGLC Trainable logistic linear classifier
% 
%   W = LOGLC(A)
%   W = A*LOGLC
% 
% INPUT
%   A   Dataset
%
% OUTPUT
%   W   Logistic linear classifier 
%
% DESCRIPTION  
% Computation of the linear classifier for the dataset A by maximizing the
% likelihood criterion using the logistic (sigmoid) function.
% This routine becomes very slow for feature sizes above 1000.
% 
% REFERENCES
% A. Webb, Statistical Pattern Recognition, John Wiley & Sons, New York, 2002.
% J. A. Anderson, Logistic discrimination, in: P. R. Krishnaiah and L. N.
% Kanal (eds.), Handbook of Statistics 2: Classification, Pattern Recognition 
% and Reduction of Dimensionality, North Holland, Amsterdam, 1982, 169--191.
%
%  SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
%  MAPPINGS, DATASETS, LDC, FISHERC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: loglc.m,v 1.8 2010/02/08 15:31:48 duin Exp $

function W = loglc(a)

		% No input data, return an untrained classifier.
	if (nargin == 0) | (isempty(a))
		W = prmapping(mfilename); 
		W = setname(W,'Logistic');
		return
	end

	islabtype(a,'crisp');
	isvaldfile(a,1,2); % at least 2 object per class, 2 classes
	a = testdatasize(a); 
	fid = []; % Progress messages default destination
	[m,k,c] = getsize(a);
	nlab = getnlab(a);
	prior = getprior(a);
	a = setprior(a,prior);
	if (c > 2)
		% Compute C classifiers: each class against all others.	
		W = mclassc(a,prmapping(mfilename));	
	else
		v = scalem(a,'variance');
		a = a*v;
		accuracy = 0.0001;		% An accuracy for the ML loop.
		x = [a,ones(m,1)];
		% A standard trick to set the labels to +1 for the first class
		% (NLAB=1) and to -1 for the second one (NLAB=2). Then, each 
		% object vector is multiplied by its new label +/-1.
		x(find(nlab==2),:) = -x(find(nlab==2),:);
		x = +x;
		alf = sum(nlab==2)/sum(nlab==1);
		weights = zeros(1,k+1);

		% Maximize the likelihood L to find WEIGHTS
		L = -inf; Lnew = -realmax;
		prwaitbar(100,'loglc: Optimizing log likelihoods',k > 100)
		d0 = 10*log(10);
		while (abs(Lnew - L) > accuracy)
			prwaitbar(100, 100-100*(log(abs(Lnew-L))-log(accuracy))/d0);
			pax = ones(m,1) ./ (1 + exp(-x*weights'));	% Estimate of P(class +1|x).
			pbx = 1 - pax;                          		% Estimate of P(class -1|x).
			L = Lnew; Lnew = sum(log(pax+realmin));     % Update likelihood.	
			p2x = sqrt(pax.*pbx); 
			y = x .* p2x(:,ones(1,k+1));
			%size(y'*y)
			weights = pbx' * x * prpinv(y'*y) + weights;
		end
		prwaitbar(0);
		
		% Define LOGLC by an affine (linear) mapping based 
		% on the [K x 1] matrix R and the offset w0.
		w0 = weights(k+1) + log(alf*prior(1)/prior(2));
		R  = weights(1:k)';
		%DXD: for a two-class classifier we have to supply two posterior
		%probabilities:
		%W  = v*affine(R,w0,a,getlablist(a));  % wrong
		W  = v*affine([R -R],[w0 -w0],a,getlablist(a)); 
		W  = setout_conv(W,1);
	end
	W = setname(W,'Logistic');

	return
	
	function [weights,L,Lnew,len2] = optimw(x,weights,L,Lnew,accuracy,len2,fid)
		% this function is never called
		[m,k] = size(x);
		while (abs(Lnew - L) > accuracy)
			pax = ones(m,1) ./ (1 + exp(-x*weights'));	% Estimate of P(class +1|x).
			pbx = 1 - pax;                          		% Estimate of P(class -1|x).
			L = Lnew; Lnew = sum(log(pax+realmin));     % Update likelihood.	
			p2x = sqrt(pax.*pbx); 
			y = x .* p2x(:,ones(1,k));
			weights = pbx' * x * prpinv(y'*y) + weights;
		end
	
	
