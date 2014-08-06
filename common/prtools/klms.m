%KLMS Karhunen Loeve Mapping, followed by scaling
% 
%   [W,FRAC] = KLMS(A,N)
%   [W,N]    = KLMS(A,FRAC)
% 
% INPUT
%   A            Dataset
%   N  or FRAC   Number of dimensions (>= 1) or fraction of variance (< 1) 
%                to retain; if > 0, perform PCA; otherwise MCA. Default: N = inf.
%
% OUTPUT
%   W            Affine Karhunen-Loeve mapping
%   FRAC or N    Fraction of variance or number of dimensions retained.
%
% DESCRIPTION
% First a Karhunen Loeve Mapping is performed (i.e. PCA or MCA on the average 
% prior-weighted class covariance matrix). The result is scaled by the mean 
% class standard deviations. For N and FRAC, see KLM.
%
% Default N: select all ('pre-whiten' the average covariance matrix, i.e.
% orthogonalize and scale). The resulting mapping has a unit average
% covariance matrix.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, KLM, PCA

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: klms.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function [w,truefrac] = klms(a,n)

		if (nargin < 2), n = []; end;
	if (nargin < 1) | (isempty(a))
		w = prmapping('klms',n);
		w = setname(w,'Scaled KL Mapping');
		return
	end
	
	[w,truefrac] = klm(a,n);        % Calculate KL mapping
	b = a*w;                        % Combine KL mapping with scaling on
	w = w*scalem(b,'c-variance');   % KL-mapped data
	w = setname(w,'Scaled KL Mapping');

	return
