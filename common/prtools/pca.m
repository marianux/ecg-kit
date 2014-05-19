%PCA Principal component analysis (PCA or MCA on overall covariance matrix)
% 
%   [W,FRAC] = PCA(A,N)
%   [W,N]    = PCA(A,FRAC)
%
% INPUT
%   A           Dataset
%   N  or FRAC  Number of dimensions (>= 1) or fraction of variance (< 1) 
%               to retain; if > 0, perform PCA; otherwise MCA. Default: N = inf.
%
% OUTPUT
%   W           Affine PCA mapping
%   FRAC or N   Fraction of variance or number of dimensions retained.
%
% DESCRIPTION
% This routine performs a principal component analysis (PCA) or minor
% component analysis (MCA) on the overall covariance matrix (weighted
% by the class prior probabilities). It finds a rotation of the dataset A to 
% an N-dimensional linear subspace such that at least (for PCA) or at most 
% (for MCA) a fraction FRAC of the total variance is preserved.
%
% PCA is applied when N (or FRAC) >= 0; MCA when N (or FRAC) < 0. If N is 
% given (abs(N) >= 1), FRAC is optimised. If FRAC is given (abs(FRAC) < 1), 
% N is optimised. 
%
% Objects in a new dataset B can be mapped by B*W, W*B or by A*PCA([],N)*B.
% Default (N = inf): the features are decorrelated and ordered, but no 
% feature reduction is performed.
%
% ALTERNATIVE
%
%   V = PCA(A,0)
% 
% Returns the cumulative fraction of the explained variance. V(N) is the 
% cumulative fraction of the explained variance by using N eigenvectors.
%
% Use KLM for a principal component analysis on the mean class covariance.
% Use FISHERM for optimizing the linear class separability (LDA).
% 
% SEE ALSO
% MAPPINGS, DATASETS, PCLDC, KLLDC, KLM, FISHERM

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: pca.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function [w,truefrac] = pca (varargin)

	prtrace(mfilename);

	[w,truefrac] = pcaklm(mfilename,varargin{:});

	return
