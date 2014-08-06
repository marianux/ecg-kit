%KLM Karhunen-Loeve Mapping (PCA or MCA of mean covariance matrix)
% 
% 	[W,FRAC] = KLM(A,N)
% 	[W,N]    = KLM(A,FRAC)
%
% INPUT
%   A	         Dataset
%   N	or FRAC  Number of dimensions (>= 1) or fraction of variance (< 1) 
%              to retain; if > 0, perform PCA; otherwise MCA.
%              Default: N = inf.
%
% OUTPUT
%   W          Affine Karhunen-Loeve mapping
%   FRAC or N  Fraction of variance or number of dimensions retained.
%
% DESCRIPTION
% The Karhunen-Loeve Mapping performs a principal component analysis
% (PCA) or minor component analysis (MCA) on the mean class covariance
% matrix (weighted by the class prior probabilities). It finds a
% rotation of the dataset A to an N-dimensional linear subspace such
% that at least (for PCA) or at most (for MCA) a fraction FRAC of the
% total variance is preserved.
%
% PCA is applied when N (or FRAC) >= 0; MCA when N (or FRAC) < 0. If N
% is given (abs(N) >= 1), FRAC is optimised. If FRAC is given
% (abs(FRAC) < 1), N is optimised. 
%
% Objects in a new dataset B can be mapped by B*W, W*B or by
% A*KLM([],N)*B. Default (N = inf): the features are decorrelated and
% ordered, but no feature reduction is performed.
%
% ALTERNATIVE
%
% 	V = KLM(A,0)
% 
% Returns the cummulative fraction of the explained variance. V(N) is
% the cumulative fraction of the explained variance by using N
% eigenvectors.
%
% Use PCA for a principal component analysis on the total data
% covariance. Use FISHERM for optimizing the linear class
% separability (LDA).
%
% This function is basically a wrapper around pcaklm.m.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, PCAKLM, PCLDC, KLLDC, PCAM, FISHERM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: klm.m,v 1.2 2006/03/08 22:06:58 duin Exp $

function [w,truefrac] = klm (varargin)

		[w,truefrac] = pcaklm(mfilename,varargin{:});

return
