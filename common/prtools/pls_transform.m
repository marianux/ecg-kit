%PLS_TRANSFORM  Partial Least Squares transformation
%
%  T = PLS_TRANSFORM(X,R)
%  T = PLS_TRANSFORM(X,R,OPTIONS)
%
% INPUT
%  X       input data matrix size [N -by- d_X], N samples, d_X variables
%  R       transformation matrix size [d_X -by- nLV]: T_new = X_new*R
%          (X_new here after preprocessing, preprocessing and un-preprocessing 
%          could be done automatically (than OPTIONS contains info about
%          preprocessing) or manually); normally, R as a field of XRes
%          output parameter of pls_train routine
%  OPTIONS structure returned by pls_train (if not supplied then will be
%          no preprocessing performed) 
%
% OUTPUT
%  T [N -by- nLV]   scores -- transformed data
%
% DESCRIPTION
% Applys PLS (Partial Least Squares) regression model
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLS_TRAIN, PLS_APPLY

% Copyright: S.Verzakov, serguei@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: pls_transform.m,v 1.1 2007/08/28 11:00:39 davidt Exp $

function T = pls_transform(X,R,Options)

if nargin < 3
  Options  = [];
end

DefaultOptions.X_centering = [];
DefaultOptions.Y_centering = [];
DefaultOptions.X_scaling = [];
DefaultOptions.Y_scaling = [];

Options = pls_updstruct(DefaultOptions, Options);

[N, d_X]    = size(X);
[d_XR, nLv] = size(R);

if d_X ~= d_XR
  error('size(X,2) must be equal to size(R,1)');
end

T = pls_prepro(X, Options.X_centering, Options.X_scaling)*R;

return;



























