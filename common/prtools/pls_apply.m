%pls_apply  Partial Least Squares (applying)
%
%  Y = pls_apply(X,B)
%  Y = pls_apply(X,B,Options)
%
% INPUT
%  X       [N -by- d_X]    the input  data matrix, N samples, d_X variables
%  B       [d_X -by- d_Y]  regression matrix: Y_new = X_new*B
%                         (X_new here after preprocessing, Y_new before
%                         un-preprocessing; preprocessing and
%                         un-preprocessing could be done automatically
%                         (than Options contains info about
%                         preprocessing) or manually)
%  Options  structure returned by pls_train (if not supplied then will
%  be no preprocessing performed)  
%
% OUTPUT
%  Y [N -by- d_Y]    the output data matrix, N samples, d_Y variables
%
% DESCRIPTION
% Applys PLS (Partial Least Squares) regression model
%
% SEE ALSO
% pls_train

% Copyright: S.Verzakov, serguei@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: pls_apply.m,v 1.1 2007/08/28 11:00:39 davidt Exp $

function Y = pls_apply(X,B,Options)

if nargin < 3
  Options  = [];
end

DefaultOptions.X_centering = [];
DefaultOptions.Y_centering = [];
DefaultOptions.X_scaling = [];
DefaultOptions.Y_scaling = [];

Options = pls_updstruct(DefaultOptions, Options);

[N, d_X]    = size(X);
[d_XB, d_Y, M] = size(B);

if d_X ~= d_XB
  error('size(X,2) must be equal to size(B,1)');
end

X = pls_prepro(X, Options.X_centering, Options.X_scaling);
Y = zeros(N,d_Y,M);
for i=1:M
  Y(:,:,i) = pls_prepro(X*B(:,:,i), Options.Y_centering, Options.Y_scaling, -1);
end

return;

