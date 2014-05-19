function c = nancov(x,varargin)
%NANCOV Covariance matrix, ignoring NaNs.
%   C = NANCOV(X), if X is a vector, returns the sample variance of the
%   values in X, treating NaNs as missing values.  For matrices, where
%   each row is an observation and each column a variable, NANCOV(X) is
%   the covariance matrix computing using rows of X that do not contain
%   any NaN values.  NANCOV(X,Y), where X and Y are matrices with
%   the same number of elements, is equivalent to NANCOV([X(:) Y(:)]). 
%   
%   NANCOV(X) or NANCOV(X,Y) normalizes by (N-1) if N>1, where N is the
%   number of observations after removing missing values.  This makes
%   NANCOV(X) the best unbiased estimate of the covariance matrix if the
%   observations are from a normal distribution. For N=1, COV normalizes
%   by N.
%
%   NANCOV(X,1) or NANCOV(X,Y,1) normalizes by N and produces the second
%   moment matrix of the observations about their mean.  NANCOV(X,Y,0) is
%   the same as NANCOV(X,Y), and NANCOV(X,0) is the same as NANCOV(X).
%
%   C = NANCOV(...,'pairwise') computes C(I,J) using rows with no NaN
%   values in columns I or J.  The result may not be a positive definite
%   matrix. C = NANCOV(...,'complete') is the default, and it omits rows
%   with any NaN values, even if they are not in column I or J.
%
%   The mean is removed from each column before calculating the
%   result.
%
%   Example:  Generate random data having non-zero covariance between
%             column 4 and the other columns.
%       x = randn(30,4);       % uncorrelated data
%       x(:,4) = sum(x,2);     % introduce correlation
%       x(2,3) = NaN;          % introduce one missing value
%       c = nancov(x)          % compute sample covariance
%
%   Class support for inputs X,Y:
%      float: double, single
%
%   See also COV, VAR, NANVAR.

%   Copyright 1984-2011 The MathWorks, Inc.
%   $Revision: 1.1.8.3 $  $Date: 2011/02/09 19:35:29 $

if nargin<1
   error(message('stats:nancov:NotEnoughInputs'));
end

% Should we ignore NaNs by complete rows or pairwise?
dopairwise = false;
if numel(varargin)>0
   temp = varargin{end};
   if ischar(temp)
      j = find(strncmpi(temp, {'pairwise' 'complete'},length(temp)));
      if isempty(j)
         error(message('stats:nancov:InvalidArg', temp));
      end
      dopairwise = (j==1);
      varargin(end) = [];
   end
end

% Should we use the mle (divide by N) or unbiased estimate (N-1)?
domle = false;
if numel(varargin)>0
   temp = varargin{end};
   if isequal(temp,0) || isequal(temp,1)
      domle = (temp==1);
      varargin(end) = [];
   end
end

if numel(varargin)>1
   error(message('stats:nancov:TooManyArgs'));
end

scalarxy = false; % nancov(scalar,scalar) is an ambiguous case
if numel(varargin)>0
   y = varargin{1};

   % Two inputs, convert to equivalent single input
   x = x(:);
   y = y(:);
   if length(x)~=length(y)
      error(message('stats:nancov:XYmismatch'));
   end
   scalarxy = isscalar(x) && isscalar(y);
   x = [x y];
elseif ndims(x)>2
   error(message('stats:nancov:InputDim'));
end

if isvector(x) && ~scalarxy
  x = x(:);
end

xnan = isnan(x);
[m,n] = size(x);

if isempty(x);
  if (m==0 && n==0)
      c = NaN(class(x));
  else
      c = NaN(n,class(x));
  end
  return;
end

if ~dopairwise || ~any(any(xnan))    % no need to do pairwise
   nanrows = any(xnan,2);
   if any(nanrows)
       x = x(~nanrows,:);
   end
   c = localcov(x,domle);
else                                         % pairwise with some NaNs
   % Compute variance using complete data separately by column
   c = zeros(n,class(x));
   x(xnan) = 0;
   colsize = sum(~xnan,1);
   xmean = sum(x,1) ./ max(1,colsize);
   xmean(colsize==0) = NaN;
   xc = x - repmat(xmean,m,1);
   xc(xnan) = 0;
   xvar = sum(xc.^2,1);
   if domle
      denom = colsize;
   else
      denom = max(0,colsize-1);
   end
   xvar(denom>1) = xvar(denom>1) ./ denom(denom>1);
   xvar(denom==0) = NaN;
   c(1:n+1:end) = xvar;

   % Now compute off-diagonal entries
   jk = 1:2;
   for j = 2:n
      jk(1) = j;
      for k=1:j-1
         jk(2) = k;
         rowsjk = ~any(xnan(:,jk),2);
         njk = sum(rowsjk);
         if njk<=1
            cjk = NaN;
         else
            cjk = localcov(x(rowsjk,jk),domle);
            cjk = cjk(1,2);
         end
         c(j,k) = cjk;
      end
   end
   c = c + tril(c,-1)';
end

% ------------------------------------------------
function [c,n] = localcov(x,domle)
%LOCALCOV Compute cov with no error checking and assuming NaNs are removed

[m,n] = size(x);
if domle
   denom = m;
else
   denom = max(0,m-1);
end

if m==1   % and doing mle, be sure to get exact 0
   c = zeros(n,class(x));
elseif denom==0
   c = NaN(n,class(x));
else
   xc = x - repmat(mean(x),m,1);
   c = xc' * xc / denom;
end

