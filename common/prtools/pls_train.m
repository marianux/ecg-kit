%pls_train Partial Least Squares (training)
%
%  [B,XRes,YRes,Options] = pls_train(X,Y)
%  [B,XRes,YRes,Options] = pls_train(X,Y,Options)
%
% INPUT
%  X   [N -by- d_X]  the training (input)  data matrix, N samples, d_X variables
%  Y   [N -by- d_Y]  the training (output) data matrix, N samples, d_Y variables
%
%   Options.
%   maxLV        maximal number of latent variables (will be corrected
%                if > rank(X)); 
%                maxLV=inf means maxLV=min(N,d_X) -- theoretical maximum
%                number of LV; 
%                by default =inf
%   method       'NIPALS' or 'SIMPLS'; by default ='SIMPLS'
%
%   X_centering  do nothing (=[] or 0), do mean centering (=nan), center
%                around some vaue v (=v); 
%                by default  =[]
%   Y_centering  do nothing (=[] or 0), do mean centering (=nan), center
%                around some vaue v (=v); 
%                by default  =[]
%   X_scaling    do nothing (=[] or 1), divide each col by std (=nan),
%                divide by some v (=v); 
%                by default =[]
%   Y_scaling    do nothing (=[] or 1), divide each col by std (=nan),
%                divide by some v (=v); 
%                by default =[]
%
% OUTPUT
%  B      [d_X -by- d_Y -by- nLV]  collection of regression matrices:
%         Y_new = X_new*B(:,:,n) represents regression on the first n
%         latent variables (X_new here after preprocessing, Y_new before
%         un-preprocessing)
%  XRes.
%   ssq  [1 -by- nLV]    the part of explaind sum of squares of
%        (preprocessed) X matrix
%   T    [N -by- nLV]    scores   (transformed (preprocessed) X)
%   R    [d_X -by- nLV]  weights  (transformation matrix)
%   P    [d_X -by- nLV]  loadings (back-transformation matrix)
%                        P(:,k) are the regression coef of
%                        (preprocessed) X on T(:,k)
%   W    [d_X -by- nLV]  local weights (local transformation matrix);
%                        ONLY FOR NIPALS
%   V    [d_X -by- nLV]  first k columns of this matrix are the
%                        orthornormal basis of the space spanned by k
%                        first columns of P; ONLY FOR SIMPLS
%  YRes.
%   ssq  [1 -by- nLV]    the part of explained sum of squares of
%                        (preprocessed) Y matrix
%   U    [N -by- nLV]    scores   (transformed (preprocessed) Y)
%   Q    [d_Y -by- nLV]  weights  (transformation matrix)
%   C    [d_Y -by- nLV]  C(:,k) are the regression coeff of
%                        (preprocessed) Y on T(:,k)
%   bin  [1 -by- nLV]    bin(k) is the regression coeff of U(:,k) on T(:,k)
%
% Options. contains the same fields as in input, but the values can be changed 
% (e.g. after mean centering X_centering = mean(X,1))
%
% DESCRIPTION
% Trains PLS (Partial Least Squares) regression model
%
% Relations between matrices (X end Y are assumed to be preprocessed):
% NIPALS:
% T = X*R (columns of T are orthogonal)
% R = W*prinv(P'*W)
% P = X'*T*prinv(T'*T)
% U = Y*Q - T*(C'*Q-tril(C'*Q))
% C = Y'*T*prinv(T'*T) = Q*diag(bin)
% bin = sqrt(diag(T'*Y*Y'*T))'*prinv(T'*T))
% B = R*C' = W*prinv(P'*W)*C'
%
% SIMPLS:
% T = X*R (columns of T are orthonormal)
% P = X'*T
% U = Y*Q
% C = Y'*T = Q*diag(bin)
% bin = sqrt(diag(T'*Y*Y'*T))'
% B = R*C'
%
% BOTH:
% T_new = X_new*R
% Y_new = X_new*B
%  
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PLS_APPLY, PLS_TRANSFORM

% Copyright: S.Verzakov, s.verzakov@ewi.tudelft.nl 
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: pls_train.m,v 1.2 2010/02/08 15:29:48 duin Exp $

function [B, XRes, YRes, Options] = pls_train(X,Y,Options)

[N_X, d_X] = size(X);
[N_Y, d_Y] = size(Y);

if N_X ~= N_Y
  error('size(X,1) must be equal to size(Y,1)');
else
  N = N_X;
end

if nargin < 3
  Options  = [];
end

DefaultOptions.X_centering = [];
DefaultOptions.Y_centering = [];
DefaultOptions.X_scaling = [];
DefaultOptions.Y_scaling = [];
DefaultOptions.maxLV = inf;
DefaultOptions.method = 'SIMPLS';

Options = pls_updstruct(DefaultOptions, Options);

if isinf(Options.maxLV)
  Options.maxLV = min(N,d_X);
elseif Options.maxLV > min(N,d_X)
  error('PLS: The number of LV(s) cannot be greater then min(N,d_X)');
end

[X, Options.X_centering, Options.X_scaling] = pls_prepro(X, Options.X_centering, Options.X_scaling);
[Y, Options.Y_centering, Options.Y_scaling] = pls_prepro(Y, Options.Y_centering, Options.Y_scaling);

ssq_X = sum(X(:).^2);
ssq_Y = sum(Y(:).^2);

B = zeros(d_X,d_Y,Options.maxLV);
XRes.ssq = zeros(1,Options.maxLV);
XRes.T   = zeros(N,Options.maxLV);
XRes.R   = zeros(d_X,Options.maxLV);
XRes.P   = zeros(d_X,Options.maxLV);
XRes.W   = [];
XRes.V   = [];

YRes.ssq = zeros(1,Options.maxLV);  
YRes.U   = zeros(N,Options.maxLV);  
YRes.Q   = zeros(d_Y,Options.maxLV);
YRes.C   = zeros(d_Y,Options.maxLV);
YRes.bin = zeros(1,Options.maxLV);  

ev = zeros(1,Options.maxLV);
nLV = Options.maxLV;
opts.disp   = 0;
opts.issym  = 1;
opts.isreal = 1;

switch upper(Options.method)
case 'NIPALS'
  XRes.W = zeros(d_X,Options.maxLV);
  for LV = 1:Options.maxLV
    S = X'*Y;
    if d_X <= d_Y
      if d_X > 1
        [w, ev(LV)] = eigs(S*S',1,'LA',opts);
      else
        w = 1;
        ev(LV) = S*S';
      end
      
      t = X*w;
  	  t2 = (t'*t);
      proj = t/t2;
  	  p = X'*proj;

      c = Y'*proj;
      bin = norm(c);  %bin = sqrt(ev)/t2:
      q = c/bin;
      u = Y*q;

    else
      if d_Y > 1
        [q, ev(LV)] = eigs(S'*S,1,'LA',opts);
      else
        q = 1;
        ev(LV) = S'*S;
      end
      u = Y*q;

      w = X'*u;
  	  w = w/norm(w);
  	  t = X*w;
   	  t2 = (t'*t);
      proj = t/t2;
  	  p = X'*proj;
      
      bin = u'*proj;  %bin = sqrt(ev)/t2:
      c = q*bin;
  	end

	  if LV == 1 
      if ev(LV) == 0
        error('PLS: Rank of the covariation matrix X''*Y is zero.');
      end  
	  elseif ev(LV) <= 1e-16*ev(1)
      nLV = LV-1;
      WarnMsg = sprintf(['\nPLS: Rank of the covariation matrix X''*Y is exausted ' ...
                         'after removing %d Latent Variable(s).\n'...
                         'Results only for the %d LV(s) will be returned'],nLV,nLV);
      if exist('prwarning')
        prwarning(1, WarnMsg);
      else
        warning(WarnMsg);
      end
      break;
	  end

    %R = W*prinv(P'*W)
    %the next portion of the code makes use of the fact taht P'*W is the upper triangle matrix:
    % 
    %if LV == 1
    %  InvPTW(1,1) = 1/(p'*w);
    %  r = w*InvPTW(1,1);
    %else
    %  InvPTW(1:LV,LV) = [(-InvPTW(1:LV-1,1:LV-1)*(XRes.P(:,1:LV-1)'*w)); 1]/(p'*w);
    %  r = [XRes.W(:,1:LV-1), w] * InvPTW(1:LV,LV);
    %end  
    % 
    %some optimization of the above code (notice that ones(m,0)*ones(0,n) == zeros(m,n)):
    %
    r = (w - (XRes.R(:,1:LV-1)*(XRes.P(:,1:LV-1)'*w)))/(p'*w);

    B(:,:,LV) = r*c';

  	XRes.ssq(:,LV) = t2*(p'*p)/ssq_X;
  	XRes.T(:,LV)   = t;
	  XRes.R(:,LV)   = r;
  	XRes.P(:,LV)   = p;
  	XRes.W(:,LV)   = w;
         
	  YRes.ssq(:,LV) = t2*(bin.^2)/ssq_Y;
  	YRes.U(:,LV)   = u;
  	YRes.Q(:,LV)   = q;
  	YRes.C(:,LV)   = c;
  	YRes.bin(:,LV) = bin;
	  
  	X = X - t*p';
  	Y = Y - t*c';
  end

  if nLV < Options.maxLV
    XRes.W = XRes.W(:,1:nLV);
  end

case 'SIMPLS'
  XRes.V = zeros(d_X,Options.maxLV);
  S = X'*Y;
  for LV = 1:Options.maxLV
	  if d_X <= d_Y
		  if d_X > 1
        [r, ev(LV)] = eigs(S*S',1,'LA',opts);
      else
        r = 1;
        ev(LV) = S*S';
      end
		  
      t = X*r;
		  norm_t = norm(t);
		  r = r/norm_t;
		  t = t/norm_t;
		  p = X'*t;

		  c = Y'*t;      %c = S'*r;
		  bin = norm(c); %bin = sqrt(ev)/norm_t:
		  q = c/bin;
		  u = Y*q;

    else
      if d_Y > 1
        [q, ev(LV)] = eigs(S'*S,1,'LA',opts);
      else
        q = 1;
        ev(LV) = S'*S;
      end
      r = S*q;
		  t = X*r;
		  norm_t = norm(t);
		  r = r/norm_t;
		  t = t/norm_t;
		  p = X'*t;

			u = Y*q;
		  bin = u'*t;    %bin = ev/norm_t:
		  c = q*bin;
	  end

	  if LV == 1
		  if ev(LV) == 0
        error('PLS: Rank of the covariation matrix X''*Y is zero.');
      else
        v = p;
      end  
	  elseif ev(LV) <= 1e-16*ev(1)
      nLV = LV-1;
      WarnMsg = sprintf(['\nPLS: Rank of the covariation matrix X''*Y is exausted ' ...
                         'after removing %d Latent Variable(s).\n'...
                         'Results only for the %d LV(s) will be returned'],nLV,nLV);
      if exist('prwarning')
        prwarning(1, WarnMsg);
      else
        warning(WarnMsg);
      end
      break;
    else
  	  v = p - XRes.V(:,1:LV-1)*(XRes.V(:,1:LV-1)'*p);
	  end

	  v = v/norm(v);

    B(:,:,LV) = r*c';	  

    XRes.ssq(:,LV) = (p'*p)/ssq_X;
	  XRes.T(:,LV)   = t;
	  XRes.R(:,LV)   = r;
	  XRes.P(:,LV)   = p;
	  XRes.V(:,LV)   = v;

	  YRes.ssq(:,LV) = (bin.^2)/ssq_Y;
	  YRes.U(:,LV)   = u;
	  YRes.Q(:,LV)   = q;
	  YRes.C(:,LV)   = c;
	  YRes.bin(:,LV) = bin;
	  
	  S = S - v*(v'*S);
  end

  if nLV < Options.maxLV
    XRes.V = XRes.V(:,1:nLV);
  end
end

if nLV < Options.maxLV
  B = B(:,:,1:nLV);

  XRes.ssq = XRes.ssq(:,1:nLV);
  XRes.T   = XRes.T(:,1:nLV);  
  XRes.R   = XRes.R(:,1:nLV);    
  XRes.P   = XRes.P(:,1:nLV);  

  YRes.ssq = YRes.ssq(:,1:nLV);
  YRes.U   = YRes.U(:,1:nLV);  
  YRes.Q   = YRes.Q(:,1:nLV);  
  YRes.C   = YRes.C(:,1:nLV);  
  YRes.bin = YRes.bin(:,1:nLV);
 
  Options.maxLV = nLV;
end

B = cumsum(B,3);

return;



function A = MakeSym(A)
A = 0.5*(A+A');
%A = max(A,A');

return

