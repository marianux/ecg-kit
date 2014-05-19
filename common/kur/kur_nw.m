function Vv = kur_nw(x0,mode)

%
%    V = kur_nw(x,mode)
%
% Computation of directions that maximize and minimize
% the kurtosis coefficient of the projections
% The values of the projections are also computed
% To be used as a subroutine of kur_rce
%
% Inputs:  x, observations (by rows)
%          mode, (see file kur_rce)
% Outputs: V, directions maximizing/minimizing kurtosis coef.
%             maximization directions first
%

% Daniel Pena/Francisco J Prieto 23/5/00

if nargin < 2,
  mode = 0;
end

% Initializations

%% Parameters (tolerances)

maxit = 100;

tol = 1.0e-5;
tol1 = 1.0e-6;
beta = 1.0e-4;
rho0 = 0.1;

[n,p0] = size(x0);

%% Initialization of vectors

Vv = [];
ti = 1;

%% Standardize data

mm = mean(x0);
S = cov(x0);
x = x0 - ones(n,1)*mm;

Rr = chol(S);
x = ((Rr')\(x'))';

% Computing directions

%% Choice of minimization/maximization

cff = 1;
if mode < 0,
  cffmax = 2;
else
  cffmax = 3;
end

%% Main loop to compute 2p directions

while (cff < cffmax),

  xx = x;
  p = p0;
  pin = p - 1;
  M = eye(p0);

  for i = 1:pin,
    if cff == 1,
      a = max_kur(xx);
    else
      a = min_kur(xx);
    end
    la = length(a);
    za = zeros(la,1); za(1) = 1;
    w = a - za; nw = w'*a;
    if abs(nw) > eps,
      Q = eye(la) - w*w'/nw;
    else
      Q = eye(la);
    end

%% Compute projected values

    Vv = [ Vv (M*a) ];

    Qp = Q(:,2:p);
    M = M*Qp;
    ti = ti + 1;

%% Reduce dimension

    Tt = xx*Q;
    xx = Tt(:,(2:p));

    p = p - 1;

  end

%% Compute last projection

  Vv = [ Vv M ];
  ti = ti + 1;

%% Proceed to minimization

  cff = cff + 1;

end

% Undo standardization transformation

Vv = Rr\Vv;
uaux = sum(Vv.*Vv);
Vv = Vv*diag(1../sqrt(uaux));
