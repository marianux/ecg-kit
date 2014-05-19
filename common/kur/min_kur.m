function [v,f] = min_kur(x,d0,maxit,show)

%
% [v,f] = min_kur(x,d,maxit,show)
%
% Compute direction minimizing the kurtosis of the projections
% Observations passed to the routine should be standardized
% Uses Newton's method and an augmented Lagrangian merit function
% To be used as a subroutine of kur_rce
%
% Inputs:  x, observations (by rows)
%          d, initial estimate of the direction (optional)
%          maxit, a limit to the number of Newton iterations
%          show, <>0 generates output for each iteration
% Outputs: v, min. kurtosis direction (local optimizer)
%          f, min. kurtosis value
%

% Daniel Pena/Francisco J Prieto 23/5/00

% Initialization

maxitdefault = 30;

if nargin < 4,
  show = 0;
end
if nargin < 3,
  maxit = maxitdefault;
end
if maxit <= 0,
  maxit = maxitdefault;
end
if nargin < 2,
  d0 = [];
end

%% Tolerances

maxitini = 1;
tol = 1.0e-4;
tol1 = 1.0e-7;
tol2 = 1.0e-2;
beta = 1.0e-4;
rho0 = 0.1;

[n,p] = size(x);

%% Initial estimate of the direction

if length(d0) == 0,

  uv = sum((x.*x)')';
  uw = 1../(eps + sqrt(uv));
  uu = x.*(uw*ones(1,p));

  Su = cov(uu);
  [V,D] = eig(Su);

  r = [];
  for i = 1:p,
    r = [ r val_kur(x,V(:,i)) ];
  end
  [v,ik] = min(r);
  a = V(:,ik(1));

  itini = 1;
  difa = 1;
  while (itini <= maxitini)&(difa > tol2),
    z = x*a;
    zaux = z.^2;
    xaux = x'.*(ones(p,1)*zaux');
    H = 12*xaux*x;
    [V,E] = eig(H);
    [vv,iv] = min(diag(E));
    aa = V(:,iv(1));
    difa = norm(a - aa);
    a = aa;
    itini = itini + 1;
  end
else
  a = d0/norm(d0);
end

%% Values at iteration 0 for the optimization algorithm

z = x*a;
sk = sum(z.^4);
lam = 2*sk;
f = sk;
g = (4*(z.^3)'*x)';
zaux = z.^2;
xaux = x'.*(ones(p,1)*zaux');
H = 12*xaux*x;

al = 0;
it = 0;
diff = 1;
rho = rho0;
clkr = 0;

c = 0;

% Newton method starting from the initial direction

if show,
  disp(' It.      F.obj.       | g |           c         alfa       rho');
end

while 1,

%% Check termination conditions

  gl = g - 2*lam*a;

  A = 2*a';
  [Q,W] = qr(a);
  Z = Q(:,(2:p));

  if show,
    aa = sprintf('%3.0f  %12.5f %13.4e',it,f,norm(gl));
    bb = sprintf(' %13.4e %8.3f %11.2e',abs(a'*a-1),al,rho);
    disp([ aa bb ]);
  end

  crit = norm(gl) + abs(c);
  if (crit <= tol)|(it >= maxit),
    break
  end

%% Compute search direction

  Hl = H - 2*lam*eye(p,p);
  Hr = Z'*Hl*Z;
  [V,E] = eig(Hr);
  Es = max(abs(E),1.0e-4);
  Hs = V*Es*V';

  py = - c/(A*A');
  rhs = Z'*(g + H*A'*py);
  pz = - Hs\rhs;
  pp = Z*pz + py*A';

  dlam = (2*a)\(gl + H*pp);

%% Adjust penalty parameter

  f0d = gl'*pp + 2*rho*c*a'*pp - dlam*c;
  crit1 = beta*norm(pp)^2;

  if f0d > -crit1,
    rho1 = 2*(crit1 + f0d)/(eps + c^2);
    rho = max([2*rho1 1.5*rho rho0]);
    f = sk - lam*c + 0.5*rho*c^2;
    f0d = gl'*pp + 2*rho*c*a'*pp - dlam*c;
    clkr = 0;
  elseif (f0d < -1000*crit1)&(rho > rho0),
    rho1 = 2*(crit1 - gl'*pp + dlam*c)/(eps + c^2);
    if (clkr == 4)&(rho > 2*rho1),
      rho = 0.5*rho;
      f = sk - lam*c + 0.5*rho*c^2;
      f0d = gl'*pp + 2*rho*c*a'*pp - dlam*c;
      clkr = 0;
    else
      clkr = clkr + 1;
    end
  end
  if (abs(f0d) < tol1),
    break
  end

%% Line search

  al = 1;
  itbl = 0;
  while itbl < 20,
    aa = a + al*pp;
    lama = lam + al*dlam;
    zz = x*aa;
    cc = aa'*aa - 1;
    sk = sum(zz.^4);
    ff = sk - lama*cc + 0.5*rho*cc^2;

    if ff < f + 0.0001*al*f0d,
      break
    end
    al = al/2;
    itbl = itbl + 1;
  end
  if itbl >= 20,
    if show,
      disp('Error in the line search');
    end
    break
  end

%% Update values for the next iteration

  a = aa;
  lam = lama;
  z = zz;

  nmd2 = a'*a;
  c = nmd2 - 1;
  f = sk - lam*c + 0.5*rho*c^2;
  g = (4*(z.^3)'*x)';
  zaux = z.^2;
  xaux = x'.*(ones(p,1)*zaux');
  H = 12*xaux*x;

  it = it + 1;

end

% Values to be returned

v = a/norm(a);
xa = x*v;
f = sum(xa.^4)/n;
