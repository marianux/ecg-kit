function [idx,dm,mm,Ss] = kur_rce(x,mode)

%
% [idx,dm,mean,Cov] = kur_rce(x,mode)
%
% Outlier identification based on the study of univariate
% projections onto the directions maximizing and minimizing
% the kurtosis coefficient
%
% Inputs:   observations, x, (a matrix having a row for each observation)
%           mode, 0   use p maxizing and p minimizing directions
%                     iterate until no outliers are detected (default)
%                 >0  use p maxizing and p minimizing directions
%                     iterate "mode" times
%                 <0  use only maximizing directions
%                     iterate "mode" times
% Outputs:  idx, zero/one vector with ones in the positions of the
%                suspected outliers
%           dm,  robust Mahalanobis distances
%           mean, robust mean estimator
%           Ss,  robust covariance matrix estimator
%

% Daniel Pena / Francisco J Prieto 23/5/00

[n,p] = size(x);

%Quite esta restricción por obsoleta
% if (n > 750)|(p > 50),
%   disp('Data set too large for this code');
%   return
% end

if nargin < 2,
  mode = 0;
end
maxit = Inf;
if mode < 0,
  maxit = abs(mode);
  mode = -1;
elseif mode > 0,
  maxit = mode;
  mode = 1;
end

% Initialization of parameters

%% Cutoff points for projections

ctf1 = 3.05 + 0.324*p;
ctf2 = exp(1.522-0.268*log(p));
ctf3 = 3.5;
ctf4 = 3.0;

%% minimum numbers of observations

lmn1 = max(floor(0.5*(n+p+1)),2*p);
lmn2 = min(2*p,max(p+2,floor(0.25*n)));

% Removing suspect outliers

xx0 = [x (1:n)'];
xx = xx0;
nn = n;
it = 1;
while (it <= maxit),
  xx1 = xx(:,1:p);
  V = kur_nw(xx1,mode);
  pr = xx1*V;
  md = median(pr);
  prn = abs(pr - ones(nn,1)*md);
%   mad = median(prn);

  %Problema con el estimador robusto de varianza (mad) al tener muchos valores en 0 o haber distribuciones con poca varianza. 
  %Solución: Si alguna dimension es muy pequeña, tomo la media de entre los percentiles 50 y 75 de prn en esa
  %dimensión como estimacion del desvio std.
  prn_percentiles = prctile(prn, [50 95]);
  
  if( size(prn,2) == 1 )
      prn_percentiles = prn_percentiles';
  end
  
  mad = prn_percentiles(1,:);
  
  disp_ratio = prn_percentiles(2,:) ./ prn_percentiles(1,:);
  
  if( any(disp_ratio > 5) )
      % Para una distribucion normal, en 3 sigma encontramos el 99% de los
      % puntos, por este motivos pusimos el umbral en 5 para ver cuando la
      % estimacion es muy mala.
    dim_with_low_var = find(mad == 0);
    mad(dim_with_low_var) = mean( prctile( prn(:,dim_with_low_var), [50 75] ) );

    if( any(mad(dim_with_low_var) == 0) )
        
%         disp(['No se pudo estimar la varianza en las dimensiones ' num2str(dim_with_low_var) '. Posiblemente esos features estén mal calculados.' ]);
        
        %No se indica ningun outlier en este caso.
        idx = false(n,1);
        
        return
        
    end
    
  end
  
  tt = prn./(ones(nn,1)*mad);

  if p > 2,
    ctfo = ctf1 - (0:(p-2))*0.75*(ctf1-ctf2)/(p-2);
  else
    ctfo = ctf1;
  end

  if mode >= 0,
    tctf = ones(nn,1)*[ ctfo ctf2 ctf3*ones(1,p-1) ctf4 ];
  else
    if( p > 1)
        tctf = ones(nn,1)*[ ctfo ctf2 ];
    else
        tctf = ones(nn,1)*[ ctfo ];
    end
  end
  
  taux = tt./tctf;
  t = max(taux,[],2);

%   in = find(t > 1);
% Meti una ampliación del nivel de rechazo para outliers ya que parecía muy
% sensible.
  in = find(t > 5);
  nn = length(in);
  if nn == 0,
    break
  end
%   inn = find(t <= 1);
% Meti una ampliación del nivel de rechazo para outliers ya que parecía muy
% sensible.
  inn = find(t <= 5);
  nu = length(inn);

  if nu < lmn2,
    [v,ix] = sort(t);
    ixx = ix(1:lmn2);
    xx = xx(ixx,:);
    break
  end
  xx = xx(inn,:);
  [nn,pp] = size(xx);
  if nn <= lmn2,
    break
  end
  it = it + 1;
end

%% Extracting the indices of the outliers

idx = ones(n,1);
[nn,pp] = size(xx);
idx(xx(:,pp)) = zeros(nn,1);

% recheck observations and relabel them

%% Cutoffs for Mahalanobis distances, based on the Hotelling-t2
%% Fixed at 0.99 . To modifay change the value of  ctflvl

ctflvl = 0.99;
ctf = (p*(n-1)/(n-p))*finv(ctflvl,p,n-p);

%% Mahalanobis distances using center and scale estimators
%% based on good observations

sidx = sum(idx);
s1 = find(idx);
s2 = find(idx == 0);
if sidx > 0,
  xx1 = xx0(s1,:);
  xx1r = xx1(:,1:p);
end
xx2r = xx0(s2,1:p);

mm = mean(xx2r);
Ss = cov(xx2r);
if sidx > 0,
  aux1 = xx1r - ones(sidx,1)*mm;
  %Quito un poco del ruido que genera esta inversion de matriz
  warning off
  dd = Ss\aux1';
  warning on
  dms = sum((aux1.*dd')');
end

%% Ensure that at least lmn1 observations are considered good

ado = sidx + lmn1 - n;
if ado > 0,
  [dmss,idms] = sort(dms);
  idm1 = idms(1:ado);
  ido = xx0(s1,pp);
  s3 = ido(idm1);
  idx(s3) = zeros(ado,1);

  sidx = sum(idx);
  s1 = find(idx);
  s2 = find(idx == 0);
  if sidx > 0,
    xx1 = xx0(s1,:);
    xx1r = xx1(:,1:p);
  end
  xx2r = xx0(s2,1:p);

  mm = mean(xx2r);
  Ss = cov(xx2r);
  if sidx > 0,
    aux1 = xx1r - ones(sidx,1)*mm;
    dd = Ss\aux1';
    dms = sum((aux1.*dd')');
  end
end

%% Check remaining observations and relabel if appropriate

while sidx > 0,

  s1 = find(dms <= ctf);
  s2 = length(s1);
  if s2 == 0,
    break
  end
  s3 = xx1(s1,pp);
  idx(s3) = zeros(s2,1);
  sidx = sum(idx);

  s1 = find(idx);
  s2 = find(idx == 0);
  if sidx > 0,
    xx1 = xx0(s1,:);
    xx1r = xx1(:,1:p);
  end
  xx2r = xx0(s2,1:p);

  mm = mean(xx2r);
  Ss = cov(xx2r);
  if sidx > 0,
    aux1 = xx1r - ones(sidx,1)*mm;
    dd = Ss\aux1';
    dms = sum((aux1.*dd')');
  end

end

% Values to be returned

idx = idx == 1;

aux1 = x - ones(n,1)*mm;
%Quito un poco del ruido que genera esta inversion de matriz
warning off
dd = Ss\aux1';
warning on
dms = sum((aux1.*dd')');

dm = sqrt(dms);
