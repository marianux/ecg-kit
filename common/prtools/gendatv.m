%GENDATV Generation of a very large dataset
% 
% 	[A,U] = GENDATV(M,C,K,R,T)
% 	A = GENDATV(M,U)
% 
% INPUT
%   M   Number of objects to be generated (def. 10^4)
%   C   Number of classes  (def. 10^2)
%   K   Number of features (def 10)
%   R   std. dev of class distributions (def 1)
%   S   Std. dev. of class priors (def 0)
%   U   Dataset with class means, class priors and R
%
% OUTPUT
%   A   C-class dataset with M objects and K features
%   U   Dataset with class means, class priors and R
%
% DESCRIPTION
% This routine generates an arbitrarily large dataset of C classes by first
% generating the class means from a standard K-dimensional normal
% distribution. The class priors are genetated from a N(1,S) distribution,
% After that their sum is normalized to one. The classes themselves have
% spherical normal distributions with an equal standard deviation R for all
% features and all classes.
%
% In U the class means, the class priors and R are returned, such that
% similar data from the same distribution can be generated.
% 
% SEE ALSO
% DATASETS, PRDATASETS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function [a,u] = gendatv(varargin)

[m,c,k,r,s] = setdefaults(varargin,10^4,10^2,10,1,0);

if isdataset(c)
	u = c;
	r = getuser(u,'R');
  k = size(u,2);
else    
  u = gauss(c,zeros(1,k));
  u = setlabels(u,[1:c]');
  if s==0
    u = setprior(u,0);
  else
    p = 1+randn(1,c)*s;
    p(p<0) = 0;
    p = p/sum(p);
    u = setprior(u,p);
  end
  u = setuser(u,r,'R');
end
n = genclass(m,getprior(u));
a = gauss(n,u,r*eye(k));

return
