function [ranset,seed]=randomset(tot,nel,seed)

%RANDOMSET draws randomly a subsample of nel cases out of tot.
%(It is called if not all (p+1)-subsets out of n will be considered.) 
%
% Required input arguments: 
%     tot : The total number of observations to consider
%     nel : The number of observations that the subsample must contain
% Optional input arguments:    
%     seed : To define the state of the generator (default=0)
%            (0 sets the generator to its default initial state)
%
% Output arguments:
%   ranset : Random subset of nel cases out of tot.
%   seed   : The corresponding state.
% 
%
% I/O:
%    [ranset,seed]=randomset(n,n/2,0)
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
    
if nargin==2
    seed=0;
end

for j=1:nel
   [random,seed]=uniran(seed);       
   num=floor(random*tot)+1;
   if j > 1
      while any(ranset==num)
         [random,seed]=uniran(seed);       
         num=floor(random*tot)+1;
      end   
   end
   ranset(j)=num;
end