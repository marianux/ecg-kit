%GENCLASS Generate class frequency distribution
%
%  M = GENCLASS(N,P)
%
% INPUT
%  N    Number (scalar)
%  P    Prior probabilities
%
% OUTPUT
%  M    Class frequency distribution
%
% DESCRIPTION
% Generates a class frequency distribution M of N (scalar) samples
% over a set of classes with prior probabilities given by the vector P.
% The numbers of elements in P determines the number of classes and
% thereby the number of elements in M. P should be such that SUM(P) = 1. 
% If N is a vector with length C, then M=N is returned. This transparent
% behavior is implemented to avoid tests in other routines.
%
% Note that this is a random process, so M = GENCLASS(100,[0.5, 0.5]) 
% may result in M = [45 55].
%
% This routines is used in various data generation routines like
% GENDATH to determine the distribution of the objects over the classes.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: genclass.m,v 1.3 2006/09/26 12:55:32 duin Exp $

function N = genclass(N,p)
	
		
	if nargin < 2 | isempty(p)
		p = ones(1,length(N))/length(N);
	end
	c = length(p);
	if length(N) == c 
		;
	elseif length(N) > 1
		error('Mismatch in numbers of classes')
	else
		if nargin < 2 | isempty(p)
			p = repmat(1/c,1,c);
		end
		P = cumsum(p(:)');
		if abs(P(c)-1) > 1e-10
			error('Sum of class prior probabilities should be one')
    end
		X = rand(N,1);
    P = [0 P];
    Z = zeros(1,c);
    for j=1:c
      Z(j) = sum((X > P(j)) & (X <= P(j+1)));
    end
    N = Z;
% 		K = repmat(X,1,c) < repmat(P(:)',N,1);
% 		L = sum(K,1);
% 		N = L - [0 L(1:c-1)];
	end

	return
	
