%MDS_STRESS - Sammon stress between dissimilarity matrices
%
% 	E = MDS_STRESS(Q,DS,D)
%
% INPUT
% 	Q					Indicator of the Sammon stress; Q = -2,-1,0,1,2
% 	DS				Original distance matrix
% 	D 				Approximated distance matrix
%
% OUTPUT
% 	E 				Sammon stress
%
% DESCRIPTION
% Computes the Sammon stress between the original distance matrix Ds
% and the approximated distance matrix D, expressed as follows:
%
%  E = 1/(sum_{i<j} DS_{ij}^(q+2)) sum_{i<j} (DS_{ij} - D_{ij})^2 * DS_{ij}^q
%

%
% Copyright: Elzbieta Pekalska, Robert P.W. Duin, ela@ph.tn.tudelft.nl, 2000-2003
% Faculty of Applied Sciences, Delft University of Technology
%

function [e,alpha] = mds_stress (q,Ds,D,isratio)

	if nargin < 4
		isratio = 0;
	end

	[m,k]   = size(Ds);
	if any(size(D) ~= size(Ds)), 
		error ('The sizes of matrices do not match.');
	end
	mk = m*k;

	D  = +D;
	Ds = +Ds;

	% I is  the index of non-zero (> eps) values to be included 
	% for the computation of the stress

	I = 1:mk; 
	nanindex = find(isnan(Ds(:)) | isnan(D(:)));
	if ~isempty(nanindex),
		I(nanindex) = [];
	end
	O = [];
	if m == k & (length(intersect(find(D(:) < eps), 1:m+1:(mk))) == m),
		O  = 1:m+1:mk;
		Ds(O) = 1;     
		D (O) = 1;
    mm = m - 1;
	else
		mm = k;
	end

  if isratio,
    II = setdiff(I,O);
    alpha = sum((Ds(II).^q).*D(II).^2)/sum((Ds(II).^(q+1)).*D(II));
    Ds = alpha*Ds;
	else
		alpha = 1; 
	end
		
 
	c = sum(Ds(I).^(q+2)) - length(O);
	if q ~= 0,
		e = sum(Ds(I).^q .* ((Ds(I)-D(I)).^2))/c;
	else
		e = sum(((Ds(I)-D(I)).^2))/c;
	end
return; 

