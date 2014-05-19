%BHATM Bhattacharryya linear feature extraction mapping 
%
%	W = BHATM(A,N)
% 
% INPUT
% 	A		Dataset
% 	N		Number of dimensions to map to (N >= 1), or fraction of cumulative
%             contribution to retain (0 < N < 1)
%
% OUTPUT
% 	W       Bhattacharryya mapping
%
% DESCRIPTION  
% Finds a mapping of the labeled dataset A onto an N-dimensional linear
% subspace such that it maximizes the Bhattacharrryya distance between the
% classes, assuming Gaussian distributions. Only for two-class datasets.
%
% SEE ALSO
% MAPPINGS, DATASETS, FISHERM, NLFISHERM, KLM, PCA

% Copyright: F. van der Heijden(1) and Dick de Ridder(2)
% (1) Laboratory for Measurement and Instrumentation, Dept. of Electrical 
% Engineering, University of Twente, Enschede, The Netherlands 
% (2) Faculty of Electrical Engineering, Mathematics and Computer Science
% Delft University of Technology, Mekelweg 4, 2628 CD Delft, The Netherlands

% $Id: bhatm.m,v 1.5 2010/02/08 15:31:48 duin Exp $

function w = bhatm(a,n)

	
	if (nargin < 2)
		prwarning (4, 'fraction to map to not specified, assuming 0.9');
		n = 0.9; 
	end

	% If no arguments are specified, return an untrained mapping.

	if (nargin < 1) | (isempty(a))
		w = prmapping('bhatm',{n});
		w = setname(w,'Bhatta mapping');
		return
	end

    % Check input.
	a = testdatasize(a);
	islabtype(a,'crisp');
	isvaldset(a,2,2);               % At least 2 objects per class, 2 classes
	a = setprior(a,getprior(a)); % set priors to avoid unnecessary warnings
	[m,k,c] = getsize(a); 
   if (c ~= 2)
       error ('works only on two-class datasets');
   end;
    
	% Shift mean of data to origin.
	b = a*scalem(a); 

   % Get class means and covariances.
   [U,G] = meancov(b);
   G1 = G(:,:,1); G2 = G(:,:,2);		

   [E1,D1] = preig(G1); A = D1^-0.5*E1'; G2 = A*G2*A';   % Whiten the first class
   G2 = (G2+G2')/2;			                        % Enforce symmetry
   [E2,D2]=preig(G2);                                    % Decorrelate the second

   % Contributions to Bhattacharryya distance.
   d = E2'*A*(U(1,:)-U(2,:))'; D2 = diag(D2);
   Jb = 0.25*((d.^2)./(1.+D2)) + 0.5*log(0.5*(D2.^0.5 + D2.^-0.5));

   % Calculate transform matrix.
   W = E2'*D1^-0.5*E1';	   % The full transformation matrix
   [Jb,I] = sort(-Jb);        % Sort according to Jb 
   W = W(I,:);                % Same

	%DXD: If N=0 we just want to get the Jb for all dimensions
	if (n==0)
		w = cumsum(Jb)/sum(Jb);
		return
	end
	
   % If a fraction N is given, find number of dimensions.
   if (n < 1)
       Jb = cumsum(-Jb);		    % Accumulate Jb
       ind = find(Jb < n*Jb(k));	% Find important components
       if (length(ind)==0)
	       n = 1;
       else
	       n = ind(length(ind));	% D is last important component
       end
   end;
    
   % Construct affine mapping.
   rot = W(1:n,1:k)'; off = -mean(b*rot);
   w = affine(rot,off,a);
   w = setname(w,'Bhatta mapping');

return
