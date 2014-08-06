%NORMM Fixed mapping for Minkowski-P distance normalization
% 
%   B = NORMM(A,P)
%   B = A*NORMM([],P)
%   B = A*NORMM(P)
% 
% INPUT
%   A  Dataset or matrix
%   P  Order of the Minkowski distance (optional; default: 1)
%
% OUTPUT
%   B  Dataset or matrix of normalized Minkowski-P distances 
%
% DESCRIPTION
% Normalizes the distances of all objects in the dataset A such that their
% Minkowski-P distances to the origin equal one. For P=1 (default), this is 
% useful for normalizing the probabilities. For 1-dimensional datasets 
% (SIZE(A,2)=1), a second feature is added before the normalization such 
% that A(:,2)=1-A(:,1).
% 
% Note that NORMM(P) or NORMM([],P) is a fixed mapping.
% 
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% MAPPINGS, DATASETS, CLASSC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = normm(varargin)
  
	mapname = 'Normalization Mapping';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    w = define_mapping(argin,'fixed',mapname);
    
  elseif mapping_task(argin,'training')			% Compute the mapping.
    [a,p] = deal(argin{:});
		[m,k] = size(a);
	
		if (k == 1) & isdataset(a)	
			
			% Normalisation of a 1-D dataset may only be used for one-class
			% classification outcomes. So we test whether this is the case as
			% good as possible and send out warnings as well
			
			a = [a 1-a];			% Add a second feature, handy for the normalization.
			k = 2;
			prwarning(4,'SIZE(A,2)=1, so the second feature 1-A(:,1) is added.');
			
			featlist = getfeatlab(a);
			if size(featlist,1) < 2
				error('No two class-names found; probably a wrong dataset used.')
			else
				a = setfeatlab(a,featlist(1:2,:));
			end
    end

		% Compute the Minkowski-P distances to the origin.
		if (p == 1)
			dist = sum(abs(a),2);					% Simplify the computations for P=1.
		else
			dist = sum(abs(a).^p,2).^(1/p);
		end

		% For non-zero distances, scale them to 1.
		J = find(dist ~= 0);
		w = a;
		if ~isempty(J)
			w(J,:) = +a(J,:)./repmat(dist(J,1),1,k);
    end
    
	end
	
return;
