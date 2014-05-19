%FEATSELV Varying feature selection
% 
% W = FEATSELV(A)
% W = A*FEATSELV
%
% Selects all features with a non-zero variance.
% Classifiers can be trained like A*(FEATSELV*LDC([],1E-3)) to make
% use of this feature selection
% 
% SEE ALSO
% MAPPINGS, DATASETS, FEATEVAL, FEATSELO, FEATSELB, FEATSELF,
% FEATSEL, FEATSELP, FEATSELM, FEATSELI

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function w = featselv(a)

	
  % If no arguments are supplied, return an untrained mapping.

	if (nargin == 0) | (isempty(a))
    
		w = prmapping('featselv');
    
  else

	  [m,k,c] = getsize(a); 
    featlist = getfeatlab(a);
    v = std(+a);
    J = find(v > 1e-6 & v < 1e6);
    if isempty(J)
      error('All objects are equal')
    end

	  % Return the mapping found.
		w = featsel(k,J);
		if ~isempty(featlist)
			w = setlabels(w,featlist(J,:));
		end
    
  end
		
  w = setname(w,'Varying FeatSel');

return
