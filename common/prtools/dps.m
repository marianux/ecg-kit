% DPS Correntropy based, hierarchical density preserving data split
% 
%   R	      = DPS(A,LEVELS,CLASSWISE)
%   [R H]   = DPS(A,LEVELS,CLASSWISE)
%
% INPUT
%   A			    Input dataset
%   LEVELS		Number of split levels, default: 3
%   CLASSWISE	Use (1, default) or ignore (0) label information
%
% OUTPUT
%   R			Index array with rotation set with 2^LEVELS folds
%   H			Hierarchy of splits
%
% DESCRIPTION
% Density Preserving Sampling (DPS) divides the input dataset into a given
% number of folds (2^LEVELS) by maximizing the correntropy between the folds
% and can be used as an alternative for cross-validation. The procedure is
% deterministic, so unlike cross-validation it does not need to be repeated.
%
% REFERENCE
% M. Budka, B. Gabrys, Correntropy-based density-preserving data sampling
% as an alternative to standard cross-validation, IJCNN2010, 1-8
% http://www.budka.co.uk/
%
% SEE ALSO
% DATASETS, CROSSVAL
  
function [R H] = dps(A,levels,classwise)
	
 	if (nargin<3) || isempty(classwise), classwise = 1; end
 	if (nargin<2) || isempty(levels), levels = 3; end
	
	islabtype(A,'crisp');

	H = zeros(levels,size(A,1));

	idxs = cell(1,levels+1);
	idxs{1} = {1:size(A,1)};
	for i = 1:levels
		for j = 1:2^(i-1)
			t = helper(A(idxs{i}{j},:),classwise);
			idxs{i+1}{2*j-1} = idxs{i}{j}(t{1});
			idxs{i+1}{2*j} = idxs{i}{j}(t{2});
		end
		for j = 1:length(idxs{i+1})
			H(i,idxs{i+1}{j}) = j;
		end
	end
	
	R = H(end,:);

end

function idxs = helper(A,classwise)

	% if the classes aro too small to be divided further, switch to classless mode
	if classwise && all(classsizes(A)>2), c = getsize(A,3);
	else classwise = 0;	c = 1;
	end
	 
	siz = zeros(2,1);
	idx = cell(1,c);

	for i = 1:c
 		if classwise, [B BI] = seldat(A,i);					% select i-th class
		else B = A; BI = (1:size(A,1));
		end
		
		m = length(BI);
		mask = true(1,m);									% mask is used for counting remaining objects
		D = +distm(B) + diag(inf(m,1));						% working distance matrix
		Dorg = D;											% original distance matrix
		idx{i} = nan(2,ceil(m/2));
		
		for j = 1:floor(m/2)
			[mD,I] = min(D,[],1);							% \
			[mmD,J] = min(mD);								%   find two closest objects
			I = I(J(1)); J = J(1);							% /

			mask(I) = 0; mask(J) = 0;						% mark them as used
			
			% split the objects to maximally increase coverage of both subsets
			if (mean(Dorg(I,idx{i}(1,1:j-1))) + mean(Dorg(J,idx{i}(2,1:j-1))) < ...
				mean(Dorg(I,idx{i}(2,1:j-1))) + mean(Dorg(J,idx{i}(1,1:j-1))))
				idx{i}(1,j) = J;
				idx{i}(2,j) = I;
			else  
				idx{i}(1,j) = I;
				idx{i}(2,j) = J;
			end
			
			% remove used objects from the distance matrix
			D(I,:) = inf; D(:,I) = inf;
			D(J,:) = inf; D(:,J) = inf;
		end
		if isempty(j), j = 0; end							% in case the loop is not entered at all
		
		% odd number of points in class
		if sum(mask)>0
			I = find(mask);
			if siz(1)<siz(2)
				idx{i}(1,end) = I;
			elseif siz(1)>siz(2)
				idx{i}(2,end) = I;
			else
				if (mean(Dorg(I,idx{i}(1,1:j))) < mean(Dorg(I,idx{i}(2,1:j))))
					idx{i}(2,j+1) = I;
				else
					idx{i}(1,j+1) = I;
				end
			end
		end
		
		% convert indexes from class-specific to dataset-specific
		idx{i}(1,~isnan(idx{i}(1,:))) = BI(idx{i}(1,~isnan(idx{i}(1,:))));
		idx{i}(2,~isnan(idx{i}(2,:))) = BI(idx{i}(2,~isnan(idx{i}(2,:))));

		% update fold sizes
		siz(1) = siz(1) + sum(~isnan(idx{i}(1,:)));
		siz(2) = siz(2) + sum(~isnan(idx{i}(2,:)));

	end

	idx = cell2mat(idx);
	
	idxs = cell(1,2);
	idxs{1} = idx(1,~isnan(idx(1,:)));
	idxs{2} = idx(2,~isnan(idx(2,:)));
	
end
