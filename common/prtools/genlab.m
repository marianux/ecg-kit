%GENLAB Generate labels for classes
% 
%   LABELS = GENLAB(N,LABLIST)
% 
% INPUT
%   N        Number of labels to be generated
%   LABLIST  Label names (optional; default: numeric labels 1,2,3,...)
%
% OUTPUT
%   LABELS   Labels in a column vector or strinag array
%
% DESCRIPTION  
% Generate a set of labels as defined by the LABLIST. N is a number or 
% a row/column vector of the values for each class. If N is a vector, then
% the first N(i) labels get the value LABLIST(i,:). N should have as many 
% components as LABLIST. If N is a scalar, then N labels are generated for 
% each class. LABLIST is a column vector or a string array. Labels can be 
% used to construct a labeled dataset.
% 
% EXAMPLES
% Numeric labels, 1..10, 10 classes, 100 labels per class.
%    LAB1 = GENLAB(100*ones(1,10));          
% Character labels, 3 classes, 10 labels per class.
%    LAB2 = GENLAB([10 10 10], ['A';'B';'C']);
% Name labels, 2 classes, 50 labels per class. 
% The commands below are equivalent.
%    LAB3 = GENLAB([50 50], {'Apple'; 'Pear'}); 
%    LAB3 = GENLAB([50 50], ['Apple'; 'Pear ']);
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATASET

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: genlab.m,v 1.4 2009/09/23 08:33:55 duin Exp $

function labels = genlab(n,lablist)

		if (nargin == 1)
		% Create numeric labels.
		labels = [];
		lab = 1;
		if (all(n == n(1))) 		% All label categories have equal cardinalities.
			J = repmat([1:length(n)],n(1),1);
			labels = J(:);
		else
			for i = 1:length(n)
				labels = [labels; repmat(lab,n(i),1)];
				lab = lab+1;
			end
		end
    labels = int2str(labels);

	else    % LABLIST present

		% Create string or character labels.

		if iscell(lablist), lablist = lablist(:); end
		[m,ncol] = size(lablist);
		
		if m ~= length(n)
			% This is only possible when all label categories have equal cardinalities.
			if (length(n) > 1)
				error('Wrong number of labels supplied.')
			else
				J = repmat([1:m],n,1);	
				labels = lablist(J(:),:);
			end

		elseif (iscell(lablist))		% Cell array

			% We are here when e.g. GENLAB([10 10],{'Apple'; 'Pear'}) is called.
			labels = {};
			for i = 1:length(n)
				labels = [labels; repmat(lablist(i),n(i),1)];
			end
			labels = char(labels); % cell string labels are not supported anymore
			
		elseif (isstr(lablist)) 		% Character array

			% We are here when e.g. GENLAB([10 10],['Apple'; 'Pear ']) is called.
			labels = char(repmat(lablist(1,:),n(1),1));
			for i = 2:length(n)
				if (n(i) > 0)
					labels = char(labels,repmat(lablist(i,:),n(i),1));
				end
			end
			if (n(1) == 0)                  % First label category not wanted.
				labels(1,:) = [];
			end

		else						

			% Again numeric labels.
			% We are here, when e.g. GENLAB([10 10 10], [3;6;7]) is called.
			if (ncol > 1)
				error('Labels should be characters, strings or single numbers.')
			end
			labels = [];
			for i = 1:length(n)
				labels = [labels; repmat(lablist(i),n(i),1)];
			end
		end
	end

	return
