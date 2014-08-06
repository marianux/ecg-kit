%RENUMLAB Renumber labels
% 
%   [NLAB,LABLIST]        = RENUMLAB(LABELS)
%   [NLAB1,NLAB2,LABLIST] = RENUMLAB(LABELS1,LABELS2)
%
% INPUT
%   LABELS,LABELS1,LABELS2  Array of labels
%
% OUTPUT
%    NLAB,NLAB1,NLAB2       Vector of numeric labels
%    LABLIST                Unique labels
%
% DESCRIPTION 
% If a single array of labels LABELS is supplied, it is converted and
% renumbered to a vector of numeric labels NLAB. The conversion table
% LABLIST is such that LABELS = LABLIST(NLAB,:). When two arrays LABELS1
% and LABELS2 are given, they are combined into two numeric label vectors
% NLAB1 and NLAB2 with a common conversion table LABLIST.
%
% Note that numeric labels with value -inf or NaN and string labels CHAR(0)
% are interpreted as missing labels. Their entry in NLAB will be 0 and they
% will not have an entry in LABLIST.
%
% A special command is
%
%   NLAB = RENUMLAB(LABELS,LABLIST)
%
% which returns the indices of LABELS in a given LABLIST.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASET

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: renumlab.m,v 1.7 2009/01/04 22:18:57 duin Exp $

function [out1,out2,out3] = renumlab(slab1,slab2)

		
	out1 = []; out2 = []; out3 = [];

	% Obsolete call?
	if (isempty(slab1)) & (nargin == 1 | isempty(slab2))
		return
	end

	% Clean up the symbolic labels.
	slab1 = clean_lab(slab1);

	if (nargin == 1) | (isempty(slab2))

		% Find unique labels for LABLIST and indices into LABLIST for NLAB.  
		[lablist,dummy,nlab] = unique(slab1,'rows');

		% Remove "missing label", if present.
		[lablist,nlab] = remove_empty_label(lablist,nlab);

		% Set output arguments.
		out1 = nlab; out2 = lablist;

	elseif nargout > 1 | nargin == 1

		% Check whether SLAB1 and SLAB2 match in type.
		if (isstr(slab1) ~= isstr(slab2))
			error(['Label lists do not match.' ... 
				newline 'They should be both characters, strings or numbers'])
		end
		
		% Clean up SLAB2 and put the union of SLAB1 and SLAB2 into SLAB.
		slab2 = clean_lab(slab2);
		if isstr(slab1), slab = char(slab1,slab2); else, slab = [slab1;slab2]; end

		% Find unique labels for LABLIST and indices into LABLIST for NLAB.
		if (iscell(slab))
			[lablist,dummy,nlab] = unique(slab);
		else
			[lablist,dummy,nlab] = unique(slab,'rows');
		end;

		% Remove "missing label", if present.
		[lablist,nlab] = remove_empty_label(lablist,nlab);

		% Set output arguments.
		out1 = nlab(1:size(slab1,1));      % Renumbered SLAB1 labels 
		out2 = nlab(size(slab1,1)+1:end);  % Renumbered SLAB2 labels
		out3 = lablist;
		
	else % nargout == 1 & nargin > 1, call like nlab = renumlab(labels,lablist)
		
		[k2,s2] = feval(mfilename,slab2);
		[n1,n2,s12] = feval(mfilename,slab1,slab2);
		
		% This gives a headache, but it seems to work
		R = zeros(size(s12,1),1);    % zeros for all possible labels
		J = find(n2~=0);
		R(n2(J)) = J;                % substitute the existing ones
		J = find(n1~=0);
		out1 = zeros(length(n1),1);  % space for output
		out1(J) = R(n1(J));          % replace existing ones by their index in lablist
			                             % pffft !!!!

		if 0
		if size(s12,1) > size(s2,1) % some labels are not in lablist
			% This gives a headache, but it seems to work
			S = zeros(max([n1;n2]),1);  % zeros for all possible labels
			R = S;
			R(n2) = [1:length(n2)];     % substitute the existing ones
			S(n2) = n2;                 % all indices zeros, except the existing ones
			n1 = S(n1);                 % replace non-existing ones by zeros
			J = find(n1~=0);            % these are the existing ones
			out1 = zeros(length(n1),1); % space for output
			out1(J) = R(n1(J));         % replace existing ones by their index in lablist
			                            % pffft !!!!                              
		else
			% easy!
			J = find(~isnan(n2));
			[nn2,listn2] = sort(n2(J));
			out1 = zeros(length(n2),1);
			out1(J) = listn2(n1(J));
      
			%alternative version
			%what is above seems not always working
			%J = find(n2~=0);
			%[nn2,listn2] = sort(n2(J));
			%listn2 = J(listn2);

			%J = find(n1~=0);
			%out1 = zeros(length(n1),1);      
			%out1(J) = listn2(n1(J));
		end
		end
		
	end

	return

% LAB = CLEAN_LAB(LAB)
%
% Clean labels; for now, just replaces occurences of NaN in numeric labels
% by -inf (both denoting "missing labels").

function slab = clean_lab(slab)

	if (iscell(slab))        % Cell array of labels.
		slab = char(slab);
	elseif isempty(slab)
		;
	elseif (size(slab,2) == 1) & (~isstr(slab))  	% Vector of numeric labels.
		J = isnan(slab);
		slab(J) = -inf;
	elseif (isstr(slab))     % Array of string labels.
		;
	else
		error('labels should be strings or scalars')
	end

	return

% [LABLIST,NLAB] = REMOVE_EMPTY_LABEL (LABLIST,NLAB)
%
% Removes the empty first label from LABLIST and NLAB (corresponding to the 
% "missing label"), if present.

function [lablist,nlab] = remove_empty_label (lablist,nlab)

% Find occurences of '' (in cell array label lists), '\0' (in string 
% label lists) or -inf (in numeric label lists) and remove them.

	if (iscellstr(lablist))    % Cell array of labels.	
		if (strcmp(lablist{1},'')), lablist(1) = []; nlab = nlab -1; end
	elseif (isstr(lablist))    % Array of string labels.
		if (strcmp(lablist(1,:),char(0))) | isempty(deblank(lablist(1,:)))
			lablist(1,:) = []; nlab = nlab -1; 
		end
	else
		% Vector of numeric labels.
		if (lablist(1) == -inf), lablist(1) = []; nlab = nlab -1; end
	end

	return
