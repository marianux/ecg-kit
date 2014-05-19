 %STAMP_MAP Stamping and storing of mappings for fast reusage
%
%		C = STAMP_MAP(A,W) or C = A*W
%		U = STAMP_MAP(V,W) or U = V*W
%
% This routine is equivalent to MAP except that it stores previous
% results (C or U) and retrieves them when the same inputs (A and W
% or V and W) are given. This is especially good in case of combining
% classifier studies when the same mappings have to be computed multiple
% times.
%
% Note that only the highest level of calls to MAP are stored by STAMP_MAP.
% If multiple calls to lower level calls (e.g. somewhere included in a
% complicated combined classifier) are expected, they should be called
% first explicitely.
% 
% STAMP_MAP is automatically called inside MAP (and thereby by every
% overloaded * operation between datasets and mappings) controlled by
% the following settings:
%
% STAMP_MAP(0)  - disables all storage and retrieval of mappings.
%                 This should be called when no duplicates are expected
%                 as it prevents the unnecessary checking of inputs.
% STAMP_MAP(1)  - Retrieving of stored results only only.
%                 No new results are stored.
% STAMP_MAP(2)  - Enables retrieval as well as storage of results.
%
% SEE ALSO
% DATASETS, MAPPINGS, MAP 

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function c = stamp_map(a,b)

	persistent CHECK_STAMPS nstamps cell_stamp stamps
	
	if isempty(CHECK_STAMPS)
		CHECK_STAMPS = 0;
		nstamps = 0;
		stamps = zeros(1000,3);
		cell_stamp = cell(1,1000);
	end
	
	if nargin == 0
		c = CHECK_STAMPS;
	elseif nargin == 1
		if a < 3
			CHECK_STAMPS = a;
			if nargout > 0 
				c = a;
			end
		else
			disp(nstamps)
			for j=1:nstamps
				fprintf('STAMPS: %10.6e %10.6e: ',stamps(j,:));
				cell_stamp{j}
			end
		end
		
	else 
			
		if CHECK_STAMPS == 0 % checking is off, no stamp_mapping
		
			c = []; % execution call: try to pass this routine
			
		elseif iscell(a) | iscell(b) % don't handle cell constructs
			
			c = [];
		
		elseif CHECK_STAMPS == 1 | CHECK_STAMPS == 2  % checking is on
			
			sa = getstamp(a); 
			sb = getstamp(b);
			Ja = find(stamps(1:nstamps,1)==sa);
			J = [];
			if ~isempty(Ja)
				Jb = find(stamps(Ja,2)==sb);
				if ~isempty(Jb)   % found, we have it!
					J = Ja(Jb);
					if ~isempty(J) & length(J) > 1
						error('Stamps are not unique')
					end 
					c = cell_stamp{J}; % use it!
					stamps(J,3) = stamps(J,3)+1;
					%disp(['found: ' getname(a) '..,..' getname(b)])
				end
			end
			if isempty(J)                       % not found, let map do the work
				checkstamps = CHECK_STAMPS;       % Take care that we skip stamp_map
				CHECK_STAMPS = checkstamps+2;     % when we call map
				c = prmap(a,b); % during execution of map CHECK_STAMPS will be reset
				CHECK_STAMPS = checkstamps;       % so this is in fact not needed
				if CHECK_STAMPS == 2        % store as well
					%disp(['stored: ' getname(a) '..,..' getname(b)])
					nstamps = nstamps + 1;
					if nstamps > length(cell_stamp)
						stamps = [stamps; zeros(1000,3)];
						cell_stamp = [cell_stamp cell(1,1000)];
					end 
					stamps(nstamps,1:2) = [sa,sb];
					cell_stamp{nstamps} = c; 
				end
			end
			
		elseif CHECK_STAMPS == 3 | CHECK_STAMPS == 4
			
			CHECK_STAMPS = CHECK_STAMPS-2;  % reset of computation of map
			c = [];
			
		else
			
			;
		
		end
		
	end
	
return 

%GETSTAMP Compute unique stamp for data
%
% S = GETSTAMP(A)
%
%Computes 64-bit stamp S of input A

function s = getstamp(a)

s = 0;
type = class(a);
	
switch(type)
	
	case 'char'
		a = abs(a);
		n = length(a(:));
		s = s + abs(sin(1:n))*(a(:));
		
	case 'double'
		n = length(a(:));
		s = s + abs(sin(1:n)*a(:));
		
	case {'prmapping','prdataset','prdatafile'}
		a = struct2cell(struct(a));
		s = s + getstamp(a);
		
	case 'struct'
		a = struct2cell(a);
		s = s + getstamp(a);
		
	case 'cell'
		for j=1:length(a(:))
			s = s + getstamp(a{j});
		end
		
	otherwise
		try
			n = length(a(:));
			s = s + abs(sin(1:n))*abs(double(a(:)));
		catch
			;
		end 
		
end