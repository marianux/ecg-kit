%MYFIXEDMAPPING Skeleton for a user supplied fixed mapping
% 
%   W = MYFIXEDMAPPING([],PAR]
%   W = MYFIXEDMAPPING(PAR)
%
%   B = MYFIXEDMAPPING(A,PAR)
%   B = A*MYFIXEDMAPPING([],PAR]
%   B = A*MYFIXEDMAPPING(PAR)
% 
% INPUT
%   A      Dataset
%   PAR    Parameter
% 
% OUTPUT
%   W      Mapping definition
%   B      Dataset A mapped by MYFIXEDMAPPING
%
% DESCRIPTION
% This is a fake routine just offering the skeleton of a user supplied
% fixed mapping. By changing the lines indicated in the source and renaming
% the routine users may program their own mapping. The routine as it is
% just selects the features (columns of A) as defined in PAR.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function out = myfixedmapping(varargin)

  % give the name you like
  mapname = 'MyMapping';
  
  % the next statement is optional. It shifts the input parameters
  % by testing the type of the first one. Thereby it is possible to write
  % B = A*MYFIXEDMAPPING(PAR) instead of B = A*MYFIXEDMAPPING([],PAR)
	argin = shiftargin(varargin,'vector');
  
  % define the defaults for empty or not existing inputs
  argin = setdefaults(argin,[],[]);
  
  % find out from argin whether the mapping is defined 
  % (no dataset given) or executed.
  if mapping_task(argin,'definition')
    out = define_mapping(argin,'fixed',mapname);
    
  % Execute the mapping for a given dataset
  elseif mapping_task(argin,'training')			%  a mapping.
  
    % first retrieve the input parameters
    [a,par] = deal(argin{:});
		% now we have a normal call with a dataset A and parameter PAR.
    
		% Let us first check the inputs
		isdataset(a);    % returns an error if A is not a dataset
		
		% checking parameter values depends on the mapping we want to
		% implement. In our fake example it should be in the range of feature
		% numbers. 
		
		if isempty(par)  
			% select all features by defaults (i.e. do nothing)
			par = [1:size(a,2)];
		end
		
		if any(par < 1) | any(par > size(a,2))
			% Features to be selected should be in a proper range
			error('Feature number out of range')
		end
		
		% now the mapping itself should be programmed. We can simpy put
		% OUT = A(:,PAR);
		% but sometimes operations on the data stored in the dataset are needed
		% that cannot be done within PRTools. The following serves as an
		% example on how this might be done.
		data = getdata(a);  % isolate the data from the dataset
		data = data(:,par); % do whatever is needed. Replace this line by what
		                    % is needed for defining the desired operation of
												% myfixedmapping
		out = setdata(a,data); % resubstitute the resulting data into the
		                       % original data structure of A. The resulting
		                       % dataset has thereby the same name, user- and
		                       % ident-fields as input dataset.
	end

	return
