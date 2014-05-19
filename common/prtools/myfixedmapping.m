%MYFIXEDMAPPING Skeleton for a user supplied mapping
% 
%   W = MYFIXEDMAPPING([],PAR]
%   W = MYFIXEDMAPPING
%   B = A*MYFIXEDMAPPING
%   B = A*MYFIXEDMAPPING([],PAR]
%   B = MYFIXEDMAPPING(A,PAR)
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
% the routine user may program their own mapping. The routine as it is just
% selects the features (columns of A) as defined in PAR.
%
% SEE ALSO
% DATASETS, MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function OUT = myfixedmapping(A,PAR)
	% take care that there is always a second parameter, even if you don't
	% need it, PRTools will.

	% define here the default values for all parameters. Defaults might be
	% empty ([]), but the routine should then test on that if the parameter
	% is used.
	if (nargin < 2)
		PAR = [];
	end
	
	% determine the type of call
	if (nargin == 0) | (isempty(A))
		% we are here if the routine is called by W = myfixedmapping, or by 
    % W = myfixedmapping([],PAR)
		
		% this is a definition call. So we need to return a fixed mapping
		OUT = prmapping(mfilename,'fixed',{PAR});
		% the first parameter in the above call is the name of the routine that
		% should handle the mapping when it is called with data. Usually this
		% is this routine. Its name is returned by mfilename.
		% The third parameter is cell array containing all parameters of the
		% call. Defaults should be substituted.
		
		OUT = setname(OUT,'Fixed Mapping Skeleton');
		% Here just a name is supplied for your own information. It is
		% displayed when you call the routine without a ';'
		
	else
		% now we have a normal call with a dataset A and parameter PAR.
		% Let us first check the inputs
		
		isdataset(A);    % returns an error if A is not a dataset
		
		% checking parameter values depends on the mapping we want to
		% implement. In our fake example it should be in the range of feature
		% numbers. If there are still empty parameter values they should be
		% given a sensible value here.
		
		if isempty(PAR)  
			% select all features by defaults (i.e. do nothing)
			PAR = [1:size(A,2)];
		end
		
		if any(PAR < 1) | any(PAR > size(A,2))
			% Features to be selected should be in a proper range
			error('Feature number out of range')
		end
		
		% now the mapping itself should be programmed. We can simpy put
		% OUT = A(:,PAR);
		% but sometimes operations on the data stored in the dataset are needed
		% that cannot be done within PRTools. The following serves as an
		% example on how this might be done.
		data = getdata(A);  % isolate the data from the dataset
		data = data(:,PAR); % do whatever is needed. Replace this line by what
		                    % is needed for defining the desired operation of
												% myfixedmapping
		OUT = setdata(A,data); % resubstitute the resulting data into the
		                       % original data structure of A. The resulting
		                       % dataset has thereby the same name, user- and
		                       % ident-fields as input dataset.
	end

	return
