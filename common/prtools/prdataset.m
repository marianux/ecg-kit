%PRDATASET Load and convert dataset from disk
%
%  A = PRDATASET(NAME,M,N)
%
% The dataset given in NAME is loaded from a .mat file and converted
% to the current PRTools definition. Objects and features requested
% by the index vectors M and N are returned.
%
% See PRDATA for loading arbitrary data into a PRTools dataset.
% See PRDATASETS for an overview of datasets.

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

% $Id: prdataset.m,v 1.3 2009/12/21 10:35:12 duin Exp $

function a = prdataset(name,M,N,no_message)

	prtrace(mfilename);
	
  if ~isstr(name) & prversion < 5
    % needed for using new prdataset routines with prtools4
    if nargin > 1
      a = dataset(name,M);
    else
      a = dataset(name);
    end
    return
  end
  
  
  persistent FIRST; 
	if isempty(FIRST), FIRST = 1; end
	if nargin < 4, no_message = 0; end
	if nargin < 3, N = []; end
	if nargin < 2, M = []; end
	if exist([name '.mat']) ~= 2
		error([newline '---- Dataset ''' name ''' not available ----'])
	end

	s = warning;
	warning off
	b = load(name);
	warning(s);
	% Try to find out what data we actually loaded:
	names = fieldnames(b);
	eval(['a = b.' names{1} ';']);

	if ~isdataset(a) & ndims(a) > 2
		% We loaded an image?  
		a = im2feat(a);
		prwarning(2,'Assumed that a feature image is loaded')
	else
		a = dataset(a);
	end
	[m,k] = size(a);
	% Maybe a subset of the features is requested:
	if ~isempty(N) 
		if min(size(N)) ~= 1
			error([newline 'Feature indices should be stored in vector'])
		end
		if max(N) > k
			error([newline 'Dataset has not the requested features'])
		end
		a = a(:,N); 
	end

	% Maybe a subset of the objects are requested:
	if ~isempty(M) 
		if min(size(M)) ~= 1
			error([newline 'Object indices should be stored in vector'])
		end
		if max(M) > m
			error([newline 'Dataset has not the requested objects'])
		end
		a = a(M,:); 
	end

	if FIRST & ~no_message
		disp(' ')
		disp('*** You are using one of the datasets distibuted by PRTools. Most of them')
		disp('*** are publicly available on the Internet and converted to the PRTools')
		disp('*** format. See the Contents of the datasets directory (HELP PRDATASETS)')
		disp('*** or the load routine of the specific dataset to retrieve its source.')
		disp(' ')
		FIRST = 0;
	end
	return
