%IM_SKEL_MEASURE Computation by DIP_Image of skeleton-based features
%
%    F = IM_SKEL_MEASURE(A,FEATURES)
%
% INPUT
%   A        Dataset with binary object images dataset
%   FEATURES Features to be computed
%
% OUTPUT
%   F        Dataset with computed features
%
% DESCRIPTION
% The following features may be computed on the skeleton images in A:
% 'branch', 'end', 'link', 'single'. They should be combined in a cell
% array.
%
% Use FEATURES = 'all' for computing all features (default).
%
% SEE ALSO
% DATASETS, DATAFILES, IM_SKEL, DIP_IMAGE, BSKELETON

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_skel_meas(a,features)

		if nargin < 2 | isempty(features), features = 'all'; end
	if strcmp(features,'all')
		features = {'branch', 'end', 'link', 'single'};
	end

  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed');
    b = setname(b,'Skeleton features',{features});
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{features});
		b = setfeatlab(b,features);
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		b = [];
		if ~iscell(features), features = {features}; end
		for i = 1:length(features)
			if strcmp (features{i}, 'branch')
				b = [ b sum(getbranchpixel(a)) ];
			end;
			if strcmp (features{i}, 'end')  
				b = [ b sum(getendpixel(a)) ];
			end;
			if strcmp (features{i}, 'link')  
				b = [ b sum(getlinkpixel(a)) ];
			end;
			if strcmp (features{i}, 'single')
				b = [ b sum(getsinglepixel(a)) ];
			end;
		end;
	end
	
return
