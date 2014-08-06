%IM_MEASURE Fixed mapping computating features by DIP_Image
%
%		F = IM_MEASURE(A,GRAY,FEATURES)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   GRAY     Gray-valued images (matched with A, optional)
%   FEATURES Features to be computed
%
% OUTPUT
%   F        Dataset with computed features
%
% In each image of the measurement set GRAY the features given in FEATURES 
% are measured. In A a segmented version of GRAY has to be supplied.
% When no GRAY is supplied, the binary images in A are used. Only
% the largest object in each image is considered.
%
% The following features may be computed:
%    'dimension','mean','stddev','gravity','size','center','max','min',
%    'maxval','minval','feret'','inertia','ccbendingenergy'.
% Note that some features like 'mean' (mean image intensity) and 'stddev'
% (standard deviation of image intensity) are not useful for binary images.
% Run MEASUREHELP to get some information on these measures.
%
% Use FEATURES = 'all' for computing all features.
% Use MEASUREHELP for some description of the features.
% Use IM_FEATURES for a set of features not based on DIP_Image
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, MEASURE, MEASUREHELP

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_measure(a,gray,features)

	 
	checktoolbox('diplib');
	if nargin < 3  features = []; end
	if nargin < 2  gray = []; end
	%if nargin < 2 | isempty(gray), gray = a; end

	if nargin < 1 | isempty(a)
		b = prmapping(mfilename,'fixed',{gray,features});
		b = setname(b,'DIP measurements');
	elseif isdataset(a)
		if (nargin < 3 | isempty(features)) & ~isdataset(gray)
			features = gray;
			gray = a;
		end
		if ~isdataset(gray)
			error('Binary and gray images should be both datasets')
		end
		fsize = getfeatsize(a);
		if any(getfeatsize(gray) ~= fsize)
			error('Image structures of binary and gray images should be identical')
		end
		if length(fsize) == 2, fsize = [fsize 1]; end
		if size(a,1) ~= size(gray,1)
			error('Same number of binary and gray images expected')
		end
		out = [];
		binim = data2im(a);
		grim = data2im(gray);
		nim = size(a,1)*fsize(3);
		s = sprintf('Measuring %i images',nim);
		prwaitbar(nim,s);
		for i=1:size(a,1)
			for j=1:fsize(3)
				prwaitbar(nim,(i-1)*fsize(3)+j);
				f = feval(mfilename,binim(:,:,j,i),grim(:,:,j,i),features);
				if isempty(out)
					out = repmat(f(:)',[size(a,1),1,fsize(3)]);
					%out = reshape(f(:)',[size(a,1),1,fsize(3)]);
				else
					out(i,:,j) = f;
				end
			end
		end
		prwaitbar(0);
		b = setdat(a,out);
		b = setfeatsize(b,[length(f),fsize(3)]);
		b = setfeatlab(b,getfeaturelabels(features));
	elseif isdatafile(a)
		if nargin < 3 | isempty(features) & ~isdatafile(gray)
			features = gray;
			gray = a;
		end
		if ~isdatafile(gray)
			error('Binary and gray images should be both datafiles')
		end
		b = dyadic(a,mfilename,gray,{features});
		%b = setfeatlab(b,getfeaturelabels(features));
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		if isempty(features), features = 'dimension'; end
		gray = 1.0*dip_image(gray);
		%labim = label(dip_image(im_select_blob(a),'bin'));
		labim = label(dip_image(a,'bin'));
    c = measure(labim,gray,'size',[],2);
    labid = c.id;
    sz = c.size;
    [bb,mm] = max(sz);
    labid = labid(mm);
		if strcmp(features,'all')
			features = {'dimension','mean','stddev','gravity',...
			'size','center','max','min', 'maxval','minval',...
			'feret','inertia', 'ccbendingenergy'};
		end
		b = measure(labim,gray,features,labid,2);
    b = double(b);
	else
		error('Wrong input')
	end
	
return

function names = getfeaturelabels(features)

names = {};
for i=1:length(features)
	switch features{i}
	case 'dimension'
		names{end+1} = 'imagewidth';
		names{end+1} = 'imageheight';
	case 'mean'
		names{end+1} = 'mean int';
	case {'mean', 'sum'}
		names{end+1} = 'mass';
	case 'stddev'
		names{end+1} = 'standard dev.';
	case 'gravity'
		names{end+1} = 'gravity x';
		names{end+1} = 'gravity y';
	case 'size'
		names{end+1} = 'size';
	case 'center'
		names{end+1} = 'center x';
		names{end+1} = 'center y';
	case 'max'
		names{end+1} = 'max x coord';
		names{end+1} = 'max y coord';
	case 'min'
		names{end+1} = 'min x coord';
		names{end+1} = 'min y coord';
	case 'maxval'
		names{end+1} = 'max int';
	case 'minval'
		names{end+1} = 'min int';
	case 'perimeter'
		names{end+1} = 'perimeter';
	case 'feret'
		names{end+1} = 'max diameter';
		names{end+1} = 'min diameter';
		names{end+1} = 'max perp. diameter';
	case 'inertia'
		names{end+1} = 'inertia moment 1';
		names{end+1} = 'inertia moment 2';
	case 'ccbendingenergy'
		names{end+1} = 'bending energy perimeter';
	otherwise
		error('I do not know feature %s.',features{i});
	end
end

return
