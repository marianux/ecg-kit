%IM_BPROPAGATION Fixed mapping or binary propagation (DIP_Image)
%
%   B = IM_BPROPAGATION(A1,A2,N,CONNECTIVITY,EDGE_CONDITION)
%
% INPUT
%   A1       Dataset with binary object images dataset (possibly multi-band)
%            to be treated as seed for the propagation
%   A2       Dataset with binary object images dataset (possibly multi-band)
%            to be treated as mask for the propagation
%   N        Number of iterations (default inf)
%   CONNECTIVITY    See BPROPAGATION
%   EDGE_CONDITION  Value of edge, default 1
%
% OUTPUT
%   B        Dataset with propagated images
%
% DESCRIPTION
% The binary images in A1 are dilated under the condition that the result 
% stays inside the components stored in A2.
%
% EXAMPLE
% a = delft_idb; a = seldat(a,9);        delfigs
% mask = a*im_gray*im_threshold;         figure, show(mask)
% seed = mask*im_berosion;               figure, show(seed)
% cleaned = im_bpropagation(seed,mask);  figure, show(cleaned)
% showfigs
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, BPROPAGATION

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_bpropagation(a1,a2,n,connect,edgecon)

		
	if nargin < 5 | isempty(edgecon), edgecon = 1; end
	if nargin < 4 | isempty(connect), connect = 2; end
	if nargin < 3 | isempty(n), n = inf; end

	if isdataset(a1)
		if ~isdataset(a2)
			error('Seed and mask images should be both datasets')
		end
		fsize = getfeatsize(a1);
		if any(getfeatsize(a2) ~= fsize)
			error('Image structures of seed and mask images should be identical')
		end
		if length(fsize) == 2, fsize = [fsize 1]; end
		if size(a1,1) ~= size(a2,1)
			error('Same number of seed and mask images expected')
		end
		out = [];
		seed = data2im(a1);
		mask = data2im(a2);
		for i=1:size(a1,1)
			for j=1:fsize(3)
				f = feval(mfilename,seed(:,:,j,i),mask(:,:,j,i),n,connect,edgecon);
				seed(:,:,j,i) = f;
			end
		end
		b = setdat(a,im2obj);
	elseif isdatafile(a1)
		if ~isdatafile(a2)
			error('Seed and mask images should be both datafiles')
		end
		b = dyadic(a1,mfilename,a2,{n,connect,edgecon});
  elseif isa(a1,'double') | isa(a2,'dip_image') % here we have a single image
		a1 = dip_image(a1,'bin');
		a2 = dip_image(a2,'bin');
    b = bpropagation(a1,a2,n,connect,edgecon);
	end
	
return
