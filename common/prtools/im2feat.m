%IM2FEAT Convert Matlab images or datafile to dataset feature
%
%   B = IM2FEAT(IM,A)
%   B = IM2OBJ(IM,OBJSIZE)
%
% INPUT
%   IM       X*Y image, X*Y*K array of K images, or cell-array of images.  
%            The images may be given as a datafile.
%   A        Input dataset
%   OBJSIZE  Vector with desired object sizes to solve ambiguities
%
% OUTPUT
%   B        Dataset with IM added
%
% DESCRIPTION
% Add standard Matlab images, as features, to an existing dataset A. If A is
% not given, a new dataset is created. Images of type 'uint8' are converted
% to 'double' and divided by 256. The set of images IM may be given as a 
% datafile.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM2OBJ, FEATIM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: im2feat.m,v 1.4 2008/07/28 08:59:47 duin Exp $

function a = im2feat (im,a)

		
	if iscell(im)
		imsize = size(im{1}); 
	else
		imsize = size(im);   
	end
	
	nodataset = 0;
	if (nargin > 1)
		if (~isdataset(a)) & size(a,1) ~= 1
			error('second argument should be a dataset or a size vector'); 
		end
		if isdataset(a)
			objsize = getobjsize(a);
		else
			objsize = a;
			nodataset = 1;
		end
	else % no information on feature size, assume most simple solution
		nodataset = 1;
		if length(imsize) == 2
			objsize = imsize;
		else
			objsize = imsize(1:end-1);
		end
	end	
	if length(objsize) == length(imsize)
		nfeat = 1;
	elseif length(objsize) == (length(imsize)-1)
		nfeat = imsize(end);
	elseif length(objsize) == (length(imsize)-2) & imsize(3) == 1
		nfeat = imsize(end);
	else
		wrongobjsize;
	end
	
	if any(objsize ~= imsize(1:length(objsize)))
		wrongobjsize;
	end


	if (isa(im,'cell'))

		% If IM is a cell array of images, unpack it and recursively call IM2FEAT
		% to add each image.

		im = im(:);
		for i = 1:length(im)
			b = feval(mfilename,im{i});
			if ~isempty(a) & any(a.objsize ~= b.objsize)
				error('Images should have equal sizes')
			end
			a = [a b];
		end
		    
  elseif isdatafile(im)
    
		if (nargin < 2)
			a = prdataset([]);
		end
    testdatasize(im);
    a = [a feval(mfilename,data2im(im))];
      
	else

		% If IM is an image or array of images, reshape it and add it in one go.

		% Convert to double, if necessary
		if (isa(im,'uint8'))
			im = double(im)/256; 
		else
			im = double(im);
		end

		% ready for the real work, at last!
		
		if nfeat == 1
			im = im(:)';
		else
			im = shiftdim(im,ndims(im)-1);
			im = reshape(im,nfeat,prod(objsize));
		end
		
		if nodataset
			a = prdataset(im');
			a = setobjsize(a,objsize);
		else
			a = [a; im'];
		end

	end	


return


function wrongobjsize
	error('Desired object size and size of supplied image array are inconsistent')
return
