%IM2OBJ Convert Matlab images or datafile to dataset object
%
%  B = IM2OBJ(IM,A)
%  B = IM2OBJ(IM,FEATSIZE)
%
% INPUT
%  IM        X*Y image, X*Y*C image, X*Y*K array of K images, 
%            X*Y*C*K array of color images, or cell-array of images
%            The images may be given as a datafile.
%  A         Input dataset
%  FEATSIZE  Vector with desired feature sizes to solve ambiguities
%
% OUTPUT
%  B   Dataset with IM added
%
% DESCRIPTION
% Add standard Matlab images, as objects, to an existing dataset A. If A is 
% not given, a new dataset is created. Images of type 'uint8' are converted
% to 'double' and divided by 256. The resulting feature size is X*Y or
% X*Y*C. The set of images IM may be given as a datafile.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM2FEAT, FEATIM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: im2obj.m,v 1.5 2008/07/28 08:59:47 duin Exp $

function a = im2obj (im,a)

			
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
			featsize = getfeatsize(a);
		else
			featsize = a;
			nodataset = 1;
		end
	else % no information on feature size, assume most simple solution
		nodataset = 1;
		if length(imsize) == 2
			featsize = imsize;
		else
			featsize = imsize(1:end-1);
		end
	end	
	if length(featsize) == length(imsize)
		nobj = 1;
	elseif length(featsize) == (length(imsize)-1)
		nobj = imsize(end);
	elseif length(featsize) == (length(imsize)-2) & imsize(3) == 1
		nobj = imsize(end);
	else
		wrongfeatsize;
	end
	
	if any(featsize ~= imsize(1:length(featsize)))
		wrongfeatsize;
	end

	if (isa(im,'cell'))
		
		% If IM is a cell array of images, unpack it and recursively call IM2OBJ
		% to add each image.

		im = im(:);																	% Reshape to 1D cell array.
		for i = 1:length(im)
			b = feval(mfilename,im{i});
			if ~isempty(a) & any(a.featsize ~= b.featsize)
				error('Images should have equal sizes')
			end
			a = [a; b];
		end
		
  elseif isdatafile(im)
    
		if (nargin < 2)
			a = prdataset([]);
		end
    testdatasize(im);
    a = [a; feval(mfilename,data2im(im))];
    
	else

		% If IM is an image or array of images, reshape it and add it in one go.

		% Convert to double, if necessary
		if (isa(im,'uint8'))
			im = double(im)/256; 
		else
			im = double(im);
		end

		% ready for the real work, at last!
		
		if nobj == 1
			im = im(:)';
		else
			im = shiftdim(im,ndims(im)-1);
			im = reshape(im,nobj,prod(featsize));
		end
		
		if nodataset
			a = prdataset(im);
			a = setfeatsize(a,featsize);
		else
			a = [a; im];
		end

	end	

return

function wrongfeatsize
	error('Desired feature size and size of supplied image array are inconsistent')
return
