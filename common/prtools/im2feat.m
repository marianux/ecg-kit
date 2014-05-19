%IM2FEAT Convert Matlab images or datafile to dataset feature
%
%   B = IM2FEAT(IM,A)
%
% INPUT
%   IM    X*Y image, X*Y*K array of K images, or cell-array of images
%         The images may be given as a datafile.
%   A     Input dataset
%
% OUTPUT
%   B     Dataset with IM added
%
% DESCRIPTION
% Add standard Matlab images, as features, to an existing dataset A. If A is
% not given, a new dataset is created. Images of type 'uint8' are converted
% to 'double' and divided by 256. The set of images IM may be given as a 
% datafile.
%
% SEE ALSO
% DATASETS, DATAFILES, IM2OBJ, DATA2IM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: im2feat.m,v 1.4 2008/07/28 08:59:47 duin Exp $

function a = im2feat (im,a)

		if (isa(im,'cell'))

		% If IM is a cell array of images, unpack it and recursively call IM2FEAT
		% to add each image.

		im = im(:);
		if (nargin < 2)
			a = prdataset([]);
		end
		n = length(im);
		s = sprintf('Converting %i images: ',n);
		prwaitbar(n,s);
		for i = 1:length(im)
			prwaitbar(n,i,[s int2str(i)]);
			a = [a feval(mfilename,im{i})];
		end
		prwaitbar(0);
    
  elseif isdatafile(im)
    
		if (nargin < 2)
			a = prdataset([]);
		end
    testdatasize(im);
    a = [a feval(mfilename,data2im(im))];
      
	else

		% If IM is an image or array of images, reshape it and add it in one go.

		n = size(im,3); imsize = size(im); 
    if length(imsize) == 4
      m = imsize(4);
    else
      m = 1;
    end
    imsize = imsize(1:2);

		% Convert to double, if necessary
		if (isa(im,'uint8') | isa(im,'uint16'))
			prwarning(4,'Image is uint; converting to double and dividing by 256');
			im = double(im)/256;
		end

		% Reshape images to vectors and store bands as features in dataset.
    if m==1
		  imm = reshape(im,imsize(1)*imsize(2),n);
      objsize = imsize(1:2);
		else
			rtc = imsize(1)*imsize(2);
			imm = zeros(m*rtc,n);
			for i=1:m
				imm(rtc*(i-1)+1:rtc*(i),:) = reshape(im(:,:,:,i),rtc,n);
			end
			objsize = [imsize(1:2) m]; 
% 			imm = [];
%       for i=1:m
%         imm = [imm; reshape(im(:,:,:,i),imsize(1)*imsize(2),n)];
%       end
%       objsize = [imsize(1:2) m];
    end
		if (nargin > 1)
			if (~isa(a,'prdataset'))
				error('Second argument is not a dataset')
			end
			if (size(imm,1) ~= size(a,1))
				error('Image size and dataset object size do not match')
			end
			a = [a imm];
		else
			a = prdataset(imm);
			a = setobjsize(a,objsize);
		end
	end

return
