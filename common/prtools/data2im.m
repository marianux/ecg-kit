%DATA2IM Convert PRTools dataset or datafile to image
%
%   IM = DATA2IM(A,J)
%   IM = DATA2IM(A(J,:))
%
% INPUT
%   A     Dataset or datafile containing images
%   J     Desired images
%
% OUTPUT
%   IM    If A is dataset, IM is a X*Y*N*K matrix with K images.
%         - K is the number of images (length(J))
%         - N is the number of bands per image.
%         - N = 3 for RGB images, N = 1 for gray value images.
%         If A is a datafile, IM is a cell array of K images.
%
% DESCRIPTION
% An image, or a set of images stored in the objects or features of the
% dataset A are retrieved and returned as a 3D matrix IM. In case A is a
% datafile the images are stored in a cell array, except when a single
% image is requested.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, IM2OBJ, IM2FEAT

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: data2im.m,v 1.13 2010/02/03 13:17:17 duin Exp $

function im = data2im(a,J)

		
	if nargin > 1, a = a(J,:); end									 
	if isdatafile(a)
		m = size(a,1);
		im = cell(1,m);
		s = sprintf('Unpacking %i images: ',m);
		prwaitbar(m,s);
		for j=1:m
			prwaitbar(m,j,[s int2str(j)]);
			%im{j} = readdatafile(a,1,0);
			im{j} = feval(mfilename,prdataset(a(j,:)));
		end
		prwaitbar(0);
		if m==1
			im = im{1};
		end
		return
	end

	%a = testdatasize(a); % Oeps, datafiles are first converted to datasets
	                     % and then to images. This can be done better!
											 
	%isdataim(a);			% Assert that A contains image(s).

	data = +a;				% Extract data from dataset, for computational advantage.
	[m,k] = size(a); [objsize,featsize] = get(a,'objsize','featsize');

		
	% Reshape data into output array.

	if (isfeatim(a))	
		% A contains K images stored as features (each object is a pixel).
		if length(objsize) == 1
			im = zeros(1,objsize(1),k);
			for j = 1:k
				im(1,:,j) = reshape(data(:,j),1,objsize(1));
			end
		elseif length(objsize) == 2
			im = zeros(objsize(1),objsize(2),k);
			for j = 1:k
				im(:,:,j) = reshape(data(:,j),objsize(1),objsize(2));
			end
		elseif length(objsize) == 3
			im = zeros(objsize(1),objsize(2),k,objsize(3));
			for j = 1:k
				im(:,:,j,:) = reshape(data(:,j),objsize(1),objsize(2),objsize(3));
			end
		else
			error('Unable to handle these images')
		end
			
	else							
		
		% A contains M images stored as objects (each feature is a pixel).
		if length(featsize) == 1
			im = zeros(1,featsize(1),1,m);
			for j = 1:m
				im(1,:,1,j) = reshape(data(j,:),1,featsize(1));
			end
		elseif length(featsize) == 2
			im = zeros(featsize(1),featsize(2),1,m);
			for j = 1:m
				im(:,:,1,j) = reshape(data(j,:),featsize(1),featsize(2));
			end
		elseif length(featsize) == 3
			im = zeros(featsize(1),featsize(2),featsize(3),m);
			for j = 1:m
				im(:,:,:,j) = reshape(data(j,:),featsize(1),featsize(2),featsize(3));
      end
%     elseif length(featsize) == 4
%       if featsize(3) == 1
%         im = zeros(featsize(1),featsize(2),featsize(4),m);
%         for j = 1:m
%           im(:,:,1,j) = reshape(data(j,:),featsize(1),featsize(2),featsize(4));
%         end
%       elseif featsize(3) == 3
        
		else
			error('Unable to handle these images')
		end
		%im = squeeze(im); % some routines, like filtim, fail by squeezing 
	end

return
