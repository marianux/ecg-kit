%SELECTIM Select one or more images in multiband image or dataset
%         depeciated, use BANDSEL
%
% B = SELECTIM(A,N)
% B = A*SELECTIM([],N)
%
% INPUT
%   A          Multiband image or dataset containing multiband images
%   N          Vector or scalar pointing to desired images
%
% OUTPUT
%   B          New, reduced, multiband image or dataset
%
% This function is obsolete. Use the more general command BANDSEL instead.

function b = selectim(a,n)

	if nargin < 2, n = 1; end
	if nargin < 1 | isempty(a)
		b = prmapping(mfilename,'fixed',{n})
		b = setname(b,'SelectImage');
		return
	end
	if isdataset(a)
		im = data2im(a);
		if ndims(im) < 3 | ndims(im) > 4
			error('3D images expected')
    elseif ndims(im) == 3
      im = im(:,:,n);
    else
      im = im(:,:,n,:);
    end
    fsize = size(im);
		if isobjim(a)
			b = setdat(a,im2obj(im,fsize(1:3)));
			b = setfeatsize(b,fsize(1:3));
		else
			b = im2feat(im);
		end
	elseif isdatafile(a)
		b = a*filtm([],'selectim',n);
	elseif isa(a,'double')
		if ndims(a) ~= 3
			error('3D images expected')
		end
		b = a(:,:,n);
	else
		error('Unexpected datatype')
	end
	
return
		
