%IM_FILL_NORM Fixed mapping for normalizing image for display puproses
%
%   B = IM_FILL_NORM(A,N,BACKGROUND)
%
%Low level routine for the DATAFILE/SHOW command to display non-square
%images of the datafile A, inside square of NxN pixels. Empty areas are
%filled with gray.
%Empty parts of images are given the value BACKGROUND (default: gray (0.5));

function b = im_fill_norm(a,n,background)

if isa(a,'prdataset')
	isobjim(a);
	if isdatafile(a) & ~isempty(getuser(a,'number_bands_differs'))
		a = band2obj(a);
		outsize = 1;
	else
		outsize = [n n getfeatsize(a,3)];
	end
% 	outsize = [n n getfeatsize(a,3)];
	b = filtim(a,mfilename,{n,background},outsize);
else
	a = double(a);
	[x,y,p] = size(a);
	mx = max(a(:));
	mn = min(a(:));
	%b = ones(n,n,p);
	if x == 1
		b = imresize(a,[1 n]);
	elseif x > y
		a = imresize(a,max(round(n*[x,y]/x),[1,1]),'bilinear');
		k = size(a,2);
		s = floor((n-k)/2)+1;
		b = background*ones(n,n,p); 
		b(:,s:s+k-1,:) = a;
	else
		a = imresize(a,max(round(n*[x,y]/y),[1,1]),'bilinear');
		k = size(a,1);
		s = floor((n-k)/2)+1;
		b = background*ones(n,n,p);
		b(s:s+k-1,:,:) = a;
	end
	b = (b - mn)/(mx-mn+eps);
end
