%IM_CENTER Shift all binary images in dataset: center to center of gravity
%
%   B = IM_CENTER(A)
%   B = A*IM_CENTER
%
% The objects in the binary images are shifted such that their centers of 
% gravities are in the image center.
%
%   B = IM_CENTER(A,N)
%
% In all directions N rows and columns are added after shifting.
%
% SEE ALSO
% DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_center(a,n)

		
	if nargin < 2 | isempty(n), n= 0; end
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{n});
    b = setname(b,'Image centering');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{n});
		b = setfeatsize(b,getfeatsize(a));
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		if isa(a,'dip_image'), a = double(a); end
		a = im_box(a);
		[ty,tx] = size(a);
		[sy,sx] = size(a);
		mxy = im_mean(a);
		mx = 2*round(mxy(1)*sx)-1-sx;
		my = 2*round(mxy(2)*sy)-1-sy;
		if mx < 0
			a = [zeros(sy,-mx) a];
		elseif mx > 0
			a = [a zeros(sy,mx)];
		end
		sx = sx + abs(mx);
		if my < 0
			a = [zeros(-my,sx); a];
		elseif my > 0
			a = [a; zeros(my,sx)];
		end
		sy = sy + abs(my);
		if n > 0
			b = zeros(sy+2*n,sx+2*n);
			b(n+1:n+sy,n+1:n+sx) = a;
		else
			b = a;
		end
		[ry,rx] = size(b);
	end

return