%IM_HARRIS Harris corner detector
%
%		X = IM_HARRIS(A,N,SIGMA)
%
% INPUT
%   A      Datafile or dataset with images
%   N      Number of desired Harris points per image (default 100)
%   SIGMA  Smoothing size (default 3)
%
% OUTPUT
%   X      Dataset with a [N,3] array with for every image 
%          x, y and strength per Harris point.
%
% DESCRIPTION
% We use Kosevi's [1] software to find the corner points according to
% Harris [2]. On top of the Kosevi Harris point detector we run
% - multi-feature images (e.g. color images) are averaged
% - only points that are maximum in a K x K window are selected. 
%   If less points are found, K is iteratively reduced.
%   The initial value of K is about 4*SIGMA.
% Although SIGMA can be interpreted as scaling parameter, it might be
% better to appropriately subsample images instead of using a large SIGMA.
%
% If you use this software for publications, please refer to [1] and [2].
%
% REFERENCES
% [1] P. D. Kovesi, MATLAB and Octave Functions for Computer Vision and 
%     Image Processing, School of Computer Science & Software Engineering,
%     The University of Western Australia.   Available from:
%     <http://www.csse.uwa.edu.au/~pk/research/matlabfns/>.
% [2] C. Harris and M. Stephens, A combined corner and edge detector,
%     Proc. 4th Alvey Vision Conf., 1988, pp. 147-151.
%
% EXAMPLE
% delfigs
% a = kimia;                  % take simple shapesas example
% b = gendat(a,25)*im_gray;   % just 25 images at random
% c = data2im(b);             % convert dataset to images for display
% x = im_harris(b,15,1);      % compute maximum 15 Harris points at scale 1
% y = data2im(x);             % unpack dataset with results 
% for j=1:25                  % show results one by one
%     figure(j); imagesc(c(:,:,1,j)); colormap gray; hold on
%     scatter(y(:,1,1,j),y(:,2,1,j),'r*');
% end
% showfigs

function x = im_harris(a,n,s,k)

	  dip_imagecheck;
	
	if nargin < 3 | isempty(s), s = 3; end
	if nargin < 2 | isempty(n), n = 100; end
	
  if nargin < 1 | isempty(a)
    x = prmapping(mfilename,'fixed',{n,s});
    x = setname(x,'Harris points');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
		a = im_gray(a);
    x = filtim(a,mfilename,{n,s});
  else
		a = double(a);
		if size(a,3) ~= 1
			error('2D gray value image expected');
		end
		[mr,mc] = size(a);         % here we put a border of 15 pixels around the image
		aa = bord(a,NaN,15);       % mirroring the border to avoid artificial corners
		cc = harris(aa,s);         % at the image border, finally we take the
		c = resize(cc,15,mr,mc);   % original image size (bord and resize are in ./private)
		
		cdip = 1.0*dip_image(c);
    
		fsizes = [25,21,17,13,11,9,7,5,3];
		F = find(fsizes < 4*s);
		fsizes = fsizes(F);
		for fsize = fsizes
			cmax = double(maxf(cdip,fsize)); % here we find the points that in a fsize-window
			J = find(cmax == c);             % are a local maximum
			[Iy,Ix] = ind2sub(size(c),J);    % second attempt to get rid of image border points
			K = find(Iy > 3 & Ix >3 & Iy < (mr-2) & Ix < (mc-2));
			J = J(K);
			Z = find(cmax(J) > 0);     % find non-zeros only
			J = J(Z);
			if length(J) > n           % if we have sufficient points
				break;                   % stop
			end                        % otherwise, get more local maxima
		end
		cs = cmax(J);               
		[ss,R] = sort(-cs);          % rank Harris points
		if length(J) >= n            % find x,y coordinates of first n Harris points
			[Iy,Ix] = ind2sub(size(c),J(R(1:n)));
			x = [Ix, Iy, -ss(1:n)];
		else 
			x = zeros(n,3);
			[Iy,Ix] = ind2sub(size(c),J(R));
			x(1:length(R),:) = [Ix, Iy, -ss(1:length(R))];
		end
	end
	
return

% HARRIS - Harris corner detector
%
% Usage:                 cim = harris(im, sigma)
%                [cim, r, c] = harris(im, sigma, thresh, radius, disp)
%  [cim, r, c, rsubp, csubp] = harris(im, sigma, thresh, radius, disp)
%
% Arguments:   
%            im     - image to be processed.
%            sigma  - standard deviation of smoothing Gaussian. Typical
%                     values to use might be 1-3.
%            thresh - threshold (optional). Try a value ~1000.
%            radius - radius of region considered in non-maximal
%                     suppression (optional). Typical values to use might
%                     be 1-3.
%            disp   - optional flag (0 or 1) indicating whether you want
%                     to display corners overlayed on the original
%                     image. This can be useful for parameter tuning. This
%                     defaults to 0
%
% Returns:
%            cim    - binary image marking corners.
%            r      - row coordinates of corner points.
%            c      - column coordinates of corner points.
%            rsubp  - If five return values are requested sub-pixel
%            csubp  - localization of feature points is attempted and
%                     returned as an additional set of floating point
%                     coords. Note that you may still want to use the integer
%                     valued coords to specify centres of correlation windows
%                     for feature matching.
%
% If thresh and radius are omitted from the argument list only 'cim' is returned
% as a raw corner strength image.  You may then want to look at the values
% within 'cim' to determine the appropriate threshold value to use. Note that
% the Harris corner strength varies with the intensity gradient raised to the
% 4th power.  Small changes in input image contrast result in huge changes in
% the appropriate threshold.

% References: 
% C.G. Harris and M.J. Stephens. "A combined corner and edge detector", 
% Proceedings Fourth Alvey Vision Conference, Manchester.
% pp 147-151, 1988.
%
% Alison Noble, "Descriptions of Image Surfaces", PhD thesis, Department
% of Engineering Science, Oxford University 1989, p45.

% Copyright (c) 2002-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in 
% all copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% March    2002 - original version
% December 2002 - updated comments
% August   2005 - changed so that code calls nonmaxsuppts

function [cim, r, c, rsubp, csubp] = harris(im, sigma, thresh, radius, disp)
    
    error(nargchk(2,5,nargin));
    if nargin == 4
	disp = 0;
    end
    
    if ~isa(im,'double')
	im = double(im);
    end

    subpixel = nargout == 5;
    
		d_x = [-1 0 1; -1 0 1; -1 0 1];   % Derivative masks
		d_y = d_x';
		Ix = conv2(im, d_x, 'same');      % Image derivatives
		Iy = conv2(im, d_y, 'same');    

    % Generate Gaussian filter of size 6*sigma (+/- 3sigma) and of
    % minimum size 1x1.
    g = fspecial('gaussian',max(1,fix(6*sigma)), sigma);
      
    Ix2 = conv2(Ix.^2, g, 'same'); % Smoothed squared image derivatives
    Iy2 = conv2(Iy.^2, g, 'same');
    Ixy = conv2(Ix.*Iy, g, 'same');

    % Compute the Harris corner measure. Note that there are two measures
    % that can be calculated.  I prefer the first one below as given by
    % Nobel in her thesis (reference above).  The second one (commented out)
    % requires setting a parameter, it is commonly suggested that k=0.04 - I
    % find this a bit arbitrary and unsatisfactory. 

    cim = (Ix2.*Iy2 - Ixy.^2)./(Ix2 + Iy2 + eps); % My preferred  measure.
%   k = 0.04;
%   cim = (Ix2.*Iy2 - Ixy.^2) - k*(Ix2 + Iy2).^2; % Original Harris measure.

    if nargin > 2   % We should perform nonmaximal suppression and threshold

	if disp  % Call nonmaxsuppts to so that image is displayed
	    if subpixel
		[r,c,rsubp,csubp] = nonmaxsuppts(cim, radius, thresh, im);
	    else
		[r,c] = nonmaxsuppts(cim, radius, thresh, im);		
	    end
	else     % Just do the nonmaximal suppression
	    if subpixel
		[r,c,rsubp,csubp] = nonmaxsuppts(cim, radius, thresh);
	    else
		[r,c] = nonmaxsuppts(cim, radius, thresh);		
	    end
	end
    end
    
% NONMAXSUPPTS - Non-maximal suppression for features/corners
%
% Non maxima suppression and thresholding for points generated by a feature
% or corner detector.
%
% Usage:   [r,c] = nonmaxsuppts(cim, radius, thresh, im)
%                                                    /
%                                                optional
%
%          [r,c, rsubp, csubp] = nonmaxsuppts(cim, radius, thresh, im)
%                                                             
% Arguments:
%            cim    - corner strength image.
%            radius - radius of region considered in non-maximal
%                     suppression. Typical values to use might
%                     be 1-3 pixels.
%            thresh - threshold.
%            im     - optional image data.  If this is supplied the
%                     thresholded corners are overlayed on this
%                     image. This can be useful for parameter tuning.
% Returns:
%            r      - row coordinates of corner points (integer valued).
%            c      - column coordinates of corner points.
%            rsubp  - If four return values are requested sub-pixel
%            csubp  - localization of feature points is attempted and
%                     returned as an additional set of floating point
%                     coords. Note that you may still want to use the integer
%                     valued coords to specify centres of correlation windows
%                     for feature matching.
%

% Copyright (c) 2003-2005 Peter Kovesi
% School of Computer Science & Software Engineering
% The University of Western Australia
% http://www.csse.uwa.edu.au/
% 
% Permission is hereby granted, free of charge, to any person obtaining a copy
% of this software and associated documentation files (the "Software"), to deal
% in the Software without restriction, subject to the following conditions:
% 
% The above copyright notice and this permission notice shall be included in all
% copies or substantial portions of the Software.
%
% The Software is provided "as is", without warranty of any kind.

% September 2003  Original version
% August    2005  Subpixel localization and Octave compatibility

function [r,c, rsubp, csubp] = nonmaxsuppts(cim, radius, thresh, im)

    v = version; Octave = v(1)<'5';     % Crude Octave test    
    subPixel = nargout == 4;            % We want sub-pixel locations    
    [rows,cols] = size(cim);
    
    % Extract local maxima by performing a grey scale morphological
    % dilation and then finding points in the corner strength image that
    % match the dilated image and are also greater than the threshold.
    
    sze = 2*radius+1;                   % Size of dilation mask.
    mx = ordfilt2(cim,sze^2,ones(sze)); % Grey-scale dilate.

    % Make mask to exclude points within radius of the image boundary. 
    bordermask = zeros(size(cim));
    bordermask(radius+1:end-radius, radius+1:end-radius) = 1;
    
    % Find maxima, threshold, and apply bordermask
    cimmx = (cim==mx) & (cim>thresh) & bordermask;
    
    [r,c] = find(cimmx);                % Find row,col coords.

    
    if subPixel        % Compute local maxima to sub pixel accuracy  
	if ~isempty(r) % ...if we have some ponts to work with
	
	ind = sub2ind(size(cim),r,c);   % 1D indices of feature points
	w = 1;         % Width that we look out on each side of the feature
                       % point to fit a local parabola
	
	% Indices of points above, below, left and right of feature point
	indrminus1 = max(ind-w,1);
	indrplus1  = min(ind+w,rows*cols);
	indcminus1 = max(ind-w*rows,1);
	indcplus1  = min(ind+w*rows,rows*cols);
	
	% Solve for quadratic down rows
	cy = cim(ind);
	ay = (cim(indrminus1) + cim(indrplus1))/2 - cy;
	by = ay + cy - cim(indrminus1);
	rowshift = -w*by./(2*ay);       % Maxima of quadradic
	
	% Solve for quadratic across columns	
	cx = cim(ind);
	ax = (cim(indcminus1) + cim(indcplus1))/2 - cx;
	bx = ax + cx - cim(indcminus1);    
	colshift = -w*bx./(2*ax);       % Maxima of quadradic

	rsubp = r+rowshift;  % Add subpixel corrections to original row
	csubp = c+colshift;  % and column coords.
	else
        rsubp = []; csubp = [];
	end
    end
    
    if nargin==4 & ~isempty(r)     % Overlay corners on supplied image.
	if Octave
	    warning('Only able to display points under Octave');
	    if subPixel
		plot(csubp,rsubp,'r+'), title('corners detected');
	    else
		plot(c,r,'r+'), title('corners detected');
	    end
	    [rows,cols] = size(cim);	    
	    axis([1 cols 1 rows]); axis('equal'); axis('ij')
	else
	    figure(1), imshow(im,[]), hold on
	    if subPixel
		plot(csubp,rsubp,'r+'), title('corners detected');
	    else	    
		plot(c,r,'r+'), title('corners detected');
	    end
	end
    end

