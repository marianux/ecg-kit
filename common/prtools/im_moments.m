%IM_MOMENTS PRTools routine for computing central moments of object images
%
%	  M = IM_MOMENTS(A,TYPE,MOMENTS)
%	  M = A*IM_MOMENTS([],TYPE,MOMENTS)
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%   TYPE     Desired type of moments
%   MOMENTS  Desired moments
%
% OUTPUT
%   M        Dataset with moments replacing images (poosibly multi-band)
%
% DESCRIPTION
% Computes for all images in A a (1*N) vector M moments as defined by TYPE
% and MOMENTS. The following types are supported:
%
% TYPE = 'none'     Standard moments as specified in the Nx2 array MOMENTS.
%                   Moments are computed with respect to the image center.
%                   This is the default for TYPE.
%                   Default MOMENTS = [1 0; 0 1];
% TYPE = 'central'  Central moments as specified in the Nx2 array MOMENTS.
%                   Moments are computed with respect to the image mean
%                   Default MOMENTS = [2 0; 1 1; 0 2], which computes
%                   the variance in the x-direction (horizontal), the
%                   covariance between x and y and the variance in the
%                   y-direction (vertical).
% TYPE = 'scaled'   Scale-invariant moments as specified in the Nx2 array
%                   MOMENTS. Default MOMENTS = [2 0; 1 1; 0 2].
%                   After: M. Sonka et al.,
%                   Image processing, analysis and machine vision.
% TYPE = 'hu'       Calculates 7 moments of Hu, invariant to translation,
%                   rotation and scale.
%                   After: M. Sonka et al.,
%                   Image processing, analysis and machine vision.
% TYPE = 'zer'      Calculates the Zernike moments up to the order as 
%                   specified in the scalar MOMENTS (1 <= MOMENTS <= 12). 
%                   MOMENTS = 12 generates in total 47 moments.
%                   After: A. Khotanzad and Y.H. Hong, Invariant image
%                   recognition by Zernike moments, IEEE-PAMI, vol. 12,
%                   no. 5, 1990, 489-497.
%
% SEE ALSO
% DATASETS, DATAFILES

% Copyright: D. de Ridder, R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_moments(a,type,mom)

	
	if nargin < 3, mom = []; end
	if nargin < 2 | isempty(type), type = 'none'; end
	
	
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{type,mom});
    b = setname(b,'Image moments');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtim(a,mfilename,{type,mom});
  elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
		if isa(a,'dip_image'), a = double(a); end
		switch type
		case {'none'}
			if isempty(mom)
				mom = [1 0; 0 1];
			end
			b = moments(a,mom(:,1),mom(:,2),0,0);
		case {'central'}
			if isempty(mom)
				mom = [2 0; 1 1; 0 2];
			end
			b = moments(a,mom(:,1),mom(:,2),1,0);		
		case {'scaled'}
			if isempty(mom)
				mom = [2 0; 1 1; 0 2];
			end
			b = moments(a,mom(:,1)',mom(:,2)',1,1);		
		case {'hu' 'Hu'}
			b = hu_moments(a);
		case {'zer' 'zernike' 'Zernike'}
			if isempty(mom)
				mom = 12;
			end
			b = zernike_moments(a,mom);
		otherwise
			error('Moments should be of type none, central, scaled, hu or zer')
		end
  else
    error('Illegal datatype for input')
  end
		
return
		
% M = MOMENTS (IM, P, Q, CENTRAL, SCALED)
%
% Calculates moments of order (P+Q) (can be arrays of indentical length)
% on image IM. If CENTRAL is set to 1 (default: 0), returns translation-
% invariant moments; if SCALED is set to 1 (default: 0), returns scale-
% invariant moments.
%
% After: M. Sonka et al., Image processing, analysis and machine vision.

function m = moments (im,p,q,central,scaled)

	if (nargin < 5), scaled = 0; 	end;
	if (nargin < 4), central = 0; end;
  if (nargin < 3)
  	error ('Insufficient number of parameters.');
  end;
   
	if (length(p) ~= length(q))
  	error ('Arrays P and Q should have equal length.');
  end;
   
  if (scaled & ~central)
  	error ('Scale-invariant moments should always be central.');
  end;

	% xx, yy are grids with co-ordinates
  [xs,ys] = size(im);
  [xx,yy] = meshgrid(-(ys-1)/2:1:(ys-1)/2,-(xs-1)/2:1:(xs-1)/2);
   
	if (central)
      
  	% Calculate zeroth and first order moments
	  m00 = sum(sum(im));
	  m10 = sum(sum(im.*xx));
	  m01 = sum(sum(im.*yy));
      
    % This gives the center of gravity
    xc  = m10/m00;
    yc  = m01/m00;
      
    % Subtract this from the grids to center the object
    xx  = xx - xc;
    yy  = yy - yc;
      
  end;
   
  % Calculate moment(s) (p,q).
  for i = 1:length(p)
		m(i) = sum(sum((xx.^p(i)).*(yy.^q(i)).*im));
  end;
   
  if (scaled)
      
  	c = 1 + (p+q)/2;
      
    % m00 should be known, as scaled moments are always central
    m = m ./ (m00.^c);
      
	end;
	      
return;

% M = HU_MOMENTS (IM)
%
% Calculates 7 moments of Hu on image IM, invariant to translation, 
% rotation and scale.
%
% After: M. Sonka et al., Image processing, analysis and machine vision.

function m = hu_moments (im)

	p = [ 1 0 2 1 2 0 3 ];
	q = [ 1 2 0 2 1 3 0 ];

  n = moments(im,p,q,1,1);
   
  m(1) = n(2) + n(3);
  m(2) = (n(3) - n(2))^2   + 4*n(1)^2;
  m(3) = (n(7) - 3*n(4))^2 + (3*n(5) - n(6))^2;
  m(4) = (n(7) +   n(4))^2 + (  n(5) + n(6))^2;
  m(5) = (  n(7) - 3*n(4)) * (n(7) + n(4)) * ...
           (  (n(7) + n(4))^2 - 3*(n(5) + n(6))^2) + ...
         (3*n(5) -   n(6)) * (n(5) + n(6)) * ...
           (3*(n(7) + n(4))^2 -   (n(5) + n(6))^2);
  m(6) = (n(3) - n(2)) * ((n(7) + n(4))^2 - (n(5) + n(6))^2) + ...
          4*n(1) * (n(7)+n(4)) * (n(5)+n(6));      
  m(7) = (3*n(5) -   n(6)) * (n(7) + n(4)) * ...
           (  (n(7) + n(4))^2 - 3*(n(5) + n(6))^2) - ...
         (  n(7) - 3*n(4)) * (n(5) + n(6)) * ...
           (3*(n(7) + n(4))^2 -   (n(5) + n(6))^2);
           
return;

% M = ZERNIKE_MOMENTS (IM, ORDER)
%
% Calculates Zernike moments up to and including ORDER (<= 12) on image IM.
% Default: ORDER = 12.

function m = zernike_moments (im, order)

  if (nargin < 2),             order = 12;                      end;
  if (order < 1 | order > 12), error ('order should be 1..12'); end;

  % xx, yy are grids with co-ordinates

  [xs,ys] = size(im);
  [xx,yy] = meshgrid(-(ys-1)/2:1:(ys-1)/2,-(xs-1)/2:1:(xs-1)/2);

  % Calculate center of mass and distance of any pixel to it

  m  = moments (im,[0 1 0],[0 0 1],0,0);
  xc = m(2)/m(1); yc = m(3)/m(1);
  xx = xx - xc; yy = yy - yc;

  len     = sqrt(xx.^2+yy.^2);
  max_len = max(max(len));

  % Map pixels to unit circle; prevent divide by zero.

  rho        = len/max_len;
  rho_tmp    = rho; rho_tmp(find(rho==0)) = 1;
  theta      = acos((xx/max_len)./rho_tmp);

  % Flip angle for pixels above center of mass

  yneg            = length(find(yy(:,1)<0));
  %disp(find(yy(:,1)<0)')
  %disp([size(xx),size(yy),size(theta),yneg])
  %disp(' ')
  
  %theta(:,1:yneg) = 2*pi - theta(:,1:yneg);
  theta(1:yneg,:) = 2*pi - theta(1:yneg,:);

  % Calculate coefficients

  c = zeros(order,order);
  s = zeros(order,order);

  i = 1;
  for n = 2:order
    for l = n:-2:0
      r    = polynomial (n,l,rho);
      c    = sum(sum(r.*cos(l*theta)))*((n+1)/(pi*max_len^2));
      s    = sum(sum(r.*sin(l*theta)))*((n+1)/(pi*max_len^2));
      m(i) = sqrt(c^2+s^2);
      i    = i + 1;
    end;
  end;

return

function p = polynomial (n,l,rho)

  switch (n)
    case 2, switch (l)
        case 0, p = 2*(rho.^2)-1;
        case 2, p =   (rho.^2);
      end;
    case 3, switch (l)
        case 1, p = 3*(rho.^3)-2*rho;
        case 3, p =   (rho.^3);
      end;
    case 4, switch (l)
        case 0, p = 6*(rho.^4)-6*(rho.^2)+1;
        case 2, p = 4*(rho.^4)-3*(rho.^2);
        case 4, p =   (rho.^4);
      end;
    case 5, switch (l)
        case 1, p = 10*(rho.^5)-12*(rho.^3)+3*rho;
        case 3, p =  5*(rho.^5)- 4*(rho.^3);
        case 5, p =    (rho.^5);
      end;
    case 6, switch (l)
        case 0, p = 20*(rho.^6)-30*(rho.^4)+12*(rho.^2)-1;
        case 2, p = 15*(rho.^6)-20*(rho.^4)+ 6*(rho.^2);
        case 4, p =  6*(rho.^6)- 5*(rho.^4);
        case 6, p =    (rho.^6);
      end;
    case 7, switch (l)
        case 1, p = 35*(rho.^7)-60*(rho.^5)+30*(rho.^3)-4*rho;
        case 3, p = 21*(rho.^7)-30*(rho.^5)+10*(rho.^3);
        case 5, p =  7*(rho.^7)- 6*(rho.^5);
        case 7, p =    (rho.^7);
      end;
    case 8, switch (l)
        case 0, p = 70*(rho.^8)-140*(rho.^6)+90*(rho.^4)-20*(rho.^2)+1;
        case 2, p = 56*(rho.^8)-105*(rho.^6)+60*(rho.^4)-10*(rho.^2);
        case 4, p = 28*(rho.^8)- 42*(rho.^6)+15*(rho.^4);
        case 6, p =  8*(rho.^8)-  7*(rho.^6);
        case 8, p =    (rho.^8);
      end;
    case 9, switch (l)
        case 1, p = 126*(rho.^9)-280*(rho.^7)+210*(rho.^5)-60*(rho.^3)+5*rho;
        case 3, p =  84*(rho.^9)-168*(rho.^7)+105*(rho.^5)-20*(rho.^3);
        case 5, p =  36*(rho.^9)- 56*(rho.^7)+ 21*(rho.^5);
        case 7, p =   9*(rho.^9)-  8*(rho.^7);
        case 9, p =     (rho.^9);
      end;
    case 10, switch (l)
        case  0, p = 252*(rho.^10)-630*(rho.^8)+560*(rho.^6)-210*(rho.^4)+30*(rho.^2)-1;
        case  2, p = 210*(rho.^10)-504*(rho.^8)+420*(rho.^6)-140*(rho.^4)+15*(rho.^2);
        case  4, p = 129*(rho.^10)-252*(rho.^8)+168*(rho.^6)- 35*(rho.^4);
        case  6, p =  45*(rho.^10)- 72*(rho.^8)+ 28*(rho.^6);
        case  8, p =  10*(rho.^10)-  9*(rho.^8);
        case 10, p =     (rho.^10);
      end;
    case 11, switch (l)
        case  1, p = 462*(rho.^11)-1260*(rho.^9)+1260*(rho.^7)-560*(rho.^5)+105*(rho.^3)-6*rho;
        case  3, p = 330*(rho.^11)- 840*(rho.^9)+ 756*(rho.^7)-280*(rho.^5)+ 35*(rho.^3);
        case  5, p = 165*(rho.^11)- 360*(rho.^9)+ 252*(rho.^7)- 56*(rho.^5);
        case  7, p =  55*(rho.^11)-  90*(rho.^9)+  36*(rho.^7);
        case  9, p =  11*(rho.^11)-  10*(rho.^9);
        case 11, p =     (rho.^11);
      end;
    case 12, switch (l)
        case  0, p = 924*(rho.^12)-2772*(rho.^10)+3150*(rho.^8)-1680*(rho.^6)+420*(rho.^4)-42*(rho.^2)+1;
        case  2, p = 792*(rho.^12)-2310*(rho.^10)+2520*(rho.^8)-1260*(rho.^6)+280*(rho.^4)-21*(rho.^2);
        case  4, p = 495*(rho.^12)-1320*(rho.^10)+1260*(rho.^8)- 504*(rho.^6)+ 70*(rho.^4);
        case  6, p = 220*(rho.^12)- 495*(rho.^10)+ 360*(rho.^8)-  84*(rho.^6);
        case  8, p =  66*(rho.^12)- 110*(rho.^10)+  45*(rho.^8);
        case 10, p =  12*(rho.^12)-  11*(rho.^10);
        case 12, p =     (rho.^12);
      end;
  end;

return

