%IM_THRESHOLD Fixed mapping thresholding images (DIP_Image)
%
%	B = IM_THRESHOLD(A,TYPE,PAR,INV)
%	B = A*IM_THRESHOLD([],TYPE,PAR,INV)
%	B = A*IM_THRESHOLD(TYPE,PAR,INV)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%   TYPE     Type of procedure, see below
%   PAR      Related parameter
%   INV      If INV = 1, result inverted, default INV = 0.
%
% OUTPUT
%   B        Dataset with thresholded images
%
% DESCRIPTION
% The following procedures are supported (TYPE)
% 'isodata':      Thresholding using the Isodata algorithm
% 'otsu':         See IM2BW
% 'triangle':    *Thresholding using chord method
%                 (a.k.a. skewed bi-modality, maximum distance to triangle)
%                 by Zack, Rogers and Latt (1973).
% 'background':  *Thresholding using unimodal background-symmetry method.
% 'fixed':        Thresholding at a fixed value.
% 'double':      *Thresholding between two fixed values.
% 'volume':      *Thresholding to obtain a given volume fraction.
% 'hysteresis':  *From the binary image (in>low) only those regions are
%                 selected for which at least one pixel is (in>high)
%
% The following parameters are related to these procedures (PAR):
% ('background'): Distance to the peak where we cut-off, in
%                 terms of the half-width at half the maximum.
%                 Inf selects the default value, which is 2.
% ('fixed'):      Threshold value. Inf means halfway between
%                 minimum and maximum value.
% ('double'):     Two threshold values. Inf means min+[1/3,2/3]*(max-min).
% ('volume'):     Parameter = the volume fraction (Inf means 0.5)
% ('hysteresis'): Two values: low, high (see above)
%
% * routine needs DIP_IMAGE
%
% This routine is for smart thresholding. Simple thresholding can also be
% performed by B = A > PAR, which is the same as using the type 'fixed'.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, THRESHOLD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_threshold(varargin)

	argin = shiftargin(varargin,'char');
  argin = setdefaults(argin,[],'isodata',inf,0);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image threshold');
  else
    [a,type,par,inv] = deal(argin{:});
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,{type,par,inv});
      b = setfeatsize(b,getfeatsize(a));
    elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
      switch lower(type)
        case('isodata')
          mina = min(a(:));  maxa = max(a(:));
          b = 255*(a-mina)/(maxa-mina);
          t0 = 128;
          for j=1:25
            U = b >  t0;
            L = b <= t0;
            n2 = mean(b(U)); n1 = mean(b(L));
            t1 = (n1+n2)/2;
            if t1 == t0;
              break;
            else
              t0 = t1;
            end
          end
          b = b > t0;
        case('fixed')
          b = a > par;
        case('otsu')
          b = im2bw(a);
        case {'triangle','background','double','volume','hysteresis'}
          checktoolbox('dipimage');       
          a = 1.0*dip_image(a);
          b = threshold(a,type,par);
        otherwise
          error 'Illegal thresholding type'
      end
      if inv
        b = 1-b;
      end
    end
  end

return