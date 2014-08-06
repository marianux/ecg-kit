%IM_STRETCH Fixed mapping for contrast stretching of images (DIP_Image)
%
%	B = IM_STRETCH(A,LOW,HIGH,MIN,MAX)
%	B = A*IM_STRETCH([],LOW,HIGH,MIN,MAX)
%	B = A*IM_STRETCH(LOW,HIGH,MIN,MAX)
%
% INPUT
%   A        Dataset with object images dataset (possibly multi-band)
%   LOW      Lower percentile (default 0)
%   HIGH     Highest percentile (default 100)
%   MIN      Miniumum (default 0)
%   MAX      Maximum (default 255)
%
% OUTPUT
%   B        Dataset with filtered images
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, STRETCH

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_stretch(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],0,100,0,255);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image stretch');
  else
    [a,low,high,min,max] = deal(argin{:});
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,{low,high,min,max});
    elseif isa(a,'double') | isa(a,'dip_image') % here we have a single image
      if checktoolbox('dipimage')
        a = 1.0*dip_image(a);
        b = stretch(a,low,high,min,max);
      else
        diplibwarn
        low_high = stretchlim(a,[low high]/100);
        b = imadjust(a,low_high,[]);
      end;
    end
  end
	
return
