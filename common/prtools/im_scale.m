%IM_SCALE Fixed mapping scaling binary images to a giving fraction of pixels 'on'
%
%   B = IM_SCALE(A,P)
%   B = A*IM_SCALE([],P)
%   B = A*IM_SCALE(P)
%
% B is a zoomed in / out version of A such that about a fraction
% P of the image pixels is 'on' (1).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM_BOX, IM_CENTER

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_scale(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],0.5);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image scale');
  else
    [a,p] = deal(argin{:});	
    if isdataset(a)
      error('Command cannot be used for datasets as it may change image size')
    elseif isdatafile(a)
      isobjim(a);
      b = filtim(a,mfilename,{p});
      b = setfeatsize(b,getfeatsize(a));
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      sca = sqrt(p/mean(a(:)));
      sa = size(a);
      c = imresize(double(a),round(sca*sa),'nearest');
      sc = size(c);
      d = abs(floor((sc - sa)/2));
      if sca < 1
        b = zeros(size(a));
        b(d(1)+1:d(1)+sc(1),d(2)+1:d(2)+sc(2)) = c;	
      else
        b = c(d(1)+1:d(1)+sa(1),d(2)+1:d(2)+sa(2));
      end
    end
  end

return