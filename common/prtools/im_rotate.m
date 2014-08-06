%IM_ROTATE Fixed mapping for image rotation 
%
%	B = IM_ROTATE(A,ALF)
%	B = A*IM_ROTATE([],ALF)
%	B = A*IM_ROTATE(ALF)
%
% INPUT
%   A        Dataset with object images (possibly multi-band)
%   ALF      Rotation angle (in radians), 
%            default: rotation to main axis
%
% OUTPUT
%   B        Dataset with rotated object images 
%
% DESCRIPTION
% The objects stored as images in the dataset or datafile A are rotated
% using the IMROTATE command.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_rotate(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image rotate');
  else
    [a,alf] = deal(argin{:});	
    if isa(a,'prdataset') % allows datafiles too
%       error('Command cannot be used for datasets as it may change image size')
%     elseif isdatafile(a)
      isobjim(a);
      b = filtim(a,mfilename,{alf});
      b = setfeatsize(b,getfeatsize(a));
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image

      a = double(a);

      if isempty(alf)
        m = im_moments(a,'central',[1 2 0; 1 0 2]');
        C = [m(2) m(1); m(1) m(3)];
        [E,D] = preig(C); [dummy,ind] = sort(diag(D));
        alf = atan2(E(2,ind(1)),E(1,ind(1)))+pi/2;
      end

      b = imrotate(a,alf*360/(2*pi),'bilinear','crop');

    end
    
  end

return
