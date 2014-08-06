%IM_BEROSION Fixed mapping for binary erosion (DIP_Image)
%
%	B = IM_BEROSION(A,N,CONNECTIVITY,EDGE_CONDITION)
%	B = A*IM_BEROSION([],N,CONNECTIVITY,EDGE_CONDITION)
%	B = A*IM_BEROSION(N,CONNECTIVITY,EDGE_CONDITION)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   N        Number of iterations (default 1)
%   CONNECTIVITY    See BEROSION
%   EDGE_CONDITION  Value of edge, default 1
%
% OUTPUT
%   B        Dataset with eroded images
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, BEROSION

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function b = im_berosion(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1,-2,1);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image erosion');
  else
    [a,n,connect,edgecon] = deal(argin{:});
    if isa(a,'prdataset') % allows datafiles too		
      isobjim(a);
      b = filtim(a,mfilename,{n,connect,edgecon});
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      if checktoolbox('dipimage')
        a = dip_image(a,'bin');
        b = berosion(a,n,connect,edgecon);
      else
        diplibwarn
        b = bwmorph(a,'erode',n);
      end
    else
      error('Illegal input')
    end
	end
	
return
