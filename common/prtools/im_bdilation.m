%IM_BDILATION Fixed mapping for binary dilation (DIP_Image)
%
%	B = IM_BDILATION(A,N,CONNECTIVITY,EDGE_CONDITION)
%	B = A*IM_BDILATION([],N,CONNECTIVITY,EDGE_CONDITION)
%	B = A*IM_BDILATION(N,CONNECTIVITY,EDGE_CONDITION)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   N        Number of iterations (default 1)
%   CONNECTIVITY    See BDILATION
%   EDGE_CONDITION  Value of edge, default 0
%
% OUTPUT
%   B        Dataset with dilated images
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, BDILATION

% Copyright: R.P.W. Duin, duin@ph.tn.tudelft.nl
% Faculty of Applied Physics, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

function b = im_bdilation(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1,-2,0);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image dilation');
  else
    [a,n,connect,edgecon] = deal(argin{:});
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,{n,connect,edgecon});
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      if checktoolbox('dipimage')
        a = dip_image(a,'bin');
        b = bdilation(a,n,connect,edgecon);
      else
        diplibwarn
        b = bwmorph(a,'dilate',n);
      end
    else
      error('Illegal input')
    end
	end
	
return
