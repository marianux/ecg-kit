%IM_LABEL Fixed mapping labeling binary images (DIP_Image)
%
%	B = IM_LABEL(A,CON,MIN_SIZE,MAX_SIZE)
%	B = A*IM_LABEL([],CON,MIN_SIZE,MAX_SIZE)
%	B = A*IM_LABEL(CON,MIN_SIZE,MAX_SIZE)
%
% INPUT
%   A        Dataset with binary object images dataset (possibly multi-band)
%   N        Number of iterations (default 1)
%   CON      Connectivity, see LABEL
%   MIN_SIZE Minimum size of objects to be labeled (default 0)
%   MAX_SIZE Maximum size of objects to be labeled (default 0: all)
%
% OUTPUT
%   B        Dataset with labeled images
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, DIP_IMAGE, LABEL

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_label(varargin)

	argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],2,0,0);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image labeling');
  else
    [a,connect,minsize,maxsize] = deal(argin{:});	
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,{connect,minsize,maxsize});
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      if checktoolbox('dipimage')
        a = dip_image(a,'bin');
        b = double(label(a,connect,minsize,maxsize));
      else
        diplibwarn
        b = bwlabel(a,connect*4);
      end
    end
  end
	
return
