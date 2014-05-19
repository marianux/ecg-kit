%IM_GRAY Conversion of multi-band images into gray images
%
%  B = IM_GRAY(A,V)
%  B = A*IM_GRAY([],V)
%
% INPUT
%   A    Multiband image or dataset with multi-band images as objects
%   V    Weight vector, one weight per band. Default: equal weights.
%
% OUTPUT
%   B    Output image or dataset.
%
% DESCRIPTION
% The multi-band components in the image A (3rd dimension) or in the
% objects in the dataset A are weigthed (default: equal weights) and
% averaged.
%
% SEE ALSO
% MAPPINGS, DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_gray(a,v)
	
		
  if nargin < 2, v = []; end
  if nargin < 1 | isempty(a)
    b = prmapping(mfilename,'fixed',{v});
    b = setname(b,'Color-to-gray conversion');
	elseif isa(a,'prdataset') % allows datafiles too
		isobjim(a);
    b = filtm(a,mfilename,v);
    imsize = getfeatsize(a);
    b = setfeatsize(b,imsize(1:2));
  else
		a = double(a);
    imsize = size(a);
    if isempty(v)
      if length(imsize) == 3
        b = mean(a,3);
        b = squeeze(b);
      elseif length(imsize) == 2
        b = a;
      else
        error('Illegal image size')
      end
    else
      if length(imsize) == 2
        b = a;
      else
        b = zeros(imsize(1),imsize(2),size(a,1));
        for i=1:size(a,1)
          for j=1:size(im,3)
            b(:,:,i) =b(:,:,i) + b(:,:,j,i)*v(j);
          end
        end
      end
    end
  end
return
