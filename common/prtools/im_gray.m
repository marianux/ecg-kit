%IM_GRAY Fixed mapping converting multi-band images to gray images
%
%  B = IM_GRAY(A,V)
%  B = A*IM_GRAY([],V)
%  B = A*IM_GRAY(V)
%
% INPUT
%   A    Multiband image or dataset with multi-band images as objects
%   V    Weight vector, one weight per band. Default: equal weights.
%        (weights will be normalized to sum to one)
%
% OUTPUT
%   B    Output image or dataset.
%
% DESCRIPTION
% The multi-band components in the image A (3rd dimension) or in the
% objects in the dataset A are weigthed (default: equal weights) and
% averaged.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_gray(varargin)
	
	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Color-to-gray conversion');
  else
    [a,v] = deal(argin{:});		
    if isa(a,'prdataset') % allows datafiles too
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
          if length(v) ~= imsize(3)
            error('Number of weights should be equal to number of bands')
          end
          v = v/sum(v);
          b = zeros(imsize(1),imsize(2));
          for j=1:size(a,3)
            b = b + a(:,:,j)*v(j);
          end
        end
      end
    end
  end
return
