%IM_PATCH Fixed mapping or generating patches from images
%
%   B = IM_PATCH(A,PSIZE,PNUM,TYPE)
%   B = IM_PATCH(A,PSIZE,COORD,'user')
%   B = A*IM_PATCH([],PSIZE,PNUM,TYPE)
%   B = A*IM_PATCH(PSIZE,PNUM,TYPE)
%
% INPUT
%   A        Dataset or datafile with (multi-band) object images dataset
%   PSIZE    2-dimensional patch size. If PSIZE is a scalar, square patches
%            of size PSIZE x PSIZE are generated. In case PSIZE < 1, it is
%            taken relative to the images size. In this case patches may
%            become non-square for non-square images. Default 3 x 3.
%   PNUM     Number of patches, see TYPE, default 1.
%   COORD    Given set of N x 2 image coordinates of patch centra. This may
%            be given either in pixels (ANY(COORD > 1) or relative to the 
%            image size (ALL(COORD <= 1).
%   TYPE     'syst': systematic sampling generating PNUM x PNUM patches.
%            'rand': generation of a random set of PNUM patches.
%            'user': user defined patch positions, see COORD.
%            Default: 'syst'.
%
% OUTPUT
%   B        Resulting dataset or datafile with the same number of objects
%            as in A. Single images are replaced by the patches.
%
% DESCRIPTION
% The object images (including their N bands) are sampled and windows of size
% PSIZE are generated, replacing the original image object. They are stored
% as [PSIZE(1) PSIZE(2) NUM*N] image objects, in which N is the original
% number of bands and NUM is either PNUM (for TYPE is 'rand' or 'user') or
% PNUM x PNUM (for TYPE is 'syst'). 
%
% By BAND2OBJ(B,N) individual patches can be transformed into objects.
%
% EXAMPLES
% b = im_patch(a,0.5) % generate just the central parts of the images 
% b = im_patch(a,[3 5],2) % generates 4 patches of size 3x5 in centra of
%                     % the four quadrants
% b = im_patch(a,0.1,10,'rand') % generate at random positions 10 patches
%                     %each with a linear size of 0.1 of the image size.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, BAND2OBJ

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_patch(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[3 3],1,'syst',[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image patches');
  else
    [a,psize,pnum,option,N] = deal(argin{:});	
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      if strcmp(option,'rand')
        pnum = rand(pnum,2);  % generate desired number of points
        option = 'user';
      end
      b = filtim(a,mfilename,{psize,pnum,option},0);
    else                          % here we have a single image
      a = double(a);
      asize = size(a);
      if length(psize) == 1, psize = [psize psize]; end
      if all(psize < 1), psize = psize.*asize; end
                                % user want to know how many patches
      if ~isempty(N) & N == 0   % this is just a hack, could be smarter
        if strcmp(option,'syst')
          if length(pnum) == 1, npatches = pnum*pnum;
          else npatches = prod(pnum); 
        end
        elseif strcmp(option,'user'), npatches = size(pnum,1);
        elseif strcmp(option,'random'), npatches = pnum;
        else error('Unknown option');
        end
        b = npatches;
        return
      end

      rsize = round(psize);
      if strcmp(option,'syst')
        if length(pnum) == 1, pnum = [pnum pnum]; end
        if isempty(N)
          patches = zeros(rsize(1),rsize(2),pnum(1)*pnum(2));      
        else
          patches = zeros(rsize(1),rsize(2),length(N));
        end
        for i=1:2
          if pnum(i)*psize(i) <= asize(i)
            %s(i) = asize(i)/(2*pnum(i))-psize(i)/2+0.5;
            s(i) = asize(i)/(2*pnum(i))-psize(i)/2+1;
            d(i) = asize(i)/pnum(i);
          else
            %s(i) = 0.5;
            s(i) = 1;
            d(i) = (asize(i)-psize(i))/(pnum(i)-1);
          end
        end
        if isempty(N)
          for j2 = 1:pnum(2)
            for j1 = 1:pnum(1)
              disp(size(a))
              disp(round(s(1)+(j1-1)*d(1)))
              disp(round(s(1)+(j1-1)*d(1))+rsize(1)-1)
              disp(round(s(2)+(j2-1)*d(2)))
              disp(round(s(2)+(j2-1)*d(2))+rsize(2)-1)
              aa = a(round(s(1)+(j1-1)*d(1)):round(s(1)+(j1-1)*d(1))+rsize(1)-1, ...
                  round(s(2)+(j2-1)*d(2)):round(s(2)+(j2-1)*d(2))+rsize(2)-1);
              patches(:,:,(j2-1)*pnum(1)+j1) = aa;
            end
          end
        else
          for n = 1:length(N)
            j2 = floor((N(n)-1)/pnum(1))+1;
            j1 = N(n) - pnum(1)*(j2-1);
            aa = a(round(s(1)+(j1-1)*d(1)):round(s(1)+(j1-1)*d(1))+rsize(1)-1, ...
                round(s(2)+(j2-1)*d(2)):round(s(2)+(j2-1)*d(2))+rsize(2)-1);
            patches(:,:,n) = aa;
          end
        end
      elseif strcmp(option,'rand')
        locs = rand(pnum,2);
        if ~isempty(N)
          locs = locs(N,2);
        end
        patches = feval(mfilename,a,psize,locs,option);
      elseif strcmp(option,'user')
        locs = pnum;
        if ~isempty(N)
          locs = locs(N,:);
        end
        pnum = size(locs,1);
        patches = zeros(rsize(1),rsize(2),pnum);
        if all(locs(:) <= 1) % positions given as relative coordinates
          locs = locs.*repmat([asize(1)-psize(1),asize(2)-psize(2)],pnum,1);
          locs = locs + repmat([0.5 0.5],pnum,1);
        end
        border = max(ceil(rsize/2));
        a = bord(a,NaN,border); % create mirror border to avoid problems (takes time!!)
        locs = locs + border;
        for j = 1:pnum
          s = round(locs(j,:)-psize/2);
          patches(:,:,j) = a(s(1):s(1)+rsize(1)-1,s(2):s(2)+rsize(2)-1);
        end
      end
      b = patches;	
      end
  end
	
return
