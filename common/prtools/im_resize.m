%IM_RESIZE Fixed mapping for resizing object images
%
%  B = IM_RESIZE(A,SIZE,METHOD)
%  B = A*IM_RESIZE([],SIZE,METHOD)
%  B = A*IM_RESIZE(SIZE,METHOD)
%
% INPUT
%  A       Dataset or datafile
%  SIZE    Desired size
%  METHOD  Method, see IMRESIZE
%
% OUTPUT
%  B       Dataset or datafile
%
% DESCRIPTION
% The objects stored as images in the dataset or datafile A are resized
% using the IMRESIZE command. Default METHOD is 'nearest' (nearest 
% neighbor interpolation). 
%
% A special method is 'preserve', which exactly copies the existing data, 
% cuts it, or extends it with zeros as far as appropriate. In SIZE the 
% desired output size has to be stored. Note that for proper use in 
% PRTools third size parameter of multi-band images, e.g. 3 for RGB 
% images, has to be supplied. 
%
% In case SIZE is a scalar the default METHOD is 'preserve', which 
% implies that the first SIZE samples are taken, which is useful if
% A is a 1-D signal. Otherwise the deafult METHOD is 'nearest'.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, DATAFILES, IM2OBJ, DATA2IM 

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

%DXD 24-8-2007
%  I rewrote a part of this function. Now there are default values
%  given, a bug is removed, and the identation is correct again.

function b = im_resize(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[16,16],[],[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image resize');
  else
    [a,imsize,method,par] = deal(argin{:});	
    if isempty(method)
      if length(imsize) == 1
        method = 'preserve';
      else
        method = 'nearest';
      end
    end
    
    if isa(a,'prdataset') % allows datafiles too

      isobjim(a);
      b = filtim(a,mfilename,{imsize,method},imsize); % prepare execution image by image

    else

      b = double(a);

      if strcmp(method,'preserve') | strcmp(method,'preserve_bottom') | ...
          strcmp(method,'preserve_centre') | strcmp(method,'preserve_top')
        % copy the image pixel by pixel into a larger or smaller image

        sizeb = size(b);               % original size
        if length(imsize) == 1
          imsize = [1 imsize];
        end
        if length(sizeb) == 3 & length(imsize) == 2
          imsize(3) = sizeb(3);
        end
        if length(imsize) ~= length(sizeb)
          error('Desired images size should have as many dimension as data')
        end
        sizec = min(imsize,sizeb);     % size of part to be copied
        subs = cell(1,length(sizec));  % store indices in cell array
        [subs{:}] = ind2sub(sizec,[1:prod(sizec)]);
        subsb = subs; subsc = subs;    % cell arrays for original and result

        switch method
          case 'preserve_bottom'
            delc  = imsize-sizeb;
          case 'preserve_centre'
            delc = round((imsize-sizeb)/2);
          case {'preserve_top','preserve'}
            delc = zeros(1,length(imsize));
        end

        for j=1:length(delc)            % indices for
          if delc(j) > 0                % result
            subsc{j} = subs{j} + delc(j);
          elseif delc(j) < 0            % and original
            subsb{j} = subs{j} - delc(j);
          end
        end

        Lc = sub2ind(imsize,subsc{:}); % corresponding linear coordinates of result
        Lb = sub2ind(sizeb,subsb{:});  % corresponding linear coordinates of original
        c = zeros(imsize);             % embed result in zeros
        c(Lc) =  b(Lb);                % copy
        b = c;                         % store result

      else % resizing the image using Matlab's imresize

        if length(imsize) > 1
          b = imresize(b,imsize(1:2),method);
        elseif imsize > 1 & round(imsize) == imsize
          b = imresize(a,[imsize imsize],method);
        else
          [m,n] = size(a);
          b = imresize(a,round(imsize*[m,n]),method);
        end

      end

    end
  end
  return

