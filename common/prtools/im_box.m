%IM_BOX Fixed mapping determining a rectangular enclosing a blob (0/1 image)
%
%   B = IM_BOX(A)
%   B = A*IM_BOX
%
% If A is a 0/1 image then B is the same image with all empty (0) border
% columns and rows removed.
%
%   B = IM_BOX(A,N)
%
% If A is a 0/1 image then B is the same image, but having in each direction
% N empty (0) columns and rows around the object (1).
%
%   B = IM_BOX(A,[NX1 NX2 NY1 NY2])
%
% If A is a 0/1 image then B is the same image, but having NX1, NX2 empty 
% columns (0) left, respectively right of the object (1) and NY1, NY2 empty
% rows (0) above, respectively below the object(1).
%
%   B = IM_BOX(A,N,ALF)
%
% Adds as many empty (0) columns or rows such that the aspect ratio of
% images (height/width) equals ALF. For ALF == 1, square images are returned.
% For ALF == 0, images are taken as they are and N rows and columns are
% added.
%
% EXAMPLE
% prdatafiles;            % make sure prdatafiles is in the path
% x = kimia_images;       % load kimia images
% x = x*im_box(0);        % remove all empty rows and columns
% x = x*im_box(0,1);      % add rows/columns to make images square
% x = x*im_resize([32 32]); % resample them 
% x = x*im_box(1,0);      % add rows/columns and keep image square
% % now all images are 34x34 and no object touches the border.
% show(gendat(x,4))       % show a few
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [b,J] = im_box(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],[],[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image bounding box');
  else
    [a,n,alf] = deal(argin{:});
    if isdataset(a)
      error('Command cannot be used for datasets as it may change image size')
    elseif isdatafile(a)
      isobjim(a);
      b = filtim(a,mfilename,{n,alf});
      %b = setfeatsize(b,getfeatsize(a));
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      if isa(a,'dip_image'), a = double(a); end
      if isempty(n)
        jx = find(sum(a,1) ~= 0);
        jy = find(sum(a,2) ~= 0);
        J = [min(jy),max(jy),min(jx),max(jx)];
        b = a(min(jy):max(jy),min(jx):max(jx));
      else
        if (~isempty(alf) && alf == 0)
          c = a;
        else
          c = feval(mfilename,a); 
        end
        [my,mx] = size(c);
        if length(n) == 1
        	n = [n n n n];
        elseif length(n) ~= 4
        	error('Second parameter should be scalar or vector of length 4')
        end
        b = zeros(my+n(3)+n(4),mx+n(1)+n(2));
        b(n(3)+1:n(3)+my,n(1)+1:n(1)+mx) = c;
      end
      if ~isempty(alf) && alf ~= 0
        [m,k] = size(b);
        r = round(m*alf) - k;
        if r == 0
          ;
        elseif r >= 1 % add r columns
          c = zeros(m,k+r);
          c(:,ceil(r/2):ceil(r/2)+k-1) = b;
          b = c;
        else % add rows
          r = round(k/alf) - m;
          c = zeros(m+r,k);
          c(ceil(r/2):ceil(r/2)+m-1,:) = b;
          b = c;
        end
      end
		end	 
	end

return