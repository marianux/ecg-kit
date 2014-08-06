%RESIZEM Mapping for resizing object images in datasets and datafiles
%(outdated, replaced by im_resize)
%
%  B = RESIZEM(A,SIZE,METHOD)
%  B = A*RESIZEM([],SIZE,METHOD)
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
% using the IMRESIZE command. Default METHOD is 'nearest' (nearest neighbor
% interpolation. In SIZE the desired output size has to be stored. Note
% that for proper use in PRTools third size parameter of multi-band images,
% e.g. 3 for RGB images, has to be supplied. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, DATASETS, DATAFILES, IM2OBJ, DATA2IM 

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = resizem(a,imsize,method)

if nargin < 3 | isempty(method)
  method = 'nearest';
end

if isempty(a)
  b = prmapping(mfilename,'fixed',{imsize,method});
  b = setname(b,'image resize');
  b = setsize_out(b,imsize);
  return
end

isvaldfile(a,0);
if isdataset(a)
  isobjim(a);
end
b = filtm(a,'imresize',{imsize(1:2),method},imsize);
  

