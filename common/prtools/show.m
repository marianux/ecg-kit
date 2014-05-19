%SHOW PRTools general show
%
%   H = SHOW(A,N,B)
%
% INPUT
%   A      Image
%   N      Number of images on a row
%   B      Intensity value of background (default 0.5);
%
% OUTPUT 
%   H      Graphics handle
%
% DESCRIPTION
% PRTools offers a SHOW command for variables of the data classes DATASET
% and PRDATAFILE. In order to have a simliar command for images not converted
% to a DATASET this commands made availble. A should be 2D, 3D or 4D image.
%
% 2D images are fully displayed.
%
% 3D images are converted to a dataset with as many feature images as given
% in the 3rd dimension and displayed by DATASET/SHOW.
%
% 4D images with 3 bands in the 3rd dimension are converted to a dataset with 
% as many 3-color object images as are given in the 4th dimension and
% displayed by DATASET/SHOW
%
% All other 4D images are converted to a dataset with as many 2D feature
% images as given by the dimensions 3 and 4 and displayed by DATASET/SHOW.
% Unless given otherwise, N is set to size(A,3).
%
% SEE ALSO DATASET/SHOW DATAFILE/SHOW

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function h = show(a,n,background)

	
if nargin < 3, background = 0.5; end
if nargin < 2, n = []; end
a = double(a);
s = size(a);

if length(s) == 2
	if any(s==1)
		error('Image expected')
	end
	a = im2obj(a);
elseif length(s) == 3 & s(3) ~= 3
	a = im2feat(a);
	a = prdataset(a,NaN); % avoid display label image
else
	if length(s) >= 3 & s(3) == 3
		a = im2obj(a,s(1:3));
	else
		a = reshape(a,s(1),s(2),prod(s(3:end)));
		a = im2feat(a);
		a = prdataset(a,0); % avoid display label image
	end
end
if nargout > 0
	h = show(a,n,background);
else
	show(a,n,background);
end
