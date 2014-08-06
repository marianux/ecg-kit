%IM_FEATURES Fixed mapping computing image features by the Matlab image toolbox
%
%		F = IM_FEATURES(A,GRAY,FEATURES)
%
% INPUT
%   A        Dataset with binary object images dataset
%   GRAY     Gray-valued images (matched with A, optional)
%   FEATURES Features to be computed
%
% OUTPUT
%   F        Dataset with computed features
%
% This function is a *replacement* for the function IM_MEASURE that is
% based on the DipLib measure.m. This implementation is using the Matlab
% function REGIONPROPS from the image toolbox. Use HELP REGIONPROPS to
% find out which features are exactly supported.
% In each image of the measurement set GRAY the features given in FEATURES 
% are measured. In A a segmented version of GRAY has to be supplied.
% When no GRAY is supplied, the binary images in A are used.
%
% Use FEATURES = 'all' for computing all features.
%
% SEE ALSO
% DATASETS, DATAFILES, REGIONPROPS

% Copyright: D.M.J.Tax, D.M.J.Tax@prtools.org.
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_features(a,gray,features)

	 
if nargin < 3  features = []; end
if nargin < 2  gray = []; end

if nargin < 1 | isempty(a)
   b = prmapping(mfilename,'fixed',{gray,features});
   b = setname(b,'Image features');
elseif isa(a,'prdataset')
   if (nargin < 3 | isempty(features)) & ~isdataset(gray)
      features = gray;
      gray = a;
   end
   if ~strcmp(class(a),class(gray))
      error('Binary and gray images should be both datasets or datafiles')
   end
   fsize = getfeatsize(a);
   if any(getfeatsize(gray) ~= fsize)
      error('Image structures of binary and gray images should be identical')
   end
   if length(fsize) == 2, fsize = [fsize 1]; end
   if size(a,1) ~= size(gray,1)
      error('Same number of binary and gray images expected')
   end
   if isempty(features)
      features = 'all';
   end
   % first generate the correct feature labels, and a structure PROPS
   % that can be used in REGIONPROPS of Matlab:
   im = data2im(a(1,:));
   [props,featlab] = getlabfromregionprops(im,data2im(gray(1,:)),features);
   % RD 27Nov13, compute dat size before
   ftrial = getfeatfromregionprops(im,data2im(gray(1,:)),props);
   out = zeros(size(a,1),numel(ftrial));
   % Run over all images, and compute the features:
   for n = 1:size(a,1)
     out(n,:) = getfeatfromregionprops(data2im(a,n),data2im(gray,n),props);
   end 
   % store in a dataset:
   %b = prdataset(out,getlabels(a));
   b = prdataset(out);
   b = setlablist(b,getlablist(a));
   b = setnlab(b,getnlab(a));
   b = setfeatlab(b,featlab);
   b = setprior(b,getprior(a,0));
   b = setname(b,getname(a));

elseif isa(a,'double') 
   if (nargin < 3 | isempty(features)) & ~isdataset(gray)
      features = gray;
      gray = a;
   end
   [props,featlab] = getlabfromregionprops(a,gray,features);
   b = getfeatfromregionprops(a,gray,props);
else
   error('Wrong input')
end
	
return

% [props,featlab] = getlabfromregionprops(im,gray,features)
%
% Generate a structure PROPS and a cell-array FEATLAB that contains the
% features that are extracted from the images. PROPS is used in the
% Matlab function REGIONPROPS.
function [props,featlab] = getlabfromregionprops(im,gray,features)

% compute the properties:
%keyboard
thr = min(im(:)) + (max(im(:))-min(im(:)))/2;
r = regionprops((im>thr),gray,features);
props = fieldnames(r);
% we do it complicated like this, because bwprops could have been 'all',
% and by this we now get all possible (black and white) properties

% remove the things that we are not interested in:
% RD 27Nov13, PixelValues added
rmfeatures = {'PixelList', 'SubarrayIdx', 'ConvexHull','ConvexImage', ...
   'Image', 'FilledImage', 'PixelIdxList', 'Extrema', 'PixelValues'};
props = setdiff(props,rmfeatures);

% generate the feature labels
featlab = {};
nr = 0;
for i=1:length(props)
   nrf = size(getfield(r,props{i}),2);
   if nrf>1 % there are more measurements per single property
      for j=1:nrf
         nr = nr+1;
         featlab{nr} = [props{i},num2str(j)];
      end
   else % one measurement per property
      nr = nr+1;
      featlab{nr} = props{i};
   end
end
featlab = strvcat(featlab);

return

%   f = getfeatfromregionprops(im,gray,props)
%
% Compute the features F from the images IM. Which features are computed
% is defined in PROPS.
function f = getfeatfromregionprops(im,gray,props)

thr = min(im(:)) + (max(im(:))-min(im(:)))/2;
r = regionprops((im>thr),gray,props);
f = [];
for i=1:length(props)
   f = [f getfield(r,props{i})];
end

