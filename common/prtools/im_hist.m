%IM_HIST Histogramming of image dataset (datafile)
%
%    H    = IM_HIST(A,N)
%    H    = A*IM_HIST([],N)
%    H    = A*IM_HIST(N)
%    H    = IM_HIST(A,X)
%    H    = A*IM_HIST([],X)
%    H    = A*IM_HIST(X)
%   
% INPUT
%   A    Dataset or datafile
%   N    Number of bins, default 256
%   X    User defined histogram bins (centers)
%
% OUTPUT
%   H    Dataset or datafile with histogram bin frequencies
%
% DESCRIPTION
% For every object in the dataset A the set of feature values (for images 
% pixels)is mapped into a histogram, specifying for each bin the number of
% features having a value as specified for that bin. If H is converted to a
% dataset it has a feature size equal to the number of bins times the
% number of bands in the image. So color images return a histogram for
% every band.
%
% Use HISTM if the bin positions have to be determined automatically is a
% 'training' stage from some images.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, MAPPINGS, HIST, HISTM

function h = im_hist(varargin)

argin = shiftargin(varargin,'double');
[a,bins] = setdefaults(argin,[],256);
if mapping_task(argin,'definition')
  h = define_mapping(argin,'fixed','histogram');
elseif isdatafile(a)
  h = addpostproc(a,im_hist(bins));
elseif isdataset(a)
  if numel(bins) == 1
    h = filtim(a,'imhist',bins);
  else
    m = size(a,1);
    im = data2im(a(1,:));
    nbands = size(im,3);
    h = zeros(1,numel(bins),nbands,m);
    t = sprintf('Histogramming %i objects: ',m);
    prwaitbar(m,t);
    for i=1:m
      prwaitbar(m,i,[t num2str(i)]);
      im = data2im(a(i,:));
      for j=1:nbands
        band = im(:,:,j);
        h(1,:,j,i) = hist(band(:),bins);
      end
    end
    prwaitbar(0);
    h = setdata(a,im2obj(h,[1,numel(bins),nbands]));
  end
else
  error('Illegal input')
end
	
			
		
		
		
		
	
