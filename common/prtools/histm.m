%HISTM Histogramming: mapping of dataset (datafile) to histogram
%
%   W = HISTM(A,N)
%   W = A*HISTM([],N)
%   C = B*W
%
%   C = HISTM(B,X)
%   C = B*HISTM([],X)
%   
%   
%
% INPUT
%   A    Dataset or datafile for defining histogram bins (training)
%   N    Scalar defining number of histogram bins (default 10)
%   B    Dataset or datafile to be transformed into a histogram with
%        predifined bins.
%   X    User defined histogram bins (centers)
%
% OUTPUT
%   C    Dataset or datafile with histogram bin frequencies
%
% DESCRIPTION
% For every object in the dataset B the set of feature values is mapped
% into a histogram, specifying for each bin the number of features having a
% value as specified for that bin. This is particular useful if the objects
% are images and the features are pixels. In that case for every image a
% histogram is found.
%
% The dataset A may be used to find the proper histogram bins. In that case
% histograms with N bins are constructed between the minimum and maximum
% values over all objects in A.
%
% Formally HISTM([],N) is an untrained mapping, to be trained by A as the 
% dataset (datafile) A is used to determine the histogram bin centers. 
% In case the bins are given like in HISTM(B,X) then we have a trained mapping.
% Consequently, if A is a datafile then in C = A*HISTM(A,10) all objects in
% A are processed twice. Once for determining the bin positions and once for
% filling them. If appropriate a command like C = A*HISTM(A(1,:),10) is
% thereby significantly faster, as it determines the bin positions by just
% a single object.
%
% SEE ALSO
% DATASETS, DATAFILES, MAPPINGS, HIST

function w = histm(a,bins)

fixed = 0;
if nargin < 2, bins = 10; end
if nargin < 1, a = []; end
if isdouble(bins) & ~is_scalar(bins)
	fixed = 1;
end
if ismapping(bins) 
	bins = getdata(bins); %DXD: I don't know how to solve it
	fixed = 1;
end

mapname = 'histogramming';

if fixed                  %  fixed histograms
	
	if isempty(a)
		%w = prmapping(mfilename,'fixed',bins,bins(:),0,length(bins));
		%DXD make it a trained mapping, you know the output size
		%    Unfortunately, you  don't really know, because for color
		%    images you get a histogram *per* color band. Therefore I
		%    decided not to set the feature labels.
		w = prmapping(mfilename,'trained',bins,[],0,0); % size_out is unknown
		                             % as it depends on the number of bands
		w = setname(w,mapname);
	else
		n = getfeatsize(a,3);
		if isdatafile(a)
			%w = addpostproc(a,histm,length(bins)*n);
			w = addpostproc(a,histm([],bins));
		else
			[m,k] = size(a);
			fsize = getfeatsize(a);
			h = zeros(m,length(bins),n);
			for i=1:m
				im = reshape(+a(i,:),fsize);
				for j=1:n
					imj = im(:,:,j);
					hh = hist(imj(:),length(bins));
					h(i,:,j) = hh;
				end
			end
			h = reshape(h,m,length(bins)*n);
			w = setdat(a,h);
		end
	end
	
else                        % adjustable histograms
	
	if isempty(a)                  % defining
		w = prmapping(mfilename,'untrained',bins);
		w = setname(w,mapname);
	elseif ~ismapping(bins)        % training
		if bins < 3
			error('Number of histogram bins should be larger than 2')
		end
		bins = bins - 2;  % room for end bins
		if isdatafile(a)  % oeps! training from datafile
			next = 1;       % we just find minimum and maximum
			xmin = inf;     % of all feature values in all datafiles
			xmax = -inf;    
			while next > 0
				[b,next] = readdatafile(a,next);
				b = +b;
				xmin = min(min(b(:)),xmin);
				xmax = max(max(b(:)),xmax);
			end
			binsize = (xmax-xmin)/bins;  % and compute the bins
			X = [xmin+binsize/2:binsize:xmax-binsize/2]';
			n = bins;
		else
			[N,X] = hist(+a,bins);
			n = size(N,1);
		end
		X(1) = X(1)-10*eps; % forces all objects inside edges
		X = [2*X(1)-X(2) ; X ; 2*X(end)-X(end-1)];
		outsize = getfeatsize(a,3)*(n+2);
		w = prmapping(mfilename,'trained',X,[1:outsize]',0,outsize);
    % histm should accept any inputsize
		w = setname(w,mapname);
	else                           % execution, bins is a mapping
		w = feval(mfilename,a,getdata(bins)); % use fixed mapping
	end
	
end

return
			
		
		
		
		
	
