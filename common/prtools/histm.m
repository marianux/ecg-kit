%HISTM Histogramming: mapping of dataset (datafile) to histogram
%
%   W = HISTM(A,N)
%   W = A*HISTM([],N)
%   W = A*HISTM(N)
%   C = B*W
%
%   C = HISTM(B,X)
%   C = B*HISTM(X)
%
% INPUT
%   A    Dataset or datafile for defining histogram bins (training)
%   N    Scalar defining number of histogram bins (default 10)
%   B    Dataset or datafile to be transformed into a histogram with
%        predifined bins.
%   X    Vector with user defined histogram bins (centers)
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
% See IM_HIST for histogramming with known, fixed bin positions.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, MAPPINGS, HIST, IM_HIST

function out = histm(varargin)

  mapname = 'histogramming';
  argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],10);
  [a,bins] = deal(argin{:});
  if mapping_task(argin,'definition')
    if isdouble(bins) && ~isscalar(bins)
      out = define_mapping(argin,'fixed');
    else
      out = define_mapping(argin,'untrained');
    end
    out = setname(out,mapname);
    
  elseif mapping_task(argin,'training')  % training or fixed execution
    
    if isdouble(bins) && ~isscalar(bins) % fixed execution
    
      n = getfeatsize(a,3);
      if isdatafile(a)
        %w = addpostproc(a,histm,length(bins)*n);
        out = addpostproc(a,histm([],bins));
      else
        [m,k] = size(a);
        fsize = getfeatsize(a);
        h = zeros(m,length(bins),n);
        for i=1:m
          im = reshape(+a(i,:),fsize);
          for j=1:n
            imj = im(:,:,j);
            hh = hist(imj(:),bins);
            h(i,:,j) = hh;
          end
        end
        h = reshape(h,m,length(bins)*n);
        out = setdat(a,h);
        out = setfeatsize(out,[length(bins),n]);
      end
    
    else % training
    
      [a,n] = deal(argin{:});
      if n < 3
        error('Number of histogram bins should be larger than 2')
      end
      n = n - 2;  % room for end bins
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
        binsize = (xmax-xmin)/n;  % and compute the bins
        X = [xmin+binsize/2:binsize:xmax-binsize/2]';
      else
        [N,X] = hist(+a,bins);
        n = size(N,1);
      end
      X(1) = X(1)-10*eps; % forces all objects inside edges
      X = [2*X(1)-X(2) ; X ; 2*X(end)-X(end-1)];
      outsize = getfeatsize(a,3)*(n+2);
      w = prmapping(mfilename,'trained',X,[1:outsize]',0,outsize);
      % histm should accept any inputsize
      out = setname(w,mapname);
      
    end
    
  elseif mapping_task(argin,'trained execution')
    
    [a,w] = deal(argin{:});
    bins = getdata(w);
    out = feval(mfilename,a,bins);
    
  else
    
    error('Illegal input')
    
  end

return
			
		
		
		
		
	
