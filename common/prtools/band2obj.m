%BAND2OBJ Mapping image bands to objects
%
%   B = BAND2OBJ(A,N)
%   B = A*BAND2OBJ([],N)
%   B = A*BAND2OBJ(N)
%
% INPUT
%   A   Dataset or datafile with multiband image objects.
%   N   Number of successive bands to be combined in an object.
%       The number of image bands in A should be multiple of N.
%       Default N = 1.
%
% OUTPUT
%   B   Output dataset or datafile.
%
% DESCRIPTION
% If the objects in a dataset or datafile A are multi-band images, e.g. RGB
% images, or the result of IM_PATCH, then the featsize of A is [C,R,L] for
% for L bands of a C x R image. This routine combines sets of N successive 
% bands as separate objects. The total number of objects is thereby
% enlarged by a factor L/N. All information of the constituting objects
% like labels, is copied to the newly created objects.
%
% Note: BAND2OBJ cannot be applied to datafiles for which already a
% bandselection (BANDSEL) has been defined.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, BANDSEL, IM2OBJ, DATA2IM

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = band2obj(varargin)

	argin = shiftargin(varargin,'vector');
  argin = setdefaults(argin,[],1);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Bands to Objects');
  else
    [a,N] = deal(argin{:});	
    isvaldfile(a);
    %isobjim(a); % datafiles have object images
    m = size(a,1);
    if isempty(N) % all needed for treating variable numbers of bands
      bandnames = getident(a,'bandnames');
      if ~isempty(bandnames)
        K = zeros(m,1);
        for j=1:m, K(j,:) = size(bandnames{j},1); end
        if any(K~=K(1))
          N = zeros(m,max(K));
          for j=1:m, N(j,1:K(j)) = [1:K(j)]; end
        else
          N = 1;
        end
      else
        N = 1;
      end
    end
    if size(N,1) == m  & m > 1 % per object all bands to be selected are given
      k = max(N(:));
      b = [];
      s = sprintf('Checking %i bands: ',k);
      prwaitbar(k,s);
      for j=1:k
        prwaitbar(k,j,[s int2str(j)]);
        [I,J] = find(N==j); % Finds all objects for which we have to select band j
        if ~isempty(I)
          b = [b; bandsel(a(I,:),repmat(j,length(I),1))]; % wrong !!!! bandsel expects something else, solve there!
        end
      end
      prwaitbar(0);
    elseif size(N,1) == 1 % N is number of successive objects to be combined
      fsize = getfeatsize(a);
      k = fsize(3); % We assume that all objects have the same number of bands
      if k > 1
        if N*round(k/N) ~= k
          error('Number of images bands should be multiple of N')
        end
        size1 = [size(a,1),fsize(1),fsize(2),N,fsize(3)/N];
        size2 = [size1(2:4),size1(5)*size(a,1)];
        b = im2obj(reshape((shiftdim(reshape(+a,size1),1)),size2));
        lablista = getlablist(a);
        nlaba = getnlab(a);
        nlabb = repmat(nlaba',k/N,1);
        b = prdataset(b,lablista(nlabb(:),:));
        labim = genlab(ones(1,size(a,1))*(k/N));
        curn = curlablist(b);
        b = addlabels(b,labim,'org_image');
        b = changelablist(b,curn);
        b = setident(b,getident(a(str2num(labim),:)));
      else
        b = a;
      end
    else
      error('Second parameter has wrong size')
    end
  end


return
