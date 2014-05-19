%DYADICM Dyadic dataset mapping
%
%   B = [C1,C2]*DYADICM([],P,Q,SIZE)
%   B = [C1,C2]*DYADICM([],FNAME,[],SIZE)
%   W = A*DYADICM([],{U1,U2,FNAME},[],SIZE)
%   B = A*DYADICM([],{V1,V2,FNAME},[],SIZE)
%
% INPUT
%  C1,C2   Datasets / datafiles to be combined
%  A       Input dataset used for training or classification
%  P,Q     Scalar multiplication factors (default 1) to compute P*A1+Q*A2
%          or string (name of a routine)
%  FNAME   String with the name of a function to combine two datasets or
%          mappings, 'plus' for PLUS(U1,U2)
%  U1,U2   Untrained mappings to be combined
%  V1,V2   Trained or fixe mappings to be combined
%  SIZE    Desired images size of output dataset objects
%
% OUTPUT
%  B       Dataset
%  W       Combined mapping
%
% DESCRIPTION
% This special mapping is a low-level routine to facilitate dyadic 
% operations on datafiles, datasets and mappings. Datasets to be combined 
% should have the same number of objects. Image objects should have the 
% same image size. Datafiles should be preprocessed or postprocessed 
% versions of the same original datafile. Mappings should have the same
% input and output size. 
%
% This routine has been written for use by PRTools programmers only. Users
% are discouraged to call it directly. The routine is called by the 
% dyadic operations of the classes 'prmapping' and 'prdatafile'.
%
% SEE ALSO
% DATASETS, MAPPINGS, DATAFILES

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = dyadicm(a,p,q,fsize)
  
  if nargin < 4, fsize = []; end
  if nargin < 3, q = []; end
  if nargin < 2, p = []; end
  
  if nargin < 1 | isempty(a)
    % note that this routine is sometimes externally set to a combiner
    b = prmapping(mfilename,'fixed',{p,q,fsize});
    b = setsize_out(b,prod(fsize));
    b = setname(b,'dyadicm');
    return
  end
  
  if isempty(p)
    p = 1; 
  end
      
  if iscell(a) & ismapping(a{1})
    % combining mappings. We are just here because the dyadic mapping 
    % operations like mapping/plus.m are programmed like this. 
    if isfixed(a{1})
      b = prmapping(mfilename,'fixed',{{a{1},a{2},p}});
    elseif isuntrained(a{1})
      b = prmapping(mfilename,'untrained',{{a{1},a{2},p},[],[]});
    elseif istrained(a{1})
      if isempty(fsize), fsize = getsize(a{1},2); end
      b = prmapping(mfilename,'trained',{a{1},a{2},p},getlabels(a{1}),getsize(a{1},1));
    else
      error('Wrong type')
    end
  
  elseif isdatafile(a)
    nodatafile; % forces to store this routine as postprocessing
    
  elseif iscell(p) & isuntrained(p{1}) % train mappings stored in p
    v1 = p{1};     % untrained mapping
    v2 = p{2};     % fixed or trained mapping or constant   
    p = p{3};      % operation
    v1 = a*v1;     % train first mapping
    [kin,kout] = size(v1);
    if ismapping(v2) & (isfixed(v2) | istrained(v2))
      ;            % leave it and hope for the best
    elseif ismapping(v2) & isuntrained(v2)
      v2 = a*v2;   % train it
    else           % should be scalar constant
      ;            % leave it 
    end
    % store a standard trained mapping
    b = prmapping(mfilename,'trained',{v1,v2,p},getlabels(v1),kin,kout);
    
  elseif iscell(p) % execute mappings stored in p
    v1 = p{1};     % fixed or trained mapping
    v2 = p{2};     % fixed or trained mapping or constant   
    p = p{3};      % operation
    a1 = a*v1;     % prepare datasets
    if ismapping(v2)
      a2 = a*v2;
    else
      a2 = setdata(a1,v2*ones(size(a1)));
    end
    b = feval(p,a1,a2);
    
  elseif ismapping(p) % execute standard PRTools trained mapping
    v1 = getdata(p,1);     % fixed or trained mapping
    v2 = getdata(p,2);     % fixed or trained mapping or constant   
    p = getdata(p,3);      % operation
    a1 = a*v1;     % prepare datasets
    if ismapping(v2)
      a2 = a*v2;
    else
      a2 = setdata(a1,v2*ones(size(a1)));
    end
    b = feval(p,a1,a2);
      
  else % the basic call by data and parameters
    
    % The dataset A has horizontally to be split in two datasets A1, A2. 
    % If P and Q are scalars or if Q = [], this is done half-half.
    % If P is a string (name of a routine) it is used to combine A1 and A2.
    
    [a1,a2,fsize] = split_dataset(a,p,q,fsize);

    if ischar(p) % function name in p
      b = feval(p,a1,a2);  
    else % addition
      if isempty(q), q = 1; end
      b = a1*p + a2*q;
    end
    
    if isdataset(b)
      b = setfeatsize(b,fsize);
    end
    
  end
  
return

function [a1,a2,fsize] = split_dataset(a,p,q,fsize)

%   if iscell(a)
%     
%     a1 = a{1};
%     a2 = a{2};
%     
%   %elseif ischar(p) & ~isempty(q)
% %   elseif ~isempty(q)
% %     % split as defined by q
% %     fsize = q;
% %     a1 = a(:,1:q);
% %     a2 = a(:,q+1:end);
%     
%   else
    % 50-50 split, but preserve possible image structure
    if isempty(fsize) | fsize == 0
      if isdataset(a)
        fsize = getfeatsize(a);
      else % may be a is set of images
        a = double(a);
        fsize = size(a);
      end
      while (fsize(end) == 1) & (length(fsize) > 1)
        fsize = fsize(1:end-1);
      end
      fsize(end) = fsize(end)/2;
    end
    k = size(a,2);
  
    if k ~= 2*floor(k/2)
      error('Feature size of dataset should be multiple of 2')
    end
        
    a1 = a(:,1:k/2);
    a2 = a(:,k/2+1:k);
    if isdataset(a)
      a1 = setfeatsize(a1,fsize);
      a2 = setfeatsize(a2,fsize);
    end
    
%   end
  
return
