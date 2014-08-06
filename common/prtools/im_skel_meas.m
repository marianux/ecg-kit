%IM_SKEL_MEASURE Fixed mapping measuring skeleton-based features
%
%    F = IM_SKEL_MEASURE(A,FEATURES)
%    F = A*IM_SKEL_MEASURE([],FEATURES)
%    F = A*IM_SKEL_MEASURE(FEATURES)
%
% INPUT
%   A        Dataset with binary object images dataset
%   FEATURES Features to be computed
%
% OUTPUT
%   F        Dataset with computed features
%
% DESCRIPTION
% The following features may be computed on the skeleton images in A:
% 'branch', 'end', 'link', 'single'. They may be combined in a cell
% array or given separately.
%
% Use FEATURES = 'all' for computing all features (default).
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, DATAFILES, IM_SKEL, DIP_IMAGE, BSKELETON

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function b = im_skel_meas(varargin)
	argin = shiftargin(varargin,{'char','cell'});
  argin = setdefaults(argin,[],'all');
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Skeleton features');
  else
    features = cell(1,numel(argin)-1);
    [a,features{:}] = deal(argin{:});
    if iscell(features{1})
      features = features{1};
    elseif strcmp(features{1},'all')
    	features = {'branch', 'end', 'link', 'single'};
    end
    if isa(a,'prdataset') % allows datafiles too
      isobjim(a);
      b = filtim(a,mfilename,{features});
      b = setfeatlab(b,features);
    elseif isa(a,'double') || isa(a,'dip_image') % here we have a single image
      a = double(a);
      b = [];
      if ~iscell(features), features = {features}; end
      cleanp  = bwmorph(a,'clean');  
      endp    = bwmorph(cleanp,'endpoints');  
      branchp = bwmorph(cleanp,'branchpoints'); 
      se = sum(endp(:));
      sb = sum(branchp(:));
      st = sum(a(:));
      ss = st-sum(cleanp(:));
      sl = st-se-sb-ss;
      for i = 1:length(features)
        switch features{i}
          case 'branch'
            b = [ b sb ];
          case 'end' 
            b = [ b se ];
          case 'link'
            b = [ b sl ];
          case 'single'
            cleanp = bwmorph(a,'clean');  
            b = [ b ss ];
          otherwise
            error('Unknown feature desired: %s',features{i})
        end
      end
    end
  end
	
return
