%SCATTERN Simple 2D scatterplot of dataset without axis annotation
%
%   SCATTERN(A,LEGEND)
%   A*SCATTERN([],LEGEND)
%   A*SCATTERN(LEGEND)
%
% INPUT
%   A       Dataset
%   LEGEND  Logical, true/false for including legend, default false
%
% DESCRIPTION
% A simple, unannotated 2D scatterplot is created, without any axes.
%
% SEE ALSO
% DATASETS, SCATTERD

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function varargout = scattern(varargin)

varargout = {};
argin = shiftargin(varargin,'logical');
argin = shiftargin(argin,'scalar');
argin = setdefaults(argin,[],false,[]);
if mapping_task(argin,'definition')
  % standard return, name = filename
  varargout = {define_mapping(argin,'fixed')};
elseif mapping_task(argin,'fixed execution')
  % a call like w = template(a,parsin)
  [a,plotlegend,cmap] = deal(argin{:});
  if isdataset(a), nlab = getnlab(a); else, nlab = []; end
  h = gscatter(+a(:,1),+a(:,2),nlab,cmap,repmat('.',1,4),9,'off');
  axis equal
  axis tight
  axis off
  if plotlegend
    legend(h,getlablist(a,'string'))
  end
  if nargout > 0
    varargout = {h};
  end
end
  