%OUTM Mapping that reorders outputs of input mapping
%
%   V = W*OUTM(N)
%
% INPUT
%   W   Mapping
%   N   Integer or vector with desired output order, default N = 2
%
% OUTPUT
%   V   Mapping
%
% DESCRIPTION
% In case A*W (A is a dataset or double) returns several outputs then
% A*(W*OUTM(N)) returns them in order as given by N
%
% EXAMPLE
% test2 = testc*out2;% define testc for second output par (# class errors)
% a = gendatd;       % train set
% t = gendatd;       % test set
% t*knnc(a,1)*test2  % execute, list of # errors per class
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function varargout = outm(varargin)

varargout = cell(1,nargout);
argin = setdefaults(varargin,[],2);
if (~isempty(argin{1})) && (nargout == 1) && (isa(argin{1},'double')),
  argin = {[] argin{1}};
end
if mapping_task(argin,'definition')
  varargout{1} = define_mapping(argin,'combiner','ReorderOut');
elseif mapping_task(argin,'combiner')
  varargout{1} = define_mapping({[],argin{:}},getmapping_type(argin{1}),getname(argin{1}));
else
  argout = cell(1,max(argin{3}));
  [argout{:}] = argin{1}*argin{2};
  varargout = argout(argin{3});
end