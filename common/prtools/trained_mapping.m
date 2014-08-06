%TRAINED_MAPPING Define trained mapping
%
%   W = TRAINED_MAPPING(A,DATA,DIM)
%
% INPUT
%   A      - Dataset used for training
%   DATA   - Data (cell araay or structure) to be stored in the data-field
%            of the mapping in order to transfer it to the execution part
%   DIM    - Dimensionality of output space.
%
% OUTPUT
%   W      - Mapping
%
% DESCRIPTION
% This is a simplified version of the definition of a trained mapping. It
% calls PRMAPPING and derives all needed information from the dataset A used
% for training the mapping. In DATA everything should be stored needed for
% the execution of the mapping, either in a structure or by a cell array.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% MAPPINGS, PRMAPPING, TRAINED_CLASSIFIER, DEFINE_MAPPING, MAPPING_TASK

% Copyright: Robert P.W. Duin, prtools@rduin.nl

function w = trained_mapping(varargin)

  [a,data,out_size] = setdefaults(varargin,[],[],0);
  fname = callername;
  mapname = getname(feval(fname));
  if isdataset(a)
    if out_size == 0, out_size = getsize(a,3); end
    w = prmapping(fname,'trained',data,getlablist(a),size(a,2),out_size);
  else
    if out_size == 0, out_size = 1; end
    w = prmapping(fname,'trained',data,[],size(a,2),out_size);
  end
  w = setname(w,mapname);

return

%CALLERNAME
%
%	NAME = CALLERNAME
%
% Returns the name the calling function 

function name = callername

[ss ,i] = dbstack;
if length(ss) < 3
	name = [];
else
	name = ss(3).name;
end
    