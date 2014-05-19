%TRAINED_MAPPING Define untrained or fixed mapping
%
%   W = TRAINED_CLASSIFIER(A,DATA)
%
% INPUT
%   A      - Dataset used for training
%   DATA   - Data (cell array or structure) to be stored in the data-field
%            of the mapping in order to transfer it to the execution part
%
% OUTPUT
%   W       - Classifier
%
% DESCRIPTION
% This routine serves as a simplified definition of a trained classifier.
% It sets automatically the name, the label list and the size. In DATA
% everything should be stored needed for the execution of the mapping, 
% either in a structure or by a cell array.
%
% SEE ALSO
% MAPPINGS, PRMAPPING, TRAINED_MAPPING, DEFINE_MAPPING, MAPPING_TASK

% Copyright: Robert P.W. Duin, prtools@rduin.nl

function w = trained_classifier(varargin)

  [a,data] = setdefaults(varargin);
  fname = callername;
  classfname = getname(feval(fname));
  [m,k,c] = getsize(a);
  w = prmapping(fname,'trained',data,getlablist(a),k,c);
  w = setname(w,classfname);

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
    