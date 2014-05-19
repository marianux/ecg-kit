function UpdStructure = pls_updstruct(OldStructure, NewStructure)
%
% Structure update
%
%   UpdStructure = updstruct(OldStructure, NewStructure)
%   
% In the result UpdStructure all new filed from NewStructure will be added, 
% coinciding fields will be replaced, fields which do not exist in the NewStruct
% will be preserved

% Copyright: S.Verzakov, serguei@ph.tn.tudelft.nl
% Faculty of Applied Sciences, Delft University of Technology
% P.O. Box 5046, 2600 GA Delft, The Netherlands

UpdStructure = OldStructure;
if nargin < 2 | isempty(NewStructure) | ~isstruct(NewStructure)
  return;
end

fn = fieldnames(NewStructure);
for i=1:length(fn)
  v = getfield(NewStructure,fn{i});

  if isa(v,'struct'); 
    if isfield(UpdStructure, fn{i}) 
      w = getfield(UpdStructure, fn{i});  
      if isa(w,'struct')
        v = updstruct(w,v); 
      end  
    end
  end 
  
  UpdStructure = setfield(UpdStructure,fn{i},v);
end

return;

