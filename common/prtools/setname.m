%SETNAME Fixed mapping for easy name setting
%
%   A = A*SETNAME(NAME)
%   W = W*SETNAME(NAME)
%
%Set name of dataset A or mapping W

function a = setname(varargin)

  argin = shiftargin(varargin,'char');
  argin = setdefaults(argin,[],[]);
  
  if mapping_task(argin,'definition')
    a = define_mapping(argin,'combiner');
    
  else			% Evaluate
  
    [a,name] = deal(argin{:});
    if isa(a,'prdataset') || isa(a,'prmapping')
      a = setname(a,name);
    else
      error('Illegal input')
    end
  end