%REMCLASS Fixed mapping to remove small or empty classes
%
%   B = REMCLASS(A,N)
%   B = A*REMCLASS([],N)
%   B = A*REMCLASS(N)
%
% INPUT
%   A   Dataset
%   N   Integer, maximum class size to be removed (optional; default 0)
%
% OUTPUT
%   B   Dataset
%
% DESCRIPTION
% Classes having N objects or less are removed. The corresponding objects
% are made unlabeled. Use SELDAT to remove unlabeled objects.
% In case of soft labeled objects the number of objects in A is compared
% with N. If it is smaller all object labels are removed (NaN).
%
% A*REMCLASS removes empty classes.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, SELDAT, GENDAT

function b = remclass(varargin)
  
	mapname = 'Remove classes';
  argin = shiftargin(varargin,'integer');
  argin = setdefaults(argin,[],0);
  
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed',mapname);
    
  else 			% Execution
    
    [a,n] = deal(argin{:});	
    N = classsizes(a);
    J = find(N <= n);
    L = findnlab(a,J);
    b = setnlab(a,0,L);
    b = setlablist(b);
    
  end
	
return