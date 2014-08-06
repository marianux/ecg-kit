%SELCLASSFixed mapping for selecting a single class from a dataset
%
%	[B,J] = SELCLASS(A,CLASS,NAME)
%	[B,J]  = A*SELCLASS([],CLASS,NAME)
%	[B,J]  = A*SELCLASS(CLASS,NAME)
%
% INPUT
%   A      Dataset
%   CLASS  Integer: Indices of desired classes in CLASSNAMES(A)
%          String array:  Class names
%          Cell array: Indices of desired classes in CLASSNAMES(A)
%          Default C = {}, i.e. return all classes separated out in a cell
%          array.
%   NAME   Integer: Index of desired labeling, see GETLABLISTNAMES
%          String:  Name of desired labeling, see GETLABLISTNAMES
%          Default: actual LABLIST
%	
% OUTPUT
%   B   Desired classes of the dataset A. In case CLASS is a cell array, B
%       is a cell array of the desired classes. In case CLASS is empty, B
%       is a cell array of all classes.
%   J   Indices of returned objects in dataset A: B = A(J,:). In case CLASS
%       is a cell array, J is a cell array as well.
%
% DESCRIPTION
% B is a subset of the dataset A defined by the set of classes (CLASS). In
% case of a multi-labeling system (see MULTI_LABELING) the desired CLASS
% should refer to the label list NAME.
%
% In case A is soft labeled or is a target dataset by B = SELCLASS(A,CLASS)
% the entire dataset is returned, but the labels or targets are reduced to
% the selected class (target) CLASS.
%
% In case CLASS is a cell array the outputs are organised as cell arrays.
%
% EXAMPLES
% a = gendatm; b = selclass(a,[2 3 4]); % selects 3 classes
% a = gendatm; b = a*selclass;         % returns every class in a cell
% a = gendatb; 
% a = addlabels(a,genlab(25*ones(4,1)),'4class'); % add second label list
% b = a*selclass('4class'); % returns 4 cells, preserves label list.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, SELDAT, CLASSNAMES, GENDAT, GETLABLIST, GETCLASSI, REMCLASS,
% GETLABLISTNAMES, MULTI_LABELING

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function [b,J] = selclass(varargin)
  
  argin = shiftargin(varargin,{'vector','cell','char'});
  argin = shiftargin(argin,'char',2);
  argin = setdefaults(argin,[],{},[]);
  
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    
  else			% Evaluate
  
    [a,clas,lablist] = deal(argin{:});
    isa(a,'prdataset');
    if ~isempty(lablist)
      curn = curlablist(a);
      a = changelablist(a,lablist);
      [b,J] = feval(mfilename,a,clas);
      if iscell(b)
        for n=1:numel(b)
          b{n} = changelablist(b{n},curn);
        end
      else
        b = changelablist(b,curn);
      end
    else
      if iscell(clas)
        if  isempty(clas)
          for n=1:getsize(a,3);
            clas{n} = n;
          end
        end
        if isdataset(a)
          b = prdataset;
        else
          b = prdatafile;
        end
        J = cell(1,numel(clas));
        b = cell(1,numel(clas));
        for n=1:numel(clas)
          [b{n},J{n}] = seldat(a,clas{n});
        end
      else
        [b,J] = seldat(a,clas);
      end
    end
    
  end