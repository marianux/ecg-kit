%FINDCLASSES Fixed mapping for finding object indices of all classes
%
%	 C = FINDCLASSES(A,NAME)
%	 C = A*FINDCLASSES([],NAME)
%	 C = A*FINDCLASSES(NAME)
%
% INPUT
%   A      Dataset
%   NAME   Integer: Index of desired labeling, see GETLABLISTNAMES
%          String:  Name of desired labeling, see GETLABLISTNAMES
%          Default: actual LABLIST
%	
% OUTPUT
%   C      Cell array, C{I} contains all object indices to class I.
%
% DESCRIPTION
% C is a cell array such that C{I} contains all object indices to class I.
% A(C{I},:) is a subset of A that contains class I only. Use REMCLASS to
% remove empty classes.
%
% The order of classes can be found by CLASSNAMES. If the name of the class
% is known then GETCLASSI can be used to retrieved its index.
%
% If for A multiple sets of class names (label lists, see MULTI_LABELING)
% are defined, the desired label list can be set by NAME. See also
% GETLABLISTNAMES.
%
% If A has N classes then C has N+1 cells. In the last cell the indices to
% all unlabeled samples are collected.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, SELCLASS, CLASSNAMES, GETCLASSI, REMCLASS, GETLABLISTNAMES,
% MULTI_LABELING

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function b = findclasses(varargin)
  
	argin = shiftargin(varargin,{'char','scalar'});
  argin = setdefaults(argin,[],[]);
  if mapping_task(argin,'definition')
    b = define_mapping(argin,'fixed');
    b = setname(b,'Image stretch');
  else
    [a,name] = deal(argin{:});
    
    if ~isempty(name)
      
      a = changelablist(a,name);
      b = feval(mfilename,a,clas);
      
    else
      
      nlab = getnlab(a);
      [m,k,c] = getsize(a);
      b = cell(1,c+1);
      n = [0 classsizes(a)];
      n = [m-sum(n) n(2:c+1)];
      s = zeros(1,c+1);
      for j=1:c+1
        b{j} = zeros(1,n(j));
      end

      for j=1:m
        p = nlab(j)+1;
        s(p) = s(p)+1;
        b{p}(s(p)) = j;
      end


      b_unlabeld = b(1);
      b(1:c) = b(2:c+1);
      if isempty(b_unlabeld{1})
        b = b(1:c);
      else
        b(c+1) = b_unlabeld;
      end
    end
  end