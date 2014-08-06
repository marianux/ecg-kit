%NLABCMP Compare two label lists and count the differences
% 
%   [N,C] = NLABCMP(LAB1,LAB2)
% 
% INPUT
%   LAB1, 
%   LAB2  Label lists
%
% OUTPUT
%   C     A false/true vector pointing to different/equal labels
%   N     Number of differences in LAB1 and LAB2
%
% DESCRIPTION  
% Compares two label lists and counts the disagreements between them.

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function [N,C] = nlabcmp(S1,S2)

  [m,k] = size(S1);
	[n,l] = size(S2);
  if isempty(S1)
    if isempty(S2)
      N = 0;
      C = [];
    else
      N = n;
      C = false(n,1);
    end
  else
    if isempty(S2)
      N = m;
      C = false(m,1);
    else
      if (m ~= n)
        error('Label list sizes do not match.')
      end
      if (iscell(S1) ~= iscell(S2))
        error('Label lists should be both cells, strings or numeric.')
      end

      if (iscell(S1))
        C = strcmp(S1,S2);
      elseif (all(size(S1) == size(S2)))
        C = all(S1'==S2',1)';
      else
        C = strcmp(cellstr(S1),cellstr(S2));
      end

      N = m - sum(C);
    end
  end
  
return
