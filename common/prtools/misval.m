%MISVAL Fixed mapping handling the missing values in a dataset
%
%    B = MISVAL(A,VAL)
%    B = A*MISVAL([],VAL)
%    B = A*MISVAL(VAL)
%
% INPUT
%    A    Dataset, containing NaNs (missing values)
%    VAL  String with substitution option
%         or value used for substitution
%
%    B    Dataset with NaNs substituted
%
% DESCRIPTION
%
% The following values for VAL are possible:
%   'remove'    remove objects (rows) that contain missing values (default)
%   'f-remove'  remove features (columns) that contain missing values
%   'mean'      fill the entries with the mean of their features
%   'c-mean'    fill the entries with the class mean of their features
%   <value>     fill the entries with a fixed constant
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

function [x,msg] = misval(varargin)

  mapname = 'Missing values';
	argin = shiftargin(varargin,{'scalar','char'});
  argin = setdefaults(argin,[],[]);
  
  if mapping_task(argin,'definition')
    
    x = define_mapping(argin,'fixed',mapname);
    
  elseif mapping_task(argin,'training')			% Compute mapping.
    
    [x,val] = deal(argin{:});
    x = prdataset(x);
    [m,k,c] = getsize(x);
    % Where are the offenders?
    I = isnan(x);
    % If there are missing values, go:
    if any(I(:))
      switch val
        case {'remove' 'delete'}
          J = find(sum(I,2)==0);
          x = x(J,:);
          msg = 'Objects with missing values have been removed.';
        case {'f-remove' 'f-delete'}
          J = find(sum(I,1)==0);
          x = x(:,J);
          msg = 'Features with missing values have been removed.';
        case 'mean'
          for i=1:k
            J = ~I(:,i);
            if any(I(:,i)) %is there a missing value in this feature?
              if ~any(J)
                error('Missing value cannot be filled: all values are NaN.');
              end
              mn = mean(x(J,i));
              x(find(I(:,i)),i) = mn;
            end
          end
          msg = 'Missing values have been replaced by the feature mean.';
        case 'c-mean'
          for j=1:c
            L = findnlab(x,j);
            for i=1:k
              J = ~I(L,i);
              if any(I(L,i)) %is there a missing value in this feature for this class?
                if ~any(J)
                 error('Missing value cannot be filled: all values are NaN.');
                end
                mn = mean(x(J,i));
                x(find(I(:,i)),i) = mn;
              end
            end
          end
          msg = 'Missing values have been replaced by the class feature mean.';
        otherwise
          if isstr(val)
            error('unknown option')
          end
          if ~isa(val,'double')
            error('Missing values can only be filled by scalars.');
          end
          x(I) = val;
          msg = sprintf('Missing values have been replaced by %f.',val);
      end
    end
  end
  
return
