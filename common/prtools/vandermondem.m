%VANDERMONDEM Fixedd mapping, extend data matrix
%
%    Z = VANDERMONDEM(X,N)
%    Z = X*VANDERMONDEM([],N)
%    Z = X*VANDERMONDEM(N)
%
% INPUT
%    X    Data matrix
%    N    Order of the polynomail
%
% OUTPUT
%    Z    New data matrix containing X upto order N
%
% DESCRIPTION
% Construct the Vandermonde matrix Z from the original data matrix X by
% including all orders upto N. Note that also order 0 is added:
%    Z = [ones  X  X^2  X^3 ... X^N]
% This construction allows for the trivial extension of linear methods
% to obtain polynomail regressions.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
%   LINEARR

% Copyright: D.M.J. Tax, D.M.J.Tax@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
function z = vandermondem(varargin)
  
	mapname = 'Vandermonde map';
  argin = shiftargin(varargin,'scalar');
  argin = setdefaults(argin,[],1);
  
  if mapping_task(argin,'definition')
    z = define_mapping(argin,'fixed',mapname);
    
  elseif mapping_task(argin,'training')	% Execute.
  
    [x,n] = deal(argin{:});
    dat = +x;
    [m,dim] = size(dat);
    I = 1:dim;
    z = ones(m,n*dim+1);
    for i=0:(n-1)
      z(:,(i+1)*dim+I) = dat.*z(:,i*dim+I);
    end
    % remove the superfluous ones:
    z(:,1:dim-1) = [];
    z = setdat(x,z);
    
  end

return
