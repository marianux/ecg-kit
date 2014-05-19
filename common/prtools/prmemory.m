%PRMEMORY Set/get size of memory usage
% 
%   N = PRMEMORY(N)
%
%   N : The desired / retrieved maximum size data of matrices (in
%       matrix elements)
% 
% DESCRIPTION
% This retoutine sets or retrieves a global variable GLOBALPRMEMORY that
% controls the maximum size of data matrices in PRTools. Routines like 
% KNNC, KNN_MAP, PARZEN_MAP and TESTP make use of it by computing
% additional loops and avoiding to define very large distance matrices.
% The default for this maximum size is set to 10000000. For most
% computers there is no need to reduce it. There is not much speed up 
% to be expected if it is enlarged.
%
% PRMEMORY gives also the number of elements for which a conversion from
% datafiles to datasets is approved.
%
% This facility is illustrated by the following example, using the
% routine PRMEM.
%
% Assume that an array of the size [M x K] has to be computed. The 
% numbers of LOOPS and ROWS are determined which are needed such that 
% ROWS*K < GLOBALPRMEMORY (a global variable that is initialized in this 
% routine, if necessary). The final number of rows for the last loop 
% is returned in LAST. 
%
% EXAMPLES
% [M,K] = size(A);
% [LOOPS,ROWS,LAST] = prmem(M,K);
% if (LOOPS == 1)
%  RESULT = < compute the result based on A >
% else
%   RESULT = 0;
%   for J =1:LOOPS
%     if (J == LOOPS), N = LAST; else N = ROWS; end
%     NN = (J-1)*ROWS;
%     RESULT = RESULT + < compute a partial result based on A(NN+1:NN+N,:) >
%	  end
% end

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

% $Id: prmemory.m,v 1.2 2007/04/13 09:30:54 duin Exp $

function n_out = prmemory(n_in)

	persistent GLOBALPRMEMORY;

	if (isempty(GLOBALPRMEMORY))
		GLOBALPRMEMORY = 50000000;
	end
	
	if nargin > 0
    if ischar(n_in)
      n_in = str2num(n_in);
    end
		GLOBALPRMEMORY = n_in;
	end
	
	if nargout > 0
		n_out = GLOBALPRMEMORY;
	elseif nargin == 0
		disp(GLOBALPRMEMORY)
	end
	
return;
