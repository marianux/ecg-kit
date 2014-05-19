%PKLIBSVC Automatic radial basis SVM, using nu_libsvc and the Parzen kernel
%
%   W = PKLIBSVC(A,ALF)
%
% INPUT
%   A   Dataset
%   ALF Parameter, default 1
%
% OUTPUT
%   W   Mapping: radial basis support vector classiifer
%
% DESCRIPTION
% This routine provides a radial basis support vector classifier based on
% NULIBSVC (which estimates NU using the leave-one-out 1NN error) and
% estimates the kernel width SIGMA by the the value found by PARZENC. The
% kernel width used is ALF*3*SQRT(2)*SIGMA. This is much faster than the
% gridsearch used by RBLIBSVC and performs about equal.
%
% SEE ALSO
% DATASETS, MAPPINGS, LIBSVC, NULIBSVC, RBLIBSVC, PARZENC

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com

function w = pklibsvc(a,alf)

if nargin < 2 | isempty(alf), alf = 1; end
if nargin < 1 | isempty(a)
	w = prmapping(mfilename,'untrained',alf);
  w = setname(w,'PK-LIBSVM');
	return
end

if size(a,1) <= 1000
  [v,sig] = parzenc(a);
else
  [v,sig] = parzenc(gendat(a,1000));
end
w = nulibsvc(a,proxm([],'r',alf*3*sig*sqrt(2)));

