function a = reorderdset(a,L,flag)
%REORDERDSET
%
% DESCRIPTION
% This is just needed to execute the below command outside the @dataset
% directory, as it fails there in some Matlab versions

if nargin < 3, flag = 1; end

if flag
	a(L,:) = a;
else
	a = a(L,:);
end
	
return