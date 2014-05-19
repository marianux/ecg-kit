% MAKEQLD
%
% Script to compile the mex file qld.c and qldorg.c. 

% Store the original path
oldpath = pwd;
% Move to the private dir of prtools
np = which('ldc');
[newpath,fname] = fileparts(np);
cd(fullfile(newpath,'private'));
% Compile
mex -c qldorg.c
if strcmp(computer,'PCWIN')
	mex qld.c qldorg.obj
else
	mex qld.c qldorg.o
end
% Return home
cd(oldpath);
