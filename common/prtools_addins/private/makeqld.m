mex -c qldorg.c
if strcmp(computer,'PCWIN') || strcmp(computer,'PCWIN64')
	mex qld.c qldorg.obj
else
	mex qld.c qldorg.o
end

