%MAPPINGTOOLS Macro defining some mappings:
%
%  ARGMIN, ARGMIN1, ARGMIN2, ARGMAX, ARGMAX1, ARGMAX2
%  EXPM, NEXPM,
%  ARGMAXMAX, ARGMAXMIN, ARGMINMAX, ARGMINMIN, ARGMEANMAX, ARGMEANMIN
%  MINK1, MINK2
%
% This is not a function. Calling MAPPINGTOOLS just loads the above
% mappings ready to use. They can be called as J = A*ARGMIN in which A is a
% dataset or an array of doubles. It returns a vector of object indices J
% pointing to the feature (column) with the lowest value.

argmin = mapm('min')*out2;
argmin1 = mapm('min',{[],1})*out2;
argmin2 = mapm('min',{[],2})*out2;
argmax = mapm('max')*out2;
argmax1 = mapm('max',{[],1})*out2;
argmax2 = mapm('max',{[],2})*out2;
expm = mapm('exp');
nexpm = mapm('uminus')*mapm('exp');
argmaxmax = mapm('max',{[],1})*argmax2; 
argmaxmin = mapm('max',{[],1})*argmin2; 
argminmax = mapm('min',{[],1})*argmax2; 
argminmin = mapm('min',{[],1})*argmin2; 
argmeanmax = mapm('mean',1)*argmax2; 
argmeanmin = mapm('mean',1)*argmin2; 
mink2 = mapm('distm')*mapm('sqrt');
mink1 = proxm([],'m',1)*mapex;

