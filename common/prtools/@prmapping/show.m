%SHOW Display axes of affine mappings as images, if available
%
%    SHOW(W,N,BACKGROUND)
%
% If W is a affine mapping operating in a space defined by images
% (i.e. each object in the space is an image) and the image size is
% properly stored in W (SIZE_IN), then the images corresponding to the axes
% of the affine mapping are displayed.
%
% The number of horizontal images is determined by N. If N is not given an
% approximately square window is generated.
%
% Borders between images and empty images are given the value 
% BACKGROUND (default: gray (0.5));
