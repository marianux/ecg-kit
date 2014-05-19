%SHOW Object display, as images or functions
%
%       H = SHOW(A,N,BACKGROUND)
%
% Displays all objects stored in the datafile A. If the objects
% are images the standard Matlab imagesc-command is used for
% automatic scaling. If the features are images, they are 
% displayed. In case no images are stored in A all objects are
% plotted as functions using PLOTO.
% The number of horizontal figures is determined by N. If N is not
% given an approximately square window is generated.
% Borders between images and empty images are given the value 
% BACKGROUND (default: gray (0.5));
%
% Note that A should be defined by the dataset command, such that
% OBJSIZE or FEATSIZE contains the image size in order to detect 
% and reconstruct the images.
%
% SEE ALSO
% DATASETS, DATA2IM, IMAGE, PLOTO.
