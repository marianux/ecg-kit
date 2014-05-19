%HIST Display feature histograms
%
%     HANDLE = HIST(A,P,NX)
%     HANDLE = HIST(+A)
%
% INPUT
%  A       Dataset
%  P       # of bins
%  NX      # of histograms displayed in a row
%
% OUTPUT
%  HANDLE  handle of subplot
% 
% DESCRIPTION
% For all feature (columns) of A a histogram is plot using P bins.
% These histograms are plot as subplots in a single figure, displaying
% NX histograms in a row. In HANDLE the handles of the subplots are 
% returned.
%
% Note that this routine is not a true overload of the HIST command.
% Use HIST(+A) if that is desired.
