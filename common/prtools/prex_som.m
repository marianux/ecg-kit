%PREX_SOM   PRTools example on training SelfOrganizing Maps
%
% Show the training and plotting of 1- or 1-D Self-Organizing Maps.
%
help prex_som

delfigs
echo on
								% Set size of the SOM
	k = [5 1];
								% Set the number of iterations
   nrruns = [20 40 40];
								% Set desired learning rates
	eta = [0.5 0.1 0.1];
	                     % Set the neighborhood widths
   h = [0.6 0.2 0.01];
								% Generate one banana class:
	A = gendatb([100,100]);
	A = seldat(A,1);
								 % Train a 1D SOM:
  W = som(A,k);  % May take some time
								 % Show the results in a scatter plot
	figure(1); clf;
  scatterd(A); hold on;
	plotsom(W);
	title('One-dimensional SOM');
	drawnow
	
								 % Train a 2D SOM:
	k = [5 5];
  W = som(A,k);  % Will take some time
								 % Show the results in a scatter plot
	figure(2); clf;
   scatterd(A); hold on;
	plotsom(W);
	title('Two-dimensional SOM');
	drawnow
	
echo off
showfigs
