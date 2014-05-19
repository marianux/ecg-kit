%PREX_CONFMAT PRTools example on confusion matrix, scatterplot and gridsize
%
% Prtools example code to show the use of confusion matrix, 
% scatterplot and gridsize.

help prex_confmat

delfigs
echo on
                % Load 8-class 2D problem
  randn('state',1);
  rand('state',1);
  a = gendatm;
                % Compute the Nearest Mean Classifier
  w = nmc(a);
                % Scatterplot
  figure;
  gridsize(30);
  scatterd(a,'legend');
                % Plot the classifier
  plotc(w);
  title([getname(a) ', Gridsize 30']);
                % Set higher gridsize
  gridsize(100);
  figure;
  scatterd(a,'legend');
  plotc(w);
  title([getname(a) ', Gridsize 100']);

                % Classify training set
  d = a*w;
                % Look at the confusion matrix and compare it to the scatterplot
  confmat(d);

echo off
showfigs
c = num2str(gridsize);
disp(' ')
disp('   Classifier plots are inaccurate for small gridsizes. The standard');
disp('   value of 30 is chosen because of the speed, but it is too low to');
disp('   ensure good plots. Other gridsizes may be set by gridsize(n).')
disp('   Compare the two figures and appreciate the difference.')

