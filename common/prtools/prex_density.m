%PREX_DENSITY Various density plots
%
% Prtools example to show the use of density estimators and how to
% visualize them.

help prex_density

delfigs
figure
echo on

  % Generate one-class data
  a = gencirc(200);

  % Parzen desity estimation
  w = parzenm(a); 

  % scatterplot
  subplot(2,2,1);
  scatterd(a,[10,5]);
  plotm(w); 

  title('Parzen Density')

  % 3D density plot
  subplot(2,2,2);
  scatterd(a,[10,5]);
  plotm(w,3);

  % Mixture of Gaussians (5)
  w = gaussm(a,5);

  % scatterplot
  subplot(2,2,3);
  scatterd(a,[10,5]);
  plotm(w);
  title('Mixture of 5 Gaussians')

  % 3D density plot
  subplot(2,2,4);
  scatterd(a,[10,5]);
  plotm(w,3);
  drawnow

  disp([newline '  Study figure at full screen, shrink and hit return'])
  
  pause

  figure

  % Define, name and store four density esimators
  W1 = gaussm;       W1 = setname(W1,'Gaussian');
  W2 = gaussm([],2); W2 = setname(W2,'Mixture of 2 Gaussians');
  W3 = parzenm;
  W4 = knnm([],10);   W4 = setname(W4,'10-Nearest Neighbor');
  W = {W1 W2 W3 W4};

  % generate data
  a = +gendath;

  % plot densities and estimator name
  for j=1:4
    subplot(2,2,j)
    scatterd(a,[10,5])
    plotm(a*W{j})
    title([getname(W{j}) ' density estimation'])
  end

echo off
showfigs
