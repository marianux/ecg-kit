%PREX_PARZEN Parzen based denisities and classifiers
%
% PRTools example to show the differences between various ways to use the
% PARZEN procedures for estimating densities and classifiers.

help prex_parzen

delfigs
figure
echo on

  delfigs
  a = gendath;   % two normally distributed classes, different covariances
  w = a*parzenc; % Parzen classifier, single smoothing parameter optimizing
                 % the classification error
  figure(1); scatterd(a); % show scatterplot
  plotm(w);  plotc(w);    % show densities and classifier
  title('Densities and classifier by PARZENC')
  w = a*parzendc;% Parzen classifier, smoothing parameter per class
                 % optimizing class densities
  figure(2); scatterd(a); % show scatterplot
  plotm(w);  plotc(w);    % show densities and classifier
  title('Densities and classifier by PARZENDC')
  w = a*parzenm; % Parzen density, smoothing parameter per class
                 % optimizing class densities, combined to single density
  figure(3); scatterd(a); % show scatterplot
  plotm(w);  plotc(w);    % show density
  title('Density by parzenm on labeled data')
  w = +a*parzenm; % Parzen density, classes combined, so just a single
                  % smoothing parameter optimizing overall density
  figure(4); scatterd(+a);% show scatterplot
  plotm(w);               % show density
  title('Density by parzenm on unlabeled data')

echo off
showfigs
