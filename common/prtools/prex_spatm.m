%PREX_SPATM  PRTools example on spatial smoothing of image classification
%
% Uses a multi-band image which is first segmented (classified) into 3
% classes based on the EM-clustering applied to a feature space 
% defined by the bands. The result is then smoothed by using the spatial
% information on the neighboring pixels.
%

help prex_spatm

delfigs
echo on
                  % Load the image
a = emim;
figure            % Show the multi-band image

show(a);        
                  % Extract a small training set
b = gendat(a,500);
                  % Use it for finding 3 clusters
[d,w] = emclust(b,nmc,3);
                  % Classify the entire image and show it
c = a*w;
figure; classim(c);
title('Original classification')
                  % Smooth the image and
                  % Combine the spectral and spatial classifier
                  % Show it
e = spatm(c)*maxc;

figure;
classim(e);
title('Smoothed classification')
echo off
showfigs
