%PREX_MATCHLAB   PRTools example on K-MEANS clustering and matching labels
%
% Illustrates the use of K-MEANS clustering and the match of labels.

help prex_matchlab
delfigs
echo on
     rand('state',5);   % Set up the random generator (used in K-MEANS)
     a = iris;
                        % Find clusters in the Iris dataset
     J1 = kmeans(a,3); 
                        % Find about the same clusters, but they are
     J2 = kmeans(a,3); 
                        % labeled differently due to a random initialization.
     confmat(J1,J2);   
                        % Match the labels. 'Best' rotation of label names 
     [J3,C] = matchlab(J1,J2); 
                        % since the confusion matrix is now almost diagonal.
     confmat(J1,J3);   
                        % Conversion from J2 to J3: J3 = C(J2,:);
     C                 
echo off
