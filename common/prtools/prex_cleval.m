%PREX_CLEVAL   PRTools example on learning curves
%
% Presents the learning curves for Highleyman's classes
%
help prex_cleval

delfigs
echo on
                % Set desired learning sizes
  learnsize = [3 5 10 15 20 30];
                % Generate Highleyman's classes
  A = gendath([100,100]); 
                % Define classifiers (untrained)
  W = {ldc,qdc,knnc([],1),treec};
                % Average error over 10 repetitions (it may take a while)
                % Test set is the complementary part of the training set
  E = cleval(A,W,learnsize,10);
                % Output E is a structure, specially designed for plotr
  plote(rmfield(E,'apperror'))   % plot without apparent error for clarity

echo off
