%PREX_COMBINING   PRTools example on classifier combining
%
% Presents the use of various fixed combiners for some 
% classifiers on the 'difficult data'.
%
help prex_combining

echo on
                % Generate 10-dimensional data
  A = gendatd([100,100],10);
                % Select the training set of 40 = 2x20 objects
                % and the test set of 160 = 2x80 objects
  [B,C] = gendat(A,0.2);

                % Define 5 untrained classifiers, (re)set their names
                % w1 is a linear discriminant (LDC) in the space reduced by PCA  
  w1 = klm([],0.95)*ldc;
  w1 = setname(w1,'klm - ldc');
                % w2 is an LDC on the best (1-NN leave-one-out error) 3 features 
  w2 = featself([],'NN',3)*ldc;
  w2 = setname(w2,'NN-FFS - ldc');
                % w3 is an LDC on the best (LDC leave-one-out error) 3 features 
  w3 = featself([],ldc,3)*ldc;
  w3 = setname(w3,'LDC-FFS - ldc');
                % w4 is an LDC 
  w4 = ldc;
  w4 = setname(w4,'ldc');
                % w5 is a 1-NN
  w5 = knnc([],1);
  w5 = setname(w5,'1-NN');

                % Store classifiers in a cell
  W = {w1,w2,w3,w4,w5};
                % Train them all
  V = B*W;
                % Test them all
  disp([newline 'Errors for individual classifiers'])
  testc(C,V);

                % Construct combined classifier
  VALL = [V{:}];
                % Define combiners
  WC = {prodc,meanc,medianc,maxc,minc,votec};
                % Combine (result is cell array of combined classifiers)
  VC = VALL * WC;
                % Test them all
  disp([newline 'Errors for combining rules'])
  testc(C,VC)
echo off
