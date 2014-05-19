%PREX_REGR PRTools regression example
%
% Show the regression functions that are available in Prtools in a 1-D
% example.
%
help prex_regr;

% Define the dataset parameters and generate the data:
n = 25; sig = 0.2;
a = gendatsinc(n,sig); % train
b = gendatsinc(n,sig); % test

% Train several regression functions:
w1 = a*linearr([],1e-9);
w2a = a*nu_svr([],'r',2,0.01,'epsilon',0.1);
w2b = a*svmr([],0.01,'r',2);
w3 = a*ridger([],10);
w4a = a*lassor([],10);
w4b = a*lassor([],1);
w5a = a*ksmoothr([],1);
w6 = a*pinvr;
w7 = a*plsr;
w8 = a*knnr([],1);
w9 = a*gpr([],proxm([],'r',1),0.1);

% Plot the functions in the scatterplot of the data:
figure(1); clf; hold on;
scatterr(a);
plotr(w1,'b-');
plotr(w2a,'r-');
plotr(w2b,'r--');
plotr(w3,'g-');
plotr(w4a,'k-');
plotr(w4b,'k--');
plotr(w5a,'m-');
plotr(w6,'y-');
plotr(w7,'c-');
plotr(w8,'b--');
plotr(w9,'k:');

% Show the MSE results:
fprintf('                           MSE\n');
fprintf('linear regression      : %f\n', b*w1*testr);
fprintf('nu-svm regression      : %f\n', b*w2a*testr);
fprintf('svm regression         : %f\n', b*w2b*testr);
fprintf('ridge regression       : %f\n', b*w3*testr);
fprintf('lasso regression (C=10): %f\n', b*w4a*testr);
fprintf('lasso regression (C=1) : %f\n', b*w4b*testr);
fprintf('smoother regression    : %f\n', b*w5a*testr);
fprintf('pseudo-inv regression  : %f\n', b*w6*testr);
fprintf('partial least squares  : %f\n', b*w7*testr);
fprintf('kNN regression         : %f\n', b*w8*testr);
fprintf('Gaussian Process regr. : %f\n', b*w9*testr);

