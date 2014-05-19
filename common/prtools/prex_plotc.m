%PREX_PLOTC  PRTools example on the dataset scatter and classifier plot
help prex_plotc
echo on
              % Generate Higleyman data
A = gendath([100 100]); 
              % Split the data into the training and test sets
[C,D] = gendat(A,[20 20]);
              % Compute classifiers
w1 = ldc(C);        % linear
w2 = qdc(C);        % quadratic
w3 = parzenc(C);    % Parzen
w4 = dtc(C);        % decision tree
              % Compute and display errors
              % Store classifiers in a cell
W = {w1,w2,w3,w4};
              % Plot errors
disp(D*W*testc);    
              % Plot the data and classifiers
figure
              % Make a scatter-plot
scatterd(A);            
              % Plot classifiers
plotc({w1,w2,w3,w4});   
echo off
