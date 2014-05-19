%PREX_COST PRTools example on cost matrices and rejection
%
% Prtools example code to show the use of cost matrices and how
% to introduce a reject class.

help prex_cost
echo on
                % Generate a three class problem
  randn('state',1);
  rand('state',1);
  n = 30;
  class_labels = char('apple','pear','banana');
  a = [gendatb([n,n]);  gendatgauss(n,[-2 6])];
  laba = genlab([n n n],class_labels);
  a = setlabels(a,laba);
                % Compute a simple ldc
  w = ldc(a);
                % Scatterplot and classifier
  figure;
  gridsize(30);
  scatterd(a,'legend');
  plotc(w);

                % Define a classifier with a new cost matrix,
					 % which puts a high cost on misclassifying
					 % pears to apples
  cost = [0.0  1.0  1.0;
          9.0  0.0  1.0;
		    1.0  1.0  0.0];
  wc = w*classc*costm([],cost,class_labels);
  plotc(wc,'b');

          % Define a classifier with a cost matrix where
					 % an outlier class is introduced. For this an
					 % extra column in the cost matrix has to be defined.
					 % Furthermore, the class labels have to be supplied
					 % to give the new class a name.
  cost = [0.0  1.0  1.0  0.2;
          9.0  0.0  1.0  0.2;
          1.0  1.0  0.0  0.2];
  class_labels = char('apple','pear','banana','reject');
  wr = w*classc*costm([],cost,class_labels);
  plotc(wr,'--')

echo off
disp(' ')
disp('   The black decision boundary shows the standard ldc classifier');
disp('   for this data. When the misclassification cost of a pear to an');
disp('   apple is increased, we obtain the blue classifier. When on top');
disp('   of that a rejection class is introduced, we get the blue dashed');
disp('   classifier. In that case, all objects between the dashed lines');
disp('   are rejected.');
fprintf('\n');
fprintf('  Cost of basic classifier  =  %4.2f\n',...
             a*w*testcost([],cost,class_labels));
fprintf('  Cost of cost classifier   =  %4.2f\n',...
             a*wc*testcost([],cost,class_labels));
fprintf('  Cost of reject classifier =  %4.2f\n',...
             a*wr*testcost([],cost,class_labels));
