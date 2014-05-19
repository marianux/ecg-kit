%PREX_SOFT Simple example of handling soft labels in PRTools
%
% Soft labels are implemented next to the 'crisp' and 'targets' labels.
% Like 'targets' labels they are stored in the target field of a dataset.
% Their values should be between 0 and 1. For every class a soft label
% values should be given. The density based classifiers can handle soft
% labels, interpreting them as class weights for every objects in the 
% density estimation. 
%
% The posterior probabilities found by classifying objects can be
% interpreted as soft labels. They, however, sum to one (over the classes),
% while this is not necessary for training and test objects.
%
% Note that the routine CLASSSIZES returns the sum of the soft labels over
% the dataset for every class separately. In contrast to crisp labels the
% sum over the classes of the output of CLASSSIZES is not necessarily
% equal to number of objects in the dataset.
%
% The routine SELDATA(A,N) returns the entire dataset in case of a soft
% labeled dataset A for every value of N and not just class N, as all
% objects may participate in all classes.

help prex_soft;
echo on

% Generate artificial soft labeled dataset using posteriors as soft labels
a = gendath([100 100]);

% retrieve a dataset with posteriors to be used for soft labels
labels = a*qdc(a)*classc;     

% create a new dataset with soft labels
s = prdataset(+a);              

% we just need the values of 'labels' 
s = setlabtype(s,'soft',+labels); 

% give the classes a name (optional, just to show how this is done)
s = setlablist(s,{'A','B'});  

% experiment: % generate train set and test set
[train_s,test_s] = gendat(s,0.5); 

% compute classifier that outputs posteriors
w_s = parzenc(train_s)*classc;    

% apply classifier on testdata
d_s = test_s*w_s;            

% result, by default for soft labeled data
% the 'soft' test type is used in testc
testc(d_s)                      
															
% compare with crisp labeling, convert train and test set to crisp labels
train_c = setlabtype(train_s,'crisp'); 
test_c  = setlabtype(test_s,'crisp'); 

 % compute classifier 
w_c = parzenc(train_c)*classc; 

% apply classifier on testdata
d_c = test_c*w_c;                   

% result, by default for crisp labeled data
% the 'crisp' test type is used in testc
testc(d_c)                     

echo off

