function prediction = rbmPredict(model, testdata)
%Use RBM to predict discrete label for testdata

%INPUTS:
%model          ... is the model from rbmFit() consisting of W,b,c,Wc,cc
%testdata   ... binary, or in [0,1] interpreted as probabilities

%OUTPUTS:
%prediction ... the discrete labels for every class

% numclasses= size(m.Wc, 1);
% numcases= size(testdata, 1);
% F = zeros(numcases, numclasses);
% 
% %set every class bit in turn and find -free energy of the configuration
% for i=1:numclasses
%     X= zeros(numcases, numclasses);
%     X(:, i)=1;
%     F(:,i) = repmat(m.cc(i),numcases,1).*X(:,i)+ ...
%        sum(log(exp(testdata*m.W+ ...
%        X*m.Wc+repmat(m.b,numcases,1))+1),2);
% end
% 
% %take the max
% [q, predid]= max(F, [], 2);
% prediction=zeros(size(predid));
% for i=1:length(prediction) %convert back to users labels
%     prediction(i)= m.labels(predid(i));
% end

targetout = exp([testdata ones(size(testdata,1),1)]* [model.Wc model.cc']' );
% targetout = bsxfun( @rdivide, targetout, sum(targetout,2));
% [ ~, prediction] = max(targetout, [], 2);
prediction = bsxfun( @rdivide, targetout, sum(targetout,2));
