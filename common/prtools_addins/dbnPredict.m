function prediction = dbnPredict(testdata, w)
%Use DBN to predict discrete label for testdata

%INPUTS:
%m          ... is the model from dbnFit()
%testdata   ... binary, or in [0,1] interpreted as probabilities

%OUTPUTS:
%prediction ... the discrete labels for every class

models = +w;
data = +testdata;

%map input all the way to the top
for i=1:length(models)
    data = rbmVtoH(models{i}, data);
end

%and predict on the last layer
data = rbmPredict(models{end}, data);

prediction = setdat(testdata, data, w);
