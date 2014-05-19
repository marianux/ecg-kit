function [model, top] = rbmBB(X, numhid, varargin)
%Learn RBM with Bernoulli hidden and visible units
%This is not meant to be applied to image data
%code by Andrej Karpathy
%based on implementation of Kevin Swersky and Ruslan Salakhutdinov

%INPUTS: 
%X              ... data. should be binary, or in [0,1] to be interpreted 
%               ... as probabilities
%numhid         ... number of hidden layers

%additional inputs (specified as name value pairs or in struct)
%method         ... CD or SML 
%eta            ... learning rate
%momentum       ... momentum for smoothness amd to prevent overfitting
%               ... NOTE: momentum is not recommended with SML
%maxepoch       ... # of epochs: each is a full pass through train data
%avglast        ... how many epochs before maxepoch to start averaging
%               ... before. Procedure suggested for faster convergence by
%               ... Kevin Swersky in his MSc thesis
%penalty        ... weight decay factor
%batchsize      ... The number of training instances per batch
%verbose        ... For printing progress
%anneal         ... Flag. If set true, the penalty is annealed linearly
%               ... through epochs to 10% of its original value

%OUTPUTS:
%model.type     ... Type of RBM (i.e. type of its visible and hidden units)
%model.W        ... The weights of the connections
%model.b        ... The biases of the hidden layer
%model.c        ... The biases of the visible layer
%model.top      ... The activity of the top layer, to be used when training
%               ... DBN's
%errors         ... The errors in reconstruction at every epoch

%Process options
%if args are just passed through in calls they become cells
if (isstruct(varargin)) 
    args= prepareArgs(varargin{1});
else
    args= prepareArgs(varargin);
end
[   method        ...
    eta           ...
    momentum      ...
    maxepoch      ...
    avglast       ...
    penalty       ...
    batchsize     ...
    verbose       ...
    anneal        ...
    ] = process_options(args    , ...
    'method'        ,  'CD'     , ...
    'eta'           ,  0.1      , ...
    'momentum'      ,  0.5      , ...
    'maxepoch'      ,  50       , ...
    'avglast'       ,  5        , ...
    'penalty'       , 2e-4      , ...
    'batchsize'     , 100       , ...
    'verbose'       , false     , ...
    'anneal'        , false);
avgstart = maxepoch - avglast;
oldpenalty= penalty;
[N,d]=size(X);

if (verbose) 
    fprintf('Preprocessing data...\n');
end

%Create batches
numcases=N;
numdims=d;
numbatches= ceil(N/batchsize);
groups= repmat(1:numbatches, 1, batchsize);
groups= groups(1:N);
perm=randperm(N);
groups = groups(perm);
for i=1:numbatches
    batchdata{i}= X(groups==i,:);
end

%train RBM
W = 0.1*randn(numdims,numhid);
c = zeros(1,numdims);
b = zeros(1,numhid);
ph = zeros(numcases,numhid);
nh = zeros(numcases,numhid);
phstates = zeros(numcases,numhid);
nhstates = zeros(numcases,numhid);
negdata = zeros(numcases,numdims);
negdatastates = zeros(numcases,numdims);
Winc  = zeros(numdims,numhid);
binc = zeros(1,numhid);
cinc = zeros(1,numdims);
Wavg = W;
bavg = b;
cavg = c;
t = 1;
errors=zeros(1,maxepoch);

for epoch = 1:maxepoch
    
	errsum=0;
    if (anneal)
        %apply linear weight penalty decay
        penalty= oldpenalty - 0.9*epoch/maxepoch*oldpenalty;
    end
    
    for batch = 1:numbatches
		[numcases numdims]=size(batchdata{batch});
		data = batchdata{batch};

        
%         %go up
% 		ph = logistic(data*W + repmat(b,numcases,1));
% 		phstates = ph > rand(numcases,numhid);
%         if (isequal(method,'SML'))
%             if (epoch == 1 && batch == 1)
%                 nhstates = phstates;
%             end
%         elseif (isequal(method,'CD'))
%             nhstates = phstates;
%         end

        %%%%%%%%% START POSITIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        poshidprobs = 1./(1 + exp(-data*W - repmat(b,numcases,1)));    
        posprods    = data' * poshidprobs;
        poshidact   = sum(poshidprobs);
        posvisact = sum(data);

        %%%%%%%%% END OF POSITIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        poshidstates = poshidprobs > rand(numcases,numhid);

%         %go down
% 		negdata = logistic(nhstates*W' + repmat(c,numcases,1));
% 		negdatastates = negdata > rand(numcases,numdims);
%         
%         %go up one more time
% 		nh = logistic(negdatastates*W + repmat(b,numcases,1));
% 		nhstates = nh > rand(numcases,numhid);
        
        %%%%%%%%% START NEGATIVE PHASE  %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        negdata = 1./(1 + exp(-poshidstates*W' - repmat(c,numcases,1)));
        neghidprobs = 1./(1 + exp(-negdata*W - repmat(b,numcases,1)));    
        negprods  = negdata'*neghidprobs;
        neghidact = sum(neghidprobs);
        negvisact = sum(negdata); 

        %%%%%%%%% END OF NEGATIVE PHASE %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%
        err= sum(sum( (data-negdata).^2 ));
        errsum = err + errsum;

        if epoch>5,
            this_momentum=1.8*momentum;
        else
            this_momentum=momentum;
        end;

%         %update weights and biases
%         dW = (data'*ph - negdatastates'*nh);
%         dc = sum(data) - sum(negdatastates);
%         db = sum(ph) - sum(nh);
% 		Winc = momentum*Winc + eta*(dW/numcases - penalty*W);
% 		binc = momentum*binc + eta*(db/numcases);
% 		cinc = momentum*cinc + eta*(dc/numcases);
% 		W = W + Winc;
% 		b = b + binc;
% 		c = c + cinc;
        
        %%%%%%%%% UPDATE WEIGHTS AND BIASES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        Winc = this_momentum*Winc + eta*( (posprods-negprods)/numcases - penalty*W);
        binc = this_momentum*binc + (eta/numcases)*(poshidact-neghidact);
        cinc = this_momentum*cinc + (eta/numcases)*(posvisact-negvisact);

        W = W + Winc;
        c = c + cinc;
        b = b + binc;

        %%%%%%%%%%%%%%%% END OF UPDATES %%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%%% 
        
        
        
%         if (epoch > avgstart)
%             %apply averaging
% 			Wavg = Wavg - (1/t)*(Wavg - W);
% 			cavg = cavg - (1/t)*(cavg - c);
% 			bavg = bavg - (1/t)*(bavg - b);
% 			t = t+1;
% 		else
			Wavg = W;
			bavg = b;
			cavg = c;
%         end
        
        %accumulate reconstruction error
        err= sum(sum( (data-negdata).^2 ));
		errsum = err + errsum;
        
    end
    
    errors(epoch)=errsum;
    if (verbose) 
        fprintf('Ended epoch %i/%i. Reconstruction error is %f\n', ...
            epoch, maxepoch, errsum);
    end
end

model.type= 'BB';
top = logistic(X*Wavg + repmat(bavg,N,1));
model.W= Wavg;
model.b= bavg;
model.c= cavg;
