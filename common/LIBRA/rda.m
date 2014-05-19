function result=rda(x,group,varargin)

%RDA performs linear and quadratic robust discriminant analysis 
%   on the data matrix x with known group structure. It is based on the 
%   MCD estimator (see mcdcov.m), hence it has to be applied to 
%   low-dimensional data.
%
% The Robust Discriminant method is described in: 
%     Hubert, M., Van Driessen, K. (2004),
%     "Fast and Robust Discriminant Analysis," 
%     Computational Statistics and Data Analysis, 45, 301-320. 
%
% Required input arguments:
%              x : training data set (matrix of size n by p).
%          group : column vector containing the group numbers of the training
%                  set x. For the group numbers, any strict positive integer is
%                  allowed assuming that the first group is the one with the smallest group number.
%
% Optional input arguments:
%          alpha : (1-alpha) measures the fraction of outliers the MCD-algorithm should 
%                  resist. Any value between 0.5 and 1 may be specified. (default = 0.75)
%         method : String which indicates whether a 'linear' (default) or 'quadratic' 
%                  discriminant rule should be applied
%     misclassif : String which indicates how to estimate the probability of
%                  misclassification. It can be based on the  
%                  training data ('training'), a validation set ('valid'),
%                  or cross-validation ('cv'). Default is 'training'.
% membershipprob : Vector which contains the membership probability of each
%                  group (sorted by increasing group number). If no priors are given, they are estimated as the
%                  proportions of regular observations in the training set. 
%          valid : If misclassif was set to 'valid', this field should contain 
%                  the validation set (a matrix of size m by p).
%     groupvalid : If misclassif was set to 'valid', this field should contain the group numbers
%                  of the validation set (a column vector).
%     predictset : Contains a new data set (a matrix of size mp by p) from which the 
%                  class memberships are unknown and should be predicted.  
%          plots : If equal to 1, one figure is created with the training data and the
%                  MCD tolerance ellipses for each group. This plot is
%                  only available for bivariate data sets. For technical reasons, a maximum 
%                  of 6 groups is allowed. Default is one.
%        classic : If equal to one, classical linear or quadratic discriminant analysis will be performed
%                  (see also cda.m). (default = 0)
%        compare : If equal to one, the classical CDA analysis will be performed
%                  with the same weights and the same priors as the robust analysis
%                  has been performed. This is especially useful to compare the robust
%                  and classical result on the same data with the same priors. (default = 0)
%
% I/O: result=rda('alpha',0.5,'plots',0,'misclassif','training','method','linear',...
%                  'membershipprob',proportions,'valid',y,'groupvalid',groupy,'classic',0);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: out=rda(x,group,'method','linear')
%           out=rda(x,group,'plots',0)
%           out=rda(x,group,'valid',y,'groupvalid',groupy)
%
% The output is a structure containing the following fields:
%    
%    result.assignedgroup : If there is a validation set, this vector contains the assigned group numbers
%                           for the observations of the validation set. Otherwise it contains the
%                           assigned group numbers of the original observations based on the discriminant rules.
%           result.scores : If there is a validation set, this columnvector of size m contains the maximal discriminant
%                           scores for each observation from the validation set. Otherwise it is a columnvector of size n 
%                           containing the maximal discriminant scores of the training set.
%           result.method : String containing the method used to obtain the discriminant rules (either 'linear' or 'quadratic'). 
%                           This is the same as the input argument method. 
%              result.cov : If method equals 'linear', this is a matrix containing the estimated common covariance matrix.
%                           If method equals 'quadratic', it is a cell array containing the covariances per group.
%           result.center : A vector in which the rows contain the estimated centers of the groups.
%               result.rd : A vector of length n containing the robust distances of each observation from the training set 
%                           to the center of its group.
%        result.flagtrain : Observations from the training set whose robust distance exceeds a certain cut-off value
%                           can be considered as outliers and receive a flag equal to zero. 
%                           The regular observations  receive a flag 1. (See also mcdcov.m)
%        result.flagvalid : Observations from the validation set whose robust distance (to the center of their group)
%                           exceeds a certain cut-off value can be considered as outliers and receive a
%                           flag equal to zero. The regular observations receive a flag 1. 
%                           If there is no validation set, this field is equal to zero.
%     result.grouppredict : If there is a prediction set, this vector contains the assigned group numbers
%                           for the observations of the prediction set. 
%      result.flagpredict : Observations from the new data set (predict) whose robust distance (to the center of their group)
%                           exceeds a certain cut-off value can be considered as overall outliers and receive a
%                           flag equal to zero. The regular observations receive a flag 1. 
%                           If there is no prediction set, this field is equal to zero.
%   result.membershipprob : A vector with the membership probabilities. 
%       result.misclassif : String containing the method used to estimate the misclassification probabilities
%                           (same as the input argument misclassif)
% result.groupmisclasprob : A vector containing the misclassification probabilities for each group.
%   result.avemisclasprob : Overall probability of misclassification (weighted average of the misclassification
%                           probabilities over all groups).
%            result.class : 'RDA'
%          result.classic : If the input argument 'classic' is equal to one, this structure
%                           contains results of the classical discriminant analysis (see also cda.m). 
%          result.compare : If the input argument 'compare' is equal to one, this strucuture 
%                           contains results for the classical discriminant analysis with the same weights
%                           and priors as in the robust analysis.
%                result.x : The training data set (same as the input argument x).
%            result.group : The group numbers of the training set (same as the input argument group).
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets and Sabine Verboven on 01/03/2004
% Last Update: 01/07/2005
%

if nargin<2
    error('There are too few input arguments.')
end

% assigning default-values
[n,p]=size(x);
if size(group,1)~=1
    group=group';
end
if n ~= length(group)
    error('The number of observations is not the same as the length of the group vector!')
end
g=group;
countsorig=tabulate(g); %contingency table (outputmatrix with 3 colums): value - number - percentage 
[lev,levi,levj]=unique(g);
%Redefining the group number
if any(lev~= (1:length(lev)))
    lev=1:length(lev);
    g=lev(levj);
    counts=tabulate(g);
else
    counts=countsorig;
end

if ~all(counts(:,2)) %some groups have zero values, omit those groups
    disp(['Warning: group(s) ', num2str(counts(counts(:,2)==0,1)'), 'are empty']);
    empty=counts(counts(:,2)==0,:);
    counts=counts(counts(:,2)~=0,:);
else
    empty=[];
end

if any(counts(:,2)<5)%some groups have less than 5 observations
    error(['Group(s) ', num2str(counts(counts(:,2)<5,1)'), ' have less than 5 observations.']);    
end
proportions = zeros(size(counts,1),1); 
y=0; %initial values of the validation data set and its groupsvector
groupy=0;
counter=1;
default=struct('alpha',0.75,'plots',1,'misclassif','training','method','linear','membershipprob',proportions,...
    'valid',y,'groupvalid',groupy,'classic',0,'compare',0,'predictset',[]);
list=fieldnames(default);
options=default;
IN=length(list);
i=1;
%reading the user's input 
if nargin>2
    %
    %placing inputfields in array of strings
    %
    for j=1:nargin-2
        if rem(j,2)~=0
            chklist{i}=varargin{j};
            i=i+1;
        end
    end 
    %
    %Checking which default parameters have to be changed
    % and keep them in the structure 'options'.
    %
    while counter<=IN 
        index=strmatch(list(counter,:),chklist,'exact');
        if ~isempty(index) %in case of similarity
            for j=1:nargin-2 %searching the index of the accompanying field
                if rem(j,2)~=0 %fieldnames are placed on odd index
                    if strcmp(chklist{index},varargin{j})
                        I=j;
                    end
                end
            end
            options=setfield(options,chklist{index},varargin{I+1});
            index=[];
        end
        counter=counter+1;
    end
end

%Checking prior (>0 )
prior=options.membershipprob;
if size(prior,1)~=1
    prior=prior';
end
epsilon=10^-4;
if sum(prior) ~= 0
    if (any(prior < 0) | (abs(sum(prior)-1)) > epsilon)
        error('Invalid membership probabilities.')
    end
end
ng=length(proportions);
if length(prior)~=ng
    error('The number of membership probabilities is not the same as the number of groups.')
end


%%%%%%%%%%%%%%%%%%MAIN FUNCTION %%%%%%%%%%%%%%%%%%%%%
%Checking if a validation set is given
if strmatch(options.misclassif, 'valid','exact') 
    if options.valid==0
        error(['The misclassification error will be estimated through a validation set',...
                'but no validation set is given!'])
    else
        validx = options.valid;
        validgrouping = options.groupvalid;   
        if size(validx,1)~=length(validgrouping)
            error('The number of observations in the validation set is not the same as the length of its group vector!')
        end
        if size(validgrouping,1)~=1
            validgrouping = validgrouping';
        end
        countsvalidorig=tabulate(validgrouping);
        countsvalid=countsvalidorig(countsvalidorig(:,2)~=0,:);
        if size(countsvalid,1)==1
            error('The validation set must contain observations from more than one group!')
        elseif any(ismember(empty,countsvalid(:,1)))
            error(['Group(s) ' ,num2str(empty(ismember(empty,countsvalid(:,1)))), 'was/were empty in the original dataset.'])
        end
    end
elseif options.valid~=0
    validx = options.valid;
    validgrouping = options.groupvalid;
    if size(validx,1) ~= length(validgrouping)
        error('The number of observations in the validation set is not the same as the length of its group vector!')
    end
    if size(validgrouping,1)~=1
        validgrouping = validgrouping';
    end
    options.misclassif='valid';
    countsvalidorig=tabulate(validgrouping);
    countsvalid=countsvalidorig(countsvalidorig(:,2)~=0);
    if size(countsvalid,1)==1
        error('The validation set must contain more than one group!')
    elseif any(ismember(empty,countsvalid(:,1)))
        error(['Group(s) ' , num2str(empty(ismember(empty,countsvalid(:,1)))), ' was/were empty in the original dataset.'])
    end
end

%Discriminant rule based on the training set x
result1 = rawrule(x, g, prior, options.alpha, options.method);

%Apply discriminant rule on validation set
if strmatch(options.misclassif,'valid','exact') 
    result2 = rewrule(validx, result1);
    finalgroup = result2.class;
else 
    result2 = rewrule(x, result1);
    finalgroup = result2.class;
end

%Estimating the misclassification error 
switch options.misclassif
    case 'valid'   
        [v,vi,vj]=unique(validgrouping);
        %Redefining the group number
        if any(v~= (1:length(v)))
            v=1:length(v);
            validgrouping=v(vj);
        end
        if any(countsvalidorig(:,2)==0)
            empty=setdiff(countsvalidorig(find(countsvalidorig(:,2)==0),1), countsorig(find(countsorig(:,2)==0)));
            disp(['Warning: the test group(s) ' , num2str(empty), ' are empty']);
        else
            empty=[];
        end
        misclas=-ones(1,length(lev)); 
        for i=1:size(validx,1)
            if strmatch(options.method,'quadratic','exact')
                dist(i) = mahalanobis(validx(i,:), result1.center(vj(i),:),'invcov',result1.invcov{vj(i)});
            else
                dist(i) = mahalanobis(validx(i, :), result1.center(vj(i),:),'invcov',result1.invcov);
            end
        end
        weightsvalid=zeros(1,length(dist));
        weightsvalid(dist <= chi2inv(0.975, p))=1;
        for i=1:length(v)
            if ~isempty(intersect(i,v))
            misclas(i)=sum((validgrouping(weightsvalid==1)==finalgroup(weightsvalid==1)') & (validgrouping(weightsvalid==1)==repmat(lev(i),1,sum(weightsvalid))));
            ingroup(i) = sum((validgrouping(weightsvalid == 1) == repmat(lev(i),1, sum(weightsvalid))));
            misclas(i) = 1 - (misclas(i)./ingroup(i));
            end
        end
        if any(misclas==-1)
            misclas(misclas==-1)=0;
        end
        misclasprobpergroup=misclas;
        misclas=misclas.*result1.prior;
        misclasprob=sum(misclas);
    case 'training'
        for i=1:ng
            misclas(i) = sum((g(result1.weights==1)==finalgroup(result1.weights==1)')&(g(result1.weights==1)==repmat(lev(i),1,sum(result1.weights))));
            ingroup(i) = sum((g(result1.weights == 1) == repmat(lev(i),1,sum(result1.weights))));
        end
        misclas = (1 - (misclas./ingroup));
        misclasprobpergroup = misclas;
        misclas = misclas.*result1.prior;
        misclasprob = sum(misclas);
        weightsvalid=0;%only available with validation set
    case 'cv'
        finalgroup=[];
        for i=1:length(x)
            if (result1.weights(i) == 1)
                xnew=removal(x,i,0);
                groupnew=removal(group,0,i);
                functie1res = rawrule(xnew, groupnew,prior,options.alpha, options.method);
                functie2res = rewrule(x(i, :), functie1res);
                finalgroup = [finalgroup; functie2res.class(1)];
            end
        end
    for i=1:ng
        misclas(i) = sum((g(result1.weights == 1) == finalgroup') & (g(result1.weights == 1) == repmat(lev(i),1,sum(result1.weights))));
        ingroup(i) = sum(g(result1.weights == 1) == repmat(lev(i),1,sum(result1.weights)));
    end
    misclas = (1 - (misclas./ingroup));  
    misclasprobpergroup= misclas;
    misclas = misclas.* result1.prior;
    misclasprob = sum(misclas);
    weightsvalid=0; %only available with validation set
end

%classify the new observations (predict)
if ~isempty(options.predictset)
    resultpredict = rewrule(options.predictset, result1);
    finalgrouppredict = resultpredict.class;
    for i=1:size(options.predictset,1) 
        for j = 1:ng
            if strmatch(options.method,'quadratic','exact') 
                distpredict(i,j) = mahalanobis(options.predictset(i,:), result1.center(j,:),'invcov',result1.invcov{j});	
            else 
                distpredict(i,j) = mahalanobis(options.predictset(i, :), result1.center(j,:),'invcov',result1.invcov);
            end
        end
    end
    weightspredict = zeros(1,size(distpredict,1));
    weightspredict(min(distpredict,[],2) <= chi2inv(0.975, p))=1;
else
    finalgrouppredict = 0;
    weightspredict = 0;
end

if options.classic
    classicout=cda(x,g,'method',result2.method,'misclassif',options.misclassif,'membershipprob',options.membershipprob,'valid',options.valid,...
        'groupvalid',options.groupvalid,'plots',0,'predictset',options.predictset);
else
    classicout=0;
end

if options.compare
    compareout=cda(x,g,'method',result2.method,'misclassif',options.misclassif,'membershipprob',result1.prior,'valid',options.valid,...
        'groupvalid',options.groupvalid,'plots',0,'predictset',options.predictset,'weightstrain',result1.weights,'weightsvalid',weightsvalid);
else
    compareout=0;
end

%Output structure
result=struct('assignedgroup',{finalgroup'},'scores',{result2.scores'},'method',{result2.method},'cov',{result1.cov}, ...
    'center',{result1.center},'rd',{result1.dist'},'flagtrain',{result1.weights},...
    'flagvalid',weightsvalid,'grouppredict',finalgrouppredict,'flagpredict',weightspredict','membershipprob',{result1.prior},...
    'misclassif',{options.misclassif},'groupmisclasprob',{misclasprobpergroup},'avemisclasprob',{misclasprob},...
    'class',{'RDA'},'classic',{classicout},'compare',{compareout},'x',{x},'group',{group});

if size(x,2)~=2
    result=rmfield(result,{'x','group'});
end

%Plotting the output
try
    if options.plots & options.classic
        makeplot(result,'classic',1)
    elseif options.plots
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end
%--------------------------------------------------------------------------
function result=rawrule(x, g,prior, alfa, method)

%computes the discrimination rule based on the training set x.

[n,p]=size(x);	
epsilon=10^-4;
counts=tabulate(g); %contingency table (outputmatrix with 3 colums): value - number - percentage 
[lev,levi,levj]=unique(g);
if ~all(counts(:,2)) %some groups have zero values, omit those groups
    empty=counts(counts(:,2)==0,1);
else
    empty=[];
end

ng=size(counts,1);	
switch method
case 'linear' %equal covariances supposed
    [gun,gi,gj]=unique(g);
    for j=1:length(gun)
        group.mcd{j}=mcdcov(x(g==gun(j),:),'alpha',alfa,'plots',0); %covariance of group j
        group.center(j,:)=group.mcd{j}.center; %center of all groups, matrix of ng x p 
    end  
    for i=1:n
        zgeg(i,:)=x(i,:)-group.center(gj(i),:);
    end
    zmcd = mcdcov(zgeg,'alpha',alfa,'plots',0);
    zmcdcenter = zmcd.center;
    zmcdcov = zmcd.cov;
    zgeg = zgeg - repmat(zmcdcenter,length(zgeg),1);
    group.center = group.center + repmat(zmcdcenter,size(group.center,1),1);    
    dist=zeros(n,1);
    for j=1:length(gun)
        dist(g==gun(j))=mahalanobis(x(g==gun(j),:),group.center(j,:),'invcov',inv(zmcd.cov));
    end
    weights=zeros(n,1); 
    weights(dist <= chi2inv(0.975,p))=1;
    result.cov=zmcd.cov; %over all group
    result.invcov=inv(zmcd.cov);
    result.center=group.center; %all groups
    result.weights=weights';
    result.dist=dist;
    result.method=method;
case 'quadratic'
    [gun,gi,gj]=unique(g);
    xmcdweights=zeros(n,1);
    for j=1:length(gun)
        [group.mcd{j} raw{j}]=mcdcov(x(g==gun(j),:),'alpha',alfa,'plots',0);
        group.cov{j}=group.mcd{j}.cov; %covariance of group j
        group.invcov{j}=inv(group.cov{j});
        group.center(j,:)=group.mcd{j}.center; %center of all groups
        xmcdweights(g==gun(j))=raw{j}.wt;
    end  
    for i=1:n
        xdist(i)=mahalanobis(x(i,:), group.center(gj(i),:), 'invcov',group.invcov{gj(i)});
    end
    weights=xmcdweights;
    result.cov=group.cov; %per group
    result.invcov=group.invcov;
    result.center = group.center; %all groups 
    result.weights = xmcdweights';
    result.dist = xdist';
    result.method = method;      
end 

%Define the prior
if sum(prior) ~= 0
    result.prior = prior;
else    
    ngood=sum(weights);
    %regular points are kept
    ggood = g(weights==1);
    countsgood=tabulate(ggood);
    if empty
        for i=1:length(empty)
            countsgood(countsgood(:,1)==empty(i),:)= [];
        end
    end            
    if  ~any(countsgood(:,2)) 
        disp(['Warning: the group(s) ', num2str(countsgood(countsgood(:,2) == 0,1)'), 'contain only outliers']);
        countsgood=countsgood(countsgood(:,2)~=0,:);
    end
    result.prior = (countsgood(:,3)/100)';
end

%--------------------------------------------------------------------------
function result=rewrule(x, rawobject)

epsilon=10^-4;
center=rawobject.center;
covar=rawobject.cov;
invcov=rawobject.invcov;
prior=rawobject.prior;
method=rawobject.method;

if (length(prior) == 0 | length(prior) ~= size(center,1))
    error('invalid prior')
end
if sum(prior) ~= 0
    if (any(prior < 0) | (abs(sum(prior)-1)) > epsilon)
        error('invalid prior')
    end
end
ngroup=length(prior);
[n,p]=size(x);	
switch method
case 'linear'
    for j=1:ngroup 
        for i=1:n
            scores(i,j) = linclassification(x(i,:)', center(j,:)', invcov, prior(j));
        end
    end
    [maxs,maxsI] = max(scores,[],2); 
    for i=1:n
        maxscore(i,1) = scores(i,maxsI(i));
    end
    result.scores = maxscore;
    result.class = maxsI;
    result.method = method;   
case 'quadratic'
    for j=1:ngroup
        for i=1:n
            scores(i, j) = classification(x(i,:)', center(j,:)', covar{j}, invcov{j}, prior(j));
        end
    end
    [maxs,maxsI] = max(scores,[],2);
    for i=1:n
        maxscore(i,1) = scores(i,maxsI(i));
    end
    result.scores = maxscore;
    result.class = maxsI;
    result.method = method;
end


%--------------make sure the input variables are column vectors!

function out=classification(x, center, covar,invcov, priorprob)

out=-0.5*log(abs(det(covar)))-0.5*(x - center)' * invcov *(x - center)+log(priorprob);

%-------------------
function out=linclassification(x, center, invcov, priorprob)

out=center'*invcov*x - 0.5*center'*invcov*center+log(priorprob);

