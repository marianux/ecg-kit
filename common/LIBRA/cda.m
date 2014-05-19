function result=cda(x,group,varargin)

%CDA performs linear and quadratic classical discriminant analysis 
%   on the data matrix x with known group structure. 
%
% Required input arguments:
%          x : training data set (matrix of size n by p).
%      group : column vector containing the group numbers of the training
%              set x. For the group numbers, any strict positive integer is
%              allowed assuming that the first group is the one with
%              the smallest group number.
%
% Optional input arguments:
%          method : String which indicates whether a 'linear' (default) or 'quadratic' 
%                   discriminant rule should be applied
%      misclassif : String which indicates how to estimate the probability of
%                   misclassification. It can be based on the  
%                   training data ('training'), a validation set ('valid'),
%                   or cross-validation ('cv'). Default is 'training'.
%  membershipprob : Vector which contains the membership probability of each
%                   group (sorted by increasing group number). If no priors are given, they are 
%                   estimated as the proportions of observations in the training set. 
%           valid : If misclassif was set to 'valid', this field should contain 
%                   the validation set (a matrix of size m by p).
%      groupvalid : If misclassif was set to 'valid', this field should contain the group numbers
%                   of the validation set (a column vector).
%      predictset : Contains a new data set (a matrix of size mp by p) from which the 
%                   class memberships are unknown and should be predicted. 
%           plots : If equal to 1, one figure is created with the training data and the
%                   tolerance ellipses for each group. This plot is
%                   only available for bivariate data sets. For technical reasons, a maximum 
%                   of 6 groups is allowed. Default is one.
%
% Options for advanced users (input comes from the program RSIMCA.m with option 'classic' = 1):
%
%    weightstrain : The weights for the training data. Corresponds to the flagtrain from RDA. (default = 1) 
%    weightsvalid : The weights for the validation data. Corresponds to the flagvalid from RDA. (default = 1)

% I/O: result=cda('plots',0,'misclassif','training','method','linear',...
%                  'membershipprob',proportions,'valid',y,'groupvalid',groupy);
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: out=cda(x,group,'method','linear')
%           out=cda(x,group,'plots',0)
%           out=cda(x,group,'valid',y,'groupvalid',groupy)
%
% The output is a structure containing the following fields:
%     result.assignedgroup : If there is a validation set, this vector contains the assigned group numbers
%                            for the observations of the validation set. Otherwise it contains the
%                            assigned group numbers of the original observations based on the discriminant rules.
%            result.scores : If there is a validation set, this column of size m contains the maximal discriminant
%                            scores for each observation from the validation set. Otherwise it is a columnvector of 
%                            size n containing the maximal discriminant scores of the training set.
%            result.method : String containing the method used to obtain
%                            the discriminant rules (either 'linear' or 'quadratic'). This
%                            is the same as the input argument method. 
%               result.cov : If method equals 'linear', this is a matrix containing the 
%                            estimated pooled covariance matrix. 
%                            If method equals 'quadratic', it is a cell array containing the covariances per group.
%            result.center : A vector in which the rows contain the estimated centers of the groups.
%                result.md : A vector of length n containing the mahalanobis distances of each observation from the training set to the
%                            center of its group.
%         result.flagtrain : Observations from the training set whose mahalanobis distance exceeds a certain cut-off value
%                            can be considered as outliers and receive a flag equal to zero. The regular observations
%                            receive a flag 1. (See also mcdcov.m)
%         result.flagvalid : Observations from the validation set whose mahalanobis distance (to the center of their group)
%                            exceeds a certain cut-off value can be considered as outliers and receive a
%                            flag equal to zero. The regular observations receive a flag 1. 
%                            If there is no validation set, this field is equal to zero.
%      result.grouppredict : If there is a prediction set, this vector contains the assigned group numbers
%                            for the observations of the prediction set. 
%       result.flagpredict : Observations from the new data set (predict) whose robust distance (to the center of their group)
%                            exceeds a certain cut-off value can be considered as overall outliers and receive a
%                            flag equal to zero. The regular observations receive a flag 1. 
%                            If there is no prediction set, this field is
%                            equal to zero.
%    result.membershipprob : A vector with the membership probabilities. 
%        result.misclassif : String containing the method used to estimate the misclassification probabilities
%                            (same as the input argument misclassif)
%  result.groupmisclasprob : A vector containing the misclassification
%                            probabilities for each group.
%    result.avemisclasprob : Overall probability of misclassification (weighted average of the misclassification
%                            probabilities over all groups).
%             result.class : 'CDA'
%           result.classic : Is equal to 0 since this analysis is a classical analysis. 
%                 result.x : The training data set (same as the input argument x).
%             result.group : The group numbers of the training set (same as the input argument group).
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Nele Smets on 02/02/2004
% Last Update: 01/07/2005
%

if nargin<2
    error('There are too few input arguments.')
end

%assigning default-values
[n,p]=size(x);
if size(group,1)~=1
    group=group';
end
if n ~= length(group)
    error('The number of observations is not the same as the length of the group vector!')
end
g=group;
countsorig=tabulate(g); %contingency table (outputmatrix with 3 colums: value - number - percentage)
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
    disp(['Warning: group(s) ', num2str(counts(counts(:,2)==0,1)'), 'are empty.']);
    empty=counts(counts(:,2)==0,:);
    counts=counts(counts(:,2)~=0,:);
else
    empty=[];
end
if any(counts(:,2)<5)%some groups have less than 5 observations
    error(['Group(s) ', num2str(counts(counts(:,2)<5,1)'), ' have less than 5 observations.']);    
end
proportions=zeros(size(counts,1),1);
y=0; %initial values of the validation data set and its groupsvector
groupy=0;
counter=1;
weightstrain = ones(1,n);
weightsvalid = 0;
default=struct('plots',1,'misclassif','training','method','linear','membershipprob',proportions,...
    'valid',y,'groupvalid',groupy,'weightstrain',weightstrain,'weightsvalid',weightsvalid,'predictset',[]);
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
if sum(prior) ~= 0
    epsilon=10^-4;
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
        error(['Group(s)' , num2str(empty(ismember(empty,countsvalid(:,1)))), 'was/were empty in the original dataset.'])
     end    
         if (length(options.weightsvalid) == 1) | (length(options.weightsvalid)~=size(validx,1))
        options.weightsvalid = ones(size(validx,1),1);
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
     countsvalid=countsvalidorig(countsvalidorig(:,2)~=0,:);
     if size(countsvalid,1)==1
        error('The validation set must contain more than one group!')
     elseif any(ismember(empty,countsvalid(:,1)))
        error(['Group(s)' , num2str(empty(ismember(empty,countsvalid(:,1)))), 'was/were empty in the original dataset.'])
     end
         if (length(options.weightsvalid) == 1) | (length(options.weightsvalid)~=size(validx,1))
        options.weightsvalid = ones(size(validx,1),1);
    end
end

%Discriminant rule based on the training set x
result1 = rawrule(x, g,prior, options.method);

%Discriminant rule based on reweighted results 
if strmatch(options.misclassif,'valid','exact') 
    result2 = rewrule(validx, result1);
    finalgroup = result2.class;
else 
    result2 = rewrule(x, result1);
    finalgroup = result2.class;
end

%Apply discriminant rule on validation set
switch options.misclassif
    case 'valid'
        [v,vi,vj]=unique(validgrouping);
         %Redefining the group number
        if any(v~= (1:length(v)))
            v=1:length(v);
            validgrouping=v(vj);
        elseif size(validgrouping,1)~=1
            validgrouping=validgrouping';
        end
        if any(countsvalidorig(:,2)==0)
            empty=setdiff(countsvalidorig(find(countsvalidorig(:,2)==0),1), countsorig(find(countsorig(:,2)==0)));
            disp(['Warning: the test group(s) ' , num2str(empty'), ' are empty']);
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
            if  ~isempty(intersect(i,v))
            misclas(i)=sum((validgrouping(options.weightsvalid == 1)==finalgroup(options.weightsvalid == 1)')...
                & (validgrouping(options.weightsvalid == 1)==repmat(lev(i),1,sum(options.weightsvalid))));
            ingroup(i) = sum((validgrouping(options.weightsvalid == 1) == repmat(lev(i), 1,sum(options.weightsvalid))));
            misclas(i) = (1 - (misclas(i)./ingroup(i)));
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
            misclas(i) = sum((g(options.weightstrain==1)==finalgroup(options.weightstrain==1)')...
                &(g(options.weightstrain==1)==repmat(lev(i),1,sum(options.weightstrain))));
            ingroup(i) = sum((g(options.weightstrain==1) == repmat(lev(i),1,sum(options.weightstrain))));
        end
        misclas = (1 - (misclas./ingroup));
        misclasprobpergroup = misclas;
        misclas = misclas.*result1.prior;
        misclasprob = sum(misclas);
        weightsvalid=0;%only available with validation set
    case 'cv'
        finalgroup=[];
        for i=1:length(x)
            xnew=removal(x,i,0);
            groupnew=removal(group,0,i);
            functie1res = rawrule(xnew, groupnew,prior, options.method);
            functie2res = rewrule(x(i, :), functie1res);
            finalgroup = [finalgroup; functie2res.class(1)];
        end
        for i=1:ng
            misclas(i) = sum((g(options.weightstrain==1) == finalgroup(options.weightstrain==1)')... 
                & (g(options.weightstrain==1) == repmat(lev(i),1,sum(options.weightstrain))));
            ingroup(i) = sum(g(options.weightstrain==1) == repmat(lev(i),1,sum(options.weightstrain)));
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

%Output structure
result=struct('assignedgroup',{finalgroup'},'scores',{result2.scores'},'method',{result2.method},...
    'cov',{result1.cov},'center',{result1.center},'md',{result1.mahal'}, 'flagtrain',{result1.flag'},...
    'flagvalid',{weightsvalid},'grouppredict',finalgrouppredict,'flagpredict',weightspredict','membershipprob',{result1.prior},...
    'misclassif',{options.misclassif},'groupmisclasprob',{misclasprobpergroup},...
    'avemisclasprob',{misclasprob},'class',{'CDA'},'x',{x},'group',{group});

if size(x,2)~=2
    result=rmfield(result,{'x','group'});
end

%Plotting the output
try
    if options.plots
        makeplot(result);
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end

%--------------------------------------------------------------------------
function result=rawrule(x, grouping,prior, method)

%calculate the discrimination rule based on the training set x.

[n,p]=size(x);	
g=grouping;
epsilon=10^-4;
[gun,gi,gj]=unique(g);

ng=length(gun);	
switch method
    case 'linear' %equal covariances supposed
        wsum=0;
        for j=1:length(gun)
            group.cov{j}=cov(x(g==gun(j),:));
            wsum=wsum+sum(g==gun(j))*group.cov{j};
            group.center(j,:)=mean(x(g==gun(j),:)); %center of all groups, matrix of ng x p 
        end 
        covar=wsum/n;
        for i=1:n
            zgeg(i,:)=x(i,:)-group.center(gj(i),:);
        end
        dist=zeros(n,1);
        for j=1:length(gun)
            dist(g==gun(j))=mahalanobis(x(g==gun(j),:),group.center(j,:),'invcov',inv(cov(zgeg)));
        end
        weights=zeros(n,1); 
        weights(dist <= chi2inv(0.975,p))=1;
        result.cov=covar; %over all group
        result.invcov=inv(covar);
        result.flag=weights;
        result.mahal=dist;
        result.center=group.center; %all groups
        result.method=method;
    case 'quadratic'
        for j=1:length(gun)
            group.cov{j}=cov(x(g==gun(j),:)); %covariance of group j
            group.invcov{j}=inv(group.cov{j});
            group.center(j,:)=mean(x(g==gun(j),:)); %center of all groups
        end  
        for i=1:n
            xdist(i)=mahalanobis(x(i,:), group.center(gj(i),:), 'invcov',group.invcov{gj(i)});
        end
        weights=zeros(n,1); 
        weights(xdist <= chi2inv(0.975,p))=1;
        result.cov=group.cov; %per group
        result.invcov=group.invcov;
        result.center = group.center; %all groups 
        result.mahal = xdist';
        result.flag = weights;
        result.method = method;
end        
if sum(prior) == 0
    counts = tabulate(g);
    if  ~any(counts(:,2))
        disp(['Warning: the group(s) ', num2str(counts(counts(:,2) == 0,1)'), 'contain only outliers']);
        counts=counts(counts(:,2)~=0,:);
    end
    result.prior = (counts(:,3)/100)';
else
    result.prior = prior;
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
if sum(prior)~=0
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
        
    
%--------------make sure the inputvariables are columnvectors!

function out=classification(x, center, covar,invcov, priorprob)
           
out=-0.5*log(abs(det(covar)))-0.5*(x - center)' * invcov *(x - center)+log(priorprob);
                
%-------------------
function out=linclassification(x, center, invcov, priorprob)

out=center'*invcov* x - 0.5*center'*invcov*center+log(priorprob);

