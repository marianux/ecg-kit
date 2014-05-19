function result = csimca(x,group,varargin);

%CSIMCA performs the SIMCA method. This is a classification
% method on a data matrix x with a known group structure. On each group a 
% robust PCA analysis is performed. Afterwards a classification
% rule is developped to determine the assignment of new observations. 
%
% Required input arguments:
%          x : training data set (matrix of size n by p).
%      group : column vector containing the group numbers of the training
%              set x. For the group numbers, any strict positive integer is
%              allowed.
%
% Optional input arguments:
%              k : Is a vector with size equal to the number of groups, or otherwise 0. Represents the number
%                  of components to be retained in each group. (default = 0.)
%         method : Indicates which classification rule is wanted. `1' results in rule (R1)
%                  based on the scaled orthogonal and score distances. `2' corresponds with
%                  (R2) based on the squared scaled orthogonal and score distances. Default is 2. 
%          gamma : Represents the value(s) used in the classification rule: weight gamma is given to the od's,
%                  weight (1-gamma) to the sd's. (default = 0.5).
%     misclassif : String which indicates how to estimate the probability of
%                  misclassification. It can be based on the training data ('training'), 
%                  a validation set ('valid'), or cross-validation ('cv'). Default is 'training'.
% membershipprob : Vector which contains the membership probability of each
%                  group (sorted by increasing group number). These values are used to
%                  obtain the total misclassification percentage. 
%          valid : If misclassif was set to 'valid', this field should contain 
%                  the validation set (a matrix of size m by p).
%     groupvalid : If misclassif was set to 'valid', this field should contain the group numbers
%                  of the validation set (a column vector).
%     predictset : Contains a new data set (a matrix of size mp by p) from which the 
%                  class memberships are unknown and should be predicted.  
%          plots : If equal to 1, one figure is created with the training data and the
%                  boundaries for each group. This plot is
%                  only available for trivariate (or smaller) data sets. For technical reasons, a maximum 
%                  of 6 groups is allowed. Default is one.
%       plotspca : If equal to one, a score diagnostic plot is
%                  drawn (default). If 'plots' is equal to zero, this plot is suppressed.
%                  See also makeplot.m
%          labsd : The 'labsd' observations with largest score distance are
%                  labeled on the diagnostic plot. (default = 3)
%          labod : The 'labod' observations with largest orthogonal distance are
%                  labeled on the diagnostic plot. default = 3) 
%
% Options for advanced users (input comes from the program RSIMCA.m with option 'classic' = 1):
%
%   weightstrain : The weights for the training data. Corresponds to the flagtrain from RSIMCA. (default = 1) 
%   weightsvalid : The weights for the validation data. Corresponds to the flagvalid from RSIMCA. (default = 1)
%
% I/O: result=csimca(x,group,'method',1,'misclassif','training',...
%                  'membershipprob',proportions,'valid',y,'groupvalid',groupy,'plots',0);
%
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: out=csimca(x,group,'method','1')
%           out=csimca(x,group,'plots',0)
%           out=csimca(x,group,'valid',y,'groupvalid',groupy)
%
% The output is a structure containing the following fields:
%        result.assignedgroup : If there is a validation set, this vector contains the assigned group numbers
%                               for the observations of the validation set. Otherwise it contains the
%                               assigned group numbers of the original observations based on the discriminant rules.
%                result.pca   : A cell containing the results of the different PCA analysis on the training sets.
%               result.method : String containing the method used to obtain
%                               the discriminant rules (either 1 for 'R1' or 2 for 'R2'). This
%                               corresponds to the input argument method. 
%            result.flagtrain : Observations from the training set whose score distance and/or orthogonal distance
%                               exceeds a certain cut-off value can be considered as outliers and receive a flag equal 
%                               to zero. The regular observations receive a flag 1. (See also robpca.m)
%            result.flagvalid : Observations from the validation set whose score distance and/or orthogonal distance 
%                               exceeds a certain cut-off value can be considered as outliers and receive a
%                               flag equal to zero. The regular observations receive a flag 1. 
%                               If there is no validation set, this field is equal to zero.
%         result.grouppredict : If there is a prediction set, this vector contains the assigned group numbers
%                               for the observations of the prediction set. 
%          result.flagpredict : Observations from the new data set (predict) whose robust distance (to the center of their group)
%                               exceeds a certain cut-off value can be considered as overall outliers and receive a
%                               flag equal to zero. The regular observations receive a flag 1. 
%                               If there is no prediction set, this field is
%                               equal to zero.
%       result.membershipprob : A vector with the membership probabilities.  If no priors are given, they are estimated 
%                               as the proportions of observations in the training set.
%           result.misclassif : String containing the method used to estimate the misclassification probabilities
%                               (same as the input argument misclassif)
%     result.groupmisclasprob : A vector containing the misclassification probabilities for each group.
%       result.avemisclasprob : Overall probability of misclassification (weighted average of the misclassification
%                               probabilities over all groups).
%                result.class : 'CSIMCA'
%                    result.x : The training data set (same as the input argument x) (only in output when p<=3).
%                result.group : The group numbers of the training set (same as the input argument group) (only in output when p<=3).
%
% This function is part of LIBRA: the Matlab Library for Robust Analysis,
% available at: 
%              http://wis.kuleuven.be/stat/robust.html
%
% Written by Karlien Vanden Branden  
% Last Update: 05/07/2005
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
counts=tabulate(g); %contingency table (outputmatrix with 3 colums): value - number - percentage 
[lev,levi,levj]=unique(g);
if ~all(counts(:,2)) %some groups have zero values, omit those groups
    disp(['Warning: group(s) ', num2str(counts(counts(:,2)==0,1)'), 'are empty']);
    empty=counts(counts(:,2)==0,:);
    counts=counts(counts(:,2)~=0,:);
else
    empty=[];
end
ng=size(counts,1);
proportions = zeros(ng,1);
y=0; %initial values of the validation data set and its groupsvector
groupy=0;
labsd = 3;
labod = 3;
counter=1;
gamma = 0.5;
k = zeros(ng,1);
weightstrain = ones(1,n);
weightsvalid = 0;
default=struct('k',k,'method',2,'gamma',0.5,'misclassif','training',...
    'membershipprob',proportions,'valid',y,'groupvalid',groupy,'plots',1,'plotspca',0,'labsd',labsd,...
    'labod',labod,'weightstrain',weightstrain,'weightsvalid',0,'predictset',[]);
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

%Checking gamma
gamma = options.gamma;
if gamma >1 | gamma <0
    error('An inappropriate number for gamma is given. A correct value lies between 0 and 1.');
end

%Checking prior (>0 )
prior=options.membershipprob;
if size(prior,1)~=1
    prior=prior';
end
epsilon=10^-4;
if sum(prior) ~=0 & (any(prior < 0) | (abs(sum(prior)-1)) > epsilon)
   error('Invalid membership probabilities.')
end
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
        validgroup = options.groupvalid;
        if size(validx,1)~=length(validgroup)
            error('The number of observations in the validation set is not the same as the length of its group vector!')
        end
        if size(validgroup,1)~=1
            validgroup = validgroup';
        end
        countsvalid=tabulate(validgroup);
        countsvalid=countsvalid(countsvalid(:,2)~=0,:);
        if size(countsvalid,1)==1
            error('The validation set must contain observations from more than one group!')
        elseif any(ismember(empty,countsvalid(:,1)))
            error(['Group(s) ' ,num2str(empty(ismember(empty,countsvalid(:,1)))), 'was/were empty in the original dataset.'])
        end
    end
    if (length(options.weightsvalid) == 1) | (length(options.weightsvalid)~=size(validx,1))
        options.weightsvalid = ones(size(validx,1),1);
    end
elseif options.valid~=0
    validx = options.valid;
    validgroup = options.groupvalid;
    if size(validx,1) ~= length(validgroup)
        error('The number of observations in the validation set is not the same as the length of its group vector!')
    end
    if size(validgroup,1)~=1
        validgroup = validgroup';
    end
    options.misclassif='valid';
    countsvalid=tabulate(validgroup);
    countsvalid=countsvalid(countsvalid(:,2)~=0,:);
    if size(countsvalid,1)==1
        error('The validation set must contain more than one group!')
    elseif any(ismember(empty,countsvalid(:,1)))
        error(['Group(s) ' , num2str(empty(ismember(empty,countsvalid(:,1)))), ' was/were empty in the original dataset.'])
    end
    if (length(options.weightsvalid) == 1) | (length(options.weightsvalid)~=size(validx,1))
        options.weightsvalid = ones(size(validx,1),1);
    end
end

model.counts = counts(:,2);
model.x = x;
model.group = group;

labsd = floor(max(0,min(options.labsd,n)));
labod = floor(max(0,min(options.labod,n)));

%CSIMCA:
%PRINCIPAL COMPONENT ANALYSIS
%Perform PCA on each group separately:
%   a) if k is not given: decide on the optimal number of components using CV.
%   b) if k is given: perform cpca with the optimal number of components.

for iClass = 1:ng
    indexgroup = find(g==iClass);
    groupi = x(indexgroup,:); 
    if options.k(iClass) == 0
        disp(['A scree plot is drawn for group ',num2str(iClass),'.'])
    end
    model.result{iClass} = cpca(groupi,'k',options.k(iClass),'plots',options.plotspca,...
        'labsd',labsd,'labod',labod);
    model.flag(1,indexgroup) = model.result{iClass}.flag.all;
end

%CLASSIFICATION
%Discriminant rule based on the training set x
[odsc,sdsc] = testmodel(model,x);
finalgrouptrain = assigngroup(odsc,sdsc,options.method,ng,gamma); 
    
if sum(prior) == 0
    result1.prior = counts(:,3)'./100;
else
    result1.prior = prior;
end

%Compute scaled orthogonal and scaled score distances for the validation set
if strmatch(options.misclassif,'valid','exact') 
    [odsc,sdsc] = testmodel(model,validx);
    finalgroup = assigngroup(odsc,sdsc,options.method,ng,gamma);
elseif strmatch(options.misclassif,'cv','exact')  %use cv
    model.k = k;
    [odsc,sdsc] = leave1out(model);
    finalgroup = assigngroup(odsc,sdsc,options.method,ng,gamma);
end

switch options.misclassif
case 'valid'
    [v,vi,vj]=unique(validgroup);
    odscgroup = [];
    sdscgroup = [];
    for iClass = 1:ng
        indexgroup = find(validgroup == iClass);
        odscgroup = [odscgroup;odsc(indexgroup,iClass)];
        sdscgroup = [sdscgroup;sdsc(indexgroup,iClass)];
    end 
    weightsvalid=zeros(length(odscgroup),1); 
    weightsvalid(((odscgroup <= 1) & (sdscgroup <= 1)))=1;
    for igamma = 1:length(gamma)
        for iClass=1:ng 
            misclas(iClass)=sum((validgroup(options.weightsvalid==1)==finalgroup(options.weightsvalid==1,igamma)') & ...
                (validgroup(options.weightsvalid==1)==repmat(v(iClass),1,sum(options.weightsvalid))));
            ingroup(iClass) = sum((validgroup(options.weightsvalid == 1) == repmat(v(iClass),1,sum(options.weightsvalid))));
        end
        misclas = (1 - (misclas./ingroup));
        misclasprobpergroup(igamma,:)=misclas;
        misclas=misclas.*result1.prior;
        misclasprob(igamma)=sum(misclas);    
    end
case 'training'
    for igamma = 1:length(gamma)
        for iClass = 1:ng
            result1.misclas(iClass) = sum((group(options.weightstrain==1)==finalgrouptrain(options.weightstrain==1,igamma)')& ...
                (group(options.weightstrain==1)==repmat(lev(iClass),1,sum(options.weightstrain))));
            result1.ingroup(iClass) = sum((group(options.weightstrain == 1) == repmat(lev(iClass),1,sum(options.weightstrain))));
        end
        misclas = (1 - (result1.misclas./result1.ingroup));
        misclasprobpergroup(igamma,:) = misclas;
        misclas = misclas.*result1.prior;
        misclasprob(igamma) = sum(misclas);
    end
    weightsvalid=0;%only available with validation set
    finalgroup = finalgrouptrain;
case 'cv' 
    for igamma = 1:length(gamma)
        for iClass=1:ng 
            misclas(iClass)=sum((group(options.weightstrain==1)==finalgroup(options.weightstrain==1,igamma)') & ...
                (group(options.weightstrain==1)==repmat(lev(iClass),1,sum(options.weightstrain))));
            ingroup(iClass) = sum((group(options.weightstrain == 1) == repmat(lev(iClass),1,sum(options.weightstrain))));
        end
        misclas = (1 - (misclas./ingroup));
        misclasprobpergroup(igamma,:)=misclas;
        misclas=misclas.*result1.prior;
        misclasprob(igamma)=sum(misclas);    
    end
    weightsvalid=0; %only available with validation set
end

if ~isempty(options.predictset)
    [odscpredict,sdscpredict] = testmodel(model,options.predictset);
    finalgrouppredict = assigngroup(odscpredict,sdscpredict,options.method,ng,gamma)';
    weightspredict = max((odscpredict <= 1) & (sdscpredict <= 1),[],2)';
else
    finalgrouppredict = 0;
    weightspredict = 0;
end   

%Output structure
result = struct('assignedgroup',{finalgroup'},'pca',{model.result},'method',options.method,...
    'flagtrain',{model.flag},'flagvalid',{weightsvalid'},'grouppredict',finalgrouppredict,'flagpredict',weightspredict,...
    'membershipprob',{result1.prior},'misclassif',{options.misclassif},'groupmisclasprob',{misclasprobpergroup},...
    'avemisclasprob',{misclasprob},'class',{'CSIMCA'},'x',x,'group',group);

if size(x,2)>3
    result=rmfield(result,{'x','group'});
end

%Plots:
try
    if options.plots
        makeplot(result)
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end

%---------------------------
%Leave-One-Out procedure 

function [odsc,sdsc] = leave1out(model)

nClass = length(model.result);

 for iClass = 1:nClass
    index = 1;
    indexgroup = find(model.group == iClass);
    teller_if_lus = 0;
    groupi = model.x(indexgroup,:);
    for i = 1:model.counts(iClass)
        groupia = removal(groupi,i,0);
        GRes = model.result;
        GRes{iClass} = cpca(groupia,'k',model.result{iClass}.k,'plots',0);
        %Calculate for each class the sd and the od for the observation that was left out:
        for jClass = 1:nClass
            dataicentered = model.x(indexgroup(index),:)-GRes{jClass}.M;
            scorei = dataicentered*GRes{jClass}.P;
            dataitilde = scorei*GRes{jClass}.P';
            sd(indexgroup(index),jClass) = sqrt(scorei*(diag(1./GRes{jClass}.L))*scorei');
            od(indexgroup(index),jClass) = norm(dataicentered-dataitilde);
            if GRes{jClass}.cutoff.od ~= 0  
                odsc(indexgroup(index),jClass) = od(indexgroup(index),jClass)/GRes{jClass}.cutoff.od;
            else
                odsc(indexgroup(index),jClass) = 0;
            end
            sdsc(indexgroup(index),jClass) = sd(indexgroup(index),jClass)/GRes{jClass}.cutoff.sd;
        end
        index = index + 1;
    end
end

%--------------------------------
function [odsc,sdsc] = testmodel(model,validx);

%Apply the given model on the test data to obtain different 
%orthogonal distances and score distances.

nClass = length(model.result);
n = size(validx,1);

for jClass = 1:nClass
    for index = 1:n
        out{jClass} = model.result{jClass};
        dataicentered = validx(index,:)-out{jClass}.M;
        scorei = dataicentered*out{jClass}.P;
        dataitilde = scorei*out{jClass}.P';
        sd(index,jClass) = sqrt(scorei*(diag(1./out{jClass}.L))*scorei');
        od(index,jClass) = norm(dataicentered-dataitilde);
        if out{jClass}.cutoff.od ~= 0  
            odsc(index,jClass) = od(index,jClass)/out{jClass}.cutoff.od;
        else
            odsc(index,jClass) = 0;
        end
        sdsc(index,jClass) = sd(index,jClass)/out{jClass}.cutoff.sd; 
    end
end

%-------------------------
function result = assigngroup(odsc,sdsc,method,nClass,gamma);

%Obtain the group assignments for given od's and sd's.

if method == 1 
    sd = sdsc;
    od = odsc;
elseif method == 2 
    sd = sdsc.^2;
    od = odsc.^2;
end

for igamma = 1:length(gamma)
    tdist = gamma(igamma).*od + (1-gamma(igamma)).*sd;
    for i = 1:size(od,1)
        result(i,igamma) = find(tdist(i,:) == min(tdist(i,:)));
    end
end


