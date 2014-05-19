function result = rsimca(x,group,varargin)

%RSIMCA performs a robust version of the SIMCA method. This is a classification
% method on a data matrix x with a known group structure. On each group a 
% robust PCA analysis (ROBPCA) is performed. Afterwards a classification
% rule is developped to determine the assignment of new observations. Since RSIMCA
% depends on the ROBPCA method (see robpca.m) it is able to deal with high-dimensional data.
%
% The RSIMCA method is described in:
%    K. Vanden Branden and M. Hubert (2005), 
%    Robust classification in high dimensions based on the SIMCA method,
%    Chemometrics and Intelligent Laboratory Systems, 79, 10-21.
%
% Required input arguments:
%          x : training data set (matrix of size n by p).
%      group : column vector containing the group numbers of the training
%              set x. For the group numbers, any strict positive integer is
%              allowed.
%
% Optional input arguments:
%          alpha : (1-alpha) measures the fraction of outliers (in each group) the algorithm should 
%                  be able to resist. Any value between 0.5 and 1 may be specified. (default = 0.75) 
%                  If k groups are present in the data, per group a different value for alpha may be specified 
%                  (default = alpha = [0.75, ... , 0.75]).
%              k : Is a vector with size equal to the number of groups, or otherwise 0. It represents the number
%                  of components to be retained in each group. (default = 0).
%           kmax : Maximal number of principal components to compute (default = 10).
%                  If k is provided, kmax does not need to be specified, unless k is larger
%                  than 10.    
%          scree : If equal to one, a scree plot is drawn for each group. (default = 0).
%          press : If equal to one, a plot of robust press-values is drawn for each group.
%                  If k is given as input, the default value is 0, else the default value is one.
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
%                  obtain the total misclassification percentage. If no priors are given, they are 
%                  estimated as the proportions of regular observations in the training set. 
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
%    plotsrobpca : If equal to one, a robust score diagnostic plot is
%                  drawn (default). If the input argument 'classic' is equal to one, 
%                  the classical plots are drawn as well.
%                  If 'plots' is equal to zero, this plot is suppressed.
%                  See also makeplot.m
%          labsd : The 'labsd' observations with largest score distance are
%                  labeled on the diagnostic plot. (default = 3)
%          labod : The 'labod' observations with largest orthogonal distance are
%                  labeled on the diagnostic plot. default = 3) 
%            mcd : If equal to one: when the number of variables is sufficiently small,
%                  the loadings are computed as the eigenvectors of the MCD covariance matrix, 
%                  hence the function 'mcdcov.m' is automatically called. The number of 
%                  principal components is then taken as k = rank(x). (default)
%                  If equal to zero, the robpca algorithm is always applied.
%        classic : If equal to one, the classical SIMCA analysis will be performed
%                  (see also csimca.m). (default = 0)
%        compare : If equal to one, the classical SIMCA analysis will be performed
%                  with the same weights and the same priors as the robust analysis
%                  has been performed. This is especially useful to compare the robust
%                  and classical result on the same data with the same priors. (default = 0)
%         
%
% I/O: result=rsimca(x,group,'alpha',0.5,'method',1,'misclassif','training',...
%                  'membershipprob',proportions,'valid',y,'groupvalid',groupy,'plots',0,'classic',0);
%
%  The user should only give the input arguments that have to change their default value.
%  The name of the input arguments needs to be followed by their value.
%  The order of the input arguments is of no importance.
%
% Examples: out=rsimca(x,group,'method','1')
%           out=rsimca(x,group,'plots',0)
%           out=rsimca(x,group,'valid',y,'groupvalid',groupy)
%
% The output is a structure containing the following fields:
%        result.assignedgroup : If there is a validation set, this vector contains the assigned group numbers
%                               for the observations of the validation set. Otherwise it contains the
%                               assigned group numbers of the original observations based on the discriminant rules.
%                  result.pca : A cell containing the results of the different ROBPCA analysis on the training sets.
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
%       result.membershipprob : A vector with the membership probabilities. If no priors are given, they are estimated 
%                               as the proportions of regular observations in the training set.
%           result.misclassif : String containing the method used to estimate the misclassification probabilities
%                               (same as the input argument misclassif).
%     result.groupmisclasprob : A vector containing the misclassification probabilities for each group.
%       result.avemisclasprob : Overall probability of misclassification (weighted average of the misclassification
%                               probabilities over all groups).
%                result.class : 'RSIMCA'
%              result.classic : If the input argument 'classic' is equal to one, this structure
%                               contains results of the classical SIMCA analysis. 
%              result.compare : If the input argument 'compare' is equal to one, this strucuture 
%                               contains results for the classical SIMCA analysis with the same weights
%                               and priors as in the robust analysis.
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
plotsrobpca = 0;
press = 1;
scree = 0;
counter=1;
gamma = 0.5;
k = zeros(ng,1);
for iClass = 1:ng
    r(iClass,1) = rank(x(find(group == iClass),:));
    kmax(iClass,1) = min([10,floor(counts(iClass,2)/2),r(iClass)]);
end
default=struct('alpha',0.75,'k',k,'kmax',kmax,'scree',scree,'press',press,'method',2,...
    'gamma',0.5,'misclassif','training','membershipprob',proportions,'valid',y,...
    'groupvalid',groupy,'plots',1,'plotsrobpca',plotsrobpca,'labsd',labsd,'labod',labod,...
    'mcd',0,'classic',0,'compare',0,'predictset',[]);
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
if any(gamma>1) | any(gamma<0)
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
end

if length(options.alpha) == 1
    alfa = ones(ng,1)*options.alpha;
elseif length(options.alpha) ~= ng
    error('The length of alpha does not correspond with the number of groups.');
else
    alfa = options.alpha;
end
model.counts = counts(:,2);
model.x = x;
model.group = group;

%Checking input variable k & kmax:
if sum(options.kmax>0)~=ng
    mess=sprintf(['Attention (rsimca.m): The value for kmax is incorrect, kmax = ',num2str(options.kmax'),...
            '\n is smaller than 0. kmax is set to ',num2str(kmax')]);
    disp(mess)
end
if sum(k<=options.kmax) ~= ng
    error('The value for k is set too large.');
end
if sum(options.k ~= 0) == ng 
    press = 0;
    scree = 0;
else
    press = options.press;
    scree = options.scree;
end

labsd = floor(max(0,min(options.labsd,n)));
labod = floor(max(0,min(options.labod,n)));

%RSIMCA:
%PRINCIPAL COMPONENT ANALYSIS
%Perform ROBPCA on each group separately:
%   a) if k is not given: decide on the optimal number of components using CV.
%   b) if k is given: perform robpca with the optimal number of components.

for iClass = 1:ng
    indexgroup = find(g==iClass);
    groupi = x(indexgroup,:); 
    if press == 1 & scree == 1
        disp(['A press curve and a scree plot are drawn for group ',num2str(iClass),'.'])
    elseif press == 1 & scree == 0
        disp(['A press curve is drawn for group ',num2str(iClass),'.'])
    elseif press == 0 & scree == 1
        disp(['A scree plot is drawn for group ',num2str(iClass),'.'])
    end
    model.result{iClass} = robpca(groupi,'k',options.k(iClass),'kmax',options.kmax(iClass),'plots',options.plotsrobpca,...
    'alpha',alfa(iClass),'scree',scree,'press',press,'labsd',labsd,'labod',labod,'mcd',options.mcd);
 
    model.flag(1,indexgroup) = model.result{iClass}.flag.all;
    model.k(iClass) = model.result{iClass}.k;
end

%CLASSIFICATION
%Discriminant rule based on the training set x
[odsc,sdsc] = testmodel(model,x,g);
finalgrouptrain = assigngroup(odsc,sdsc,options.method,ng,gamma); 
result1.weights = model.flag;
if sum(prior) == 0
    for iClass=1:ng
        result1.ingroup(iClass) = sum((group(result1.weights == 1) == repmat(lev(iClass),1,sum(result1.weights))));
    end
    result1.prior = result1.ingroup./sum(result1.ingroup)';
else
    result1.prior = prior;
end

%Compute scaled orthogonal and scaled score distances for the validation set
if strmatch(options.misclassif,'valid','exact') 
    [odsc,sdsc] = testmodel(model,validx);
    finalgroup = assigngroup(odsc,sdsc,options.method,ng,gamma);
elseif strmatch(options.misclassif,'cv','exact')  %use cv
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
    weightsvalid=zeros(1,length(odscgroup)); 
    weightsvalid(((odscgroup <= 1) & (sdscgroup <= 1)))=1;
    for igamma = 1:length(gamma)
        for iClass=1:ng 
            misclas(iClass)=sum((validgroup(weightsvalid==1)==finalgroup(weightsvalid==1,igamma)') & ...
                (validgroup(weightsvalid==1)==repmat(v(iClass),1,sum(weightsvalid))));
            ingroup(iClass) = sum((validgroup(weightsvalid == 1) == repmat(v(iClass),1, sum(weightsvalid))));
        end
        misclas = (1 - (misclas./ingroup));
        misclasprobpergroup(igamma,:)=misclas;
        misclas=misclas.*result1.prior;
        misclasprob(igamma)=sum(misclas);    
    end
case 'training'
    for igamma = 1:length(gamma)
        for iClass = 1:ng
            result1.misclas(iClass) = sum((group(result1.weights==1)==finalgrouptrain(result1.weights==1,igamma)')& ...
                (group(result1.weights==1)==repmat(lev(iClass),1,sum(result1.weights))));
            result1.ingroup(iClass) = sum((group(result1.weights == 1) == repmat(lev(iClass),1,sum(result1.weights))));
        end
        misclas = (1 - (result1.misclas./result1.ingroup));
        misclasprobpergroup(igamma,:) = misclas;
        misclas = misclas.*result1.prior;
        misclasprob(igamma) = sum(misclas);
        weightsvalid=0;%only available with validation set
    end
    finalgroup = finalgrouptrain;
case 'cv' 
    for igamma = 1:length(gamma)
        for iClass=1:ng 
            misclas(iClass)=sum((group(result1.weights==1)==finalgroup(result1.weights==1,igamma)') & ...
                (group(result1.weights==1)==repmat(lev(iClass),1,sum(result1.weights))));
            ingroup(iClass) = sum((group(result1.weights == 1) == repmat(lev(iClass),1,sum(result1.weights))));
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
    weightspredict = (max((odscpredict <= 1) & (sdscpredict <= 1),[],2))';
else
    finalgrouppredict = 0;
    weightspredict = 0;
end   
    
if options.classic
    classicout=csimca(x,g,'k',model.k,'scree',scree,'press',press,'method',options.method,'gamma',gamma,...
        'misclassif',options.misclassif,'membershipprob',prior,'valid',options.valid,...
        'groupvalid',options.groupvalid,'plots',0,'plotspca',0,'labsd',labsd,'labod',labod,...
        'predictset',options.predictset);
else
     classicout=0;
end

if options.compare
   compareout=csimca(x,g,'k',model.k,'scree',scree,'press',press,'method',options.method,'gamma',gamma,...
        'misclassif',options.misclassif,'membershipprob',result1.prior,'valid',options.valid,...
        'groupvalid',options.groupvalid,'plots',0,'plotspca',0,'labsd',labsd,'labod',labod,...
        'weightstrain',result1.weights,'weightsvalid',weightsvalid,'predictset',options.predictset);
else
    compareout = 0;
end

%Output structure
result = struct('assignedgroup',{finalgroup'},'pca',{model.result},'method',options.method,...
    'flagtrain',{model.flag},'flagvalid',weightsvalid,'grouppredict',finalgrouppredict,'flagpredict',weightspredict,...
    'membershipprob',{result1.prior},'misclassif',{options.misclassif},'groupmisclasprob',{misclasprobpergroup},...
    'avemisclasprob',{misclasprob},'class',{'RSIMCA'},'classic',{classicout},'compare',{compareout},'x',x,'group',group);

if size(x,2)>3
    result=rmfield(result,{'x','group'});
end

%Plots:
try
    if options.plots
        makeplot(result);
    end
catch %output must be given even if plots are interrupted 
    %> delete(gcf) to get rid of the menu 
end

%---------------------------------------------
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
        H0 = GRes{iClass}.Hsubsets.H0;
        H1 = GRes{iClass}.Hsubsets.H1;
        Hfreq = GRes{iClass}.Hsubsets.Hfreq;
        Hsets = [H0;H1;Hfreq];
        Hsetsmini = RemoveObsHsets(Hsets,i);
        same.value = 0;
        if isempty(find(H0 == i))
            if teller_if_lus >= 1
                same.value = 1;
            end
            teller_if_lus = teller_if_lus + 1;
        end
        interimRes = removeObsRobpca(groupi,i,GRes{iClass}.k,...
            Hsetsmini,same);
        if isempty(find(H0 == i))
            same.res = interimRes;
        end
        GRes{iClass}.M = interimRes.muk_min_i;
        GRes{iClass}.P = interimRes.Pk_min_i;
        GRes{iClass}.L = interimRes.Lk_min_i;
        GRes{iClass}.T = (groupia - ones(model.counts(iClass)-1,1)*GRes{iClass}.M)*GRes{iClass}.P;
        GRes{iClass}.h = GRes{iClass}.h - 1;
        outDist = CompDist(groupia,GRes{iClass});
        GRes{iClass}.cutoff.od = outDist.cutoff.od;
        GRes{iClass}.cutoff.sd = outDist.cutoff.sd;
       
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
function [odsc,sdsc] = testmodel(model,validx,validgroup)

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


%-----------------------------------------------------------------------------------------------------------
function Hsets_min_i = RemoveObsHsets(Hsets,i)

% removes the right index from the $h$-subsets in Hsets to 
% obtain (h - 1)-subsets.
% every h-set is put as a row in Hsets.
% i is the index of the observation that is removed from the whole data.

for r = 1:size(Hsets,1)
    if ~isempty(find(Hsets(r,:)== i))
        Hsets_min_i(r,:) = removal(Hsets(r,:),0,find(Hsets(r,:) == i));
    else
        Hsets_min_i(r,:) = Hsets(r,1:(end-1));
    end

    for j = 1:length(Hsets_min_i(r,:))
        if Hsets_min_i(r,j) > i
            Hsets_min_i(r,j) = Hsets_min_i(r,j) - 1;
        end
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

%----------------------------

function outDist = CompDist(data,out)

% Calculates the score and orthogonal distances 
% input: data : the original data
%        out is a structure that contains the results of the PCA.


[n,p] = size(data);
r = rank(data);
k = out.k;

% Computing distances 
% Robust score distances 
out.sd=sqrt(mahalanobis(out.T,zeros(size(out.T,2),1),'cov',out.L))';
out.cutoff.sd=sqrt(chi2inv(0.975,out.k));
% Orthogonal distances 
XRc=data-repmat(out.M,n,1);
Xtilde=out.T*out.P';
Rdiff=XRc-Xtilde;
out.od = [];
for i=1:n
    out.od(i,1)=norm(Rdiff(i,:));
end
% Robust cutoff-value for the orthogonal distance
if k~=r
    [m,s]=unimcd(out.od.^(2/3),out.h);
    out.cutoff.od = sqrt(norminv(0.975,m,s).^3); 
else
    out.cutoff.od=0;
end

outDist = out;

