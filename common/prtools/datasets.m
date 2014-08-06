%DATASETS Info on the dataset class construction for PRTools
%
% This is not a command, just an information file.
%
% Datasets in PRTools are in the MATLAB language defined as objects of the
% class PRDATASET. Below, the words 'object' and 'class' are used in the pattern 
% recognition sense.
%
% A dataset is a set consisting of M objects, each described by K features. 
% In PRTools, such a dataset is represented by a M x K matrix: M rows, each
% containing an object vector of K elements. Usually, a dataset is labeled.
% An example of a definition is:
%
%  DATA = [RAND(3,2) ; RAND(3,2)+0.5];
%  LABS = ['A';'A';'A';'B';'B';'B'];
%  A = PRDATASET(DATA,LABS)
%
% which defines a [6 x 2] dataset with 2 classes.
%
% The [6 x 2] data matrix (6 objects given by 2 features) is accompanied by
% labels, assigning each of the objects to one of the two classes A and B.
% Class labels can be numbers or strings and should always be given as rows
% in the label list. A lable may also have the value NaN or may be an empty
% string, indicating an ulabeled object. If the label list is not given, 
% all objects are marked as unlabeled.
%
% Various other types of information can be stored in a dataset. The most
% simple way to get an overview is by typing:
%
%  STRUCT(A)
%
% which for the above example displays the following:
%
%         DATA: [6x2 double]
%      LABLIST: [2x1 double]
%         NLAB: [6x1 double]
%      LABTYPE: 'crisp'
%      TARGETS: []
%      FEATLAB: [2x1 double]
%      FEATDOM: {1x2 cell}
%        PRIOR: []
%         COST: []
%      OBJSIZE: 6
%     FEATSIZE: 2
%        IDENT: {6x1 cell}
%      VERSION: {1x2 cell}
%         NAME: []
%         USER: []
%
% These fields have the following meaning:
% 
% DATA     : an array containing the objects (the rows) represented by  
%            features (the columns). In the software and help-files, the number
%            of objects is usually denoted by M and the number of features is
%            denoted by K. So, DATA has the size of [M,K]. This is also defined 
%            as the size of the entire dataset.
% LABLIST  : The names of the classes, stored row-wise. These class names
%            should be integers, strings or cells of strings. Mixtures of
%            these are not supported. LABLIST has as many rows as there are 
%            classes. This number is usually denoted by C. LABLIST is
%            constructed from the set of LABELS given in the DATASET command
%            by determining the unique names while ordering them alphabetically.
% NLAB     : an [M x 1] vector of integers between 1 and C, defining for each
%            of the M objects its class. They are indexing LABLIST.
% LABTYPE  : 'CRISP', 'SOFT' or 'TARGETS' are the three possible label types.
%            In case of 'CRISP' labels, a unique class, defined by NLAB, is
%            assigned to each object, pointing to the class names given in
%            LABLIST.
%            For 'SOFT' labels, each object has a corresponding vector of C 
%            numbers between 0 and 1 indicating its membership (or confidence 
%            or posterior probability) of each of the C classes. These numbers
%            are stored in the array TARGETS of the size M x C. They don't
%            necessarily sum to one for individual row vectors.
%            Labels of type 'TARGETS' are in fact no labels, but merely target
%            vectors of length C. The values are again stored in TARGETS and
%            are not restricted in value.
% TARGETS  : [M,C] array storing the values of the soft labels or targets.
% FEATLAB  : A label list (like LABLIST) of K rows storing the names of the
%            features.
% FEATDOM  : A cell array describing for each feature its domain.
% PRIOR    : Vector of length C storing the class prior probabilities. They 
%            should sum to one. If PRIOR is empty ([]) it is assumed that the
%            class prior probabilities correspond to the class frequencies.
% COST     : Classification cost matrix. COST(I,J) are the costs
%            of classifying an object from class I as class J. Column C+1
%            generates an alternative reject class and may be omitted, 
%            yielding a size of [C,C]. An empty cost matrix, COST = [] 
%            (default) is interpreted as COST = ONES(C) - EYE(C) (identical
%            costs of misclassification).
% OBJSIZE  : The number of objects, M. In case the objects are related to a
%            n-dimensional structure, OBJSIZE is a vector of length n, storing
%            the size of this structure. For instance, if the objects are pixels
%            in a [20 x 16] image, then OBJSIZE = [20,16] and M = 320.
% FEATSIZE : The number of features, K. In case the features are related to 
%            an n-dimensional structure, FEATSIZE is a vector of length n, 
%            storing the size of this structure. For instance, if the features
%            are pixels in a [20 x 16] image, then FEATSIZE = [20,16] and 
%            K = 320.
% IDENT    : A cell array of M elements storing indicators of the M objects.
%            They are initialized by integers 1:M.
% VERSION  : Some information related to the version of PRTools used for
%            defining the dataset.
% NAME     : A character string naming the dataset, possibly used to annotate
%            related graphics.
% USER     : Free field for the user, not used by PRTools.
%
%
% The fields can be set by commands like SETDATA, SETFEATLAB, SETLABELS,
% see below for a complete list.
% Note that there is no field LABELS in the DATASET definition. Labels are
% converted to NLAB and LABLIST. The command SETLABELS however exists and
% takes care of the conversion.
%
% The data and information stored in a dataset can be retrieved as follows:
%
% - By DOUBLE(A) and by +A, the content of A.DATA is returned.
% - [N,LABLIST] = CLASSSIZES(A); 
%   It returns the numbers of objects per class and the class names stored 
%   in LABLIST.
% - By DISPLAY(A), it writes the size of the dataset, the number of classes 
%   and the label type on the terminal screen.
% - By SIZE(A), it returns the size of A.DATA: numbers of objects and features.
% - By SCATTERD(A), it makes a scatter plot of a dataset.
% - By SHOW(A), it may be used to display images that are stored as features 
%   or as objects in a dataset. 
% - By commands like: GETDATA, GETFEATLAB, etcetera, see below. With some
%   exceptions they point to a single dataset field. E.g. GETSIZE(A) returns 
%   [M,K,C]. A aet of commands does not return data, but instead they return 
%   indices to objects that have specific identifiers, labels or class indices:
%   FINDIDENT, FINDLABELS, FINDNLAB.
%
% Many standard MATLAB operations and a number of general MATLAB commands have 
% been overloaded for variables of the DATASET type.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATASET, DATA2IM, OBJ2FEAT, FEAT2OBJ, IM2FEAT, IM2OBJ, DATAIM  
% SETDATA, SETFEATLAB, SETFEATDOM, SETFEATSIZE, SETIDENT, SETLABELS, 
% SETLABLIST, SETLABTYPE, SETNAME, SETNLAB, SETOBJSIZE, SETPRIOR, SETCOST, 
% SETTARGETS, SETUSER, SETLABLISTNAMES, SETVERSION
% GETDATA, GETFEATLAB, GETFEATDOM, GETFEATSIZE, GETIDENT, GETLABELS, 
% GETLABLIST, GETLABTYPE, GETNAME, GETNLAB, GETOBJSIZE, GETPRIOR, GETCOST, 
% GETSIZE, GETTARGETS, GETUSER, GETVERSION, GETCLASSI, GETLABLISTNAMES,
% FINDIDENT, FINDLABELS, FINDNLAB

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

