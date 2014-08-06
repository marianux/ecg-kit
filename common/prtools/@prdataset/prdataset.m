%PRDATASET Dataset class constructor
%
%    A = PRDATASET(DATA,LABELS)
%
% INPUT
%   DATA    size [M,K]  a set of M datavectors (objects) of length K.
%                       a cell array of datasets will be concatenated. 
%   LABELS  size [M,N]  array with labels for the M datavectors.
%                       They should be either integers or character strings.
%                       Choose single characters for the fastest implementation.
%                       Numeric labels with value NaN or character labels
%                       with value CHAR(0) are interpreted as missing labels.
%
% OUTPUT
%    A      Dataset
%
% DESCRIPTION
% This command is the class constructor for datasets. In addition to the object labels  
% various other types of information can be stored in the fields of A.
% These fields are:
%
% DATA    size [M,K]  array (doubles) with M K-dimensional feature vectors (objects)
% FEATLAB size [K,F]  array with labels for the K features
% FEATDOM size [K]    cell array with domain description for the K features
% TARGETS size [M,C]  dataset with soft labels or targets
% PRIOR   size [C,1]  prior probabilities for each of the C classes
%                     - PRIOR = 0: all classes have equal probability 1/C
%                     - PRIOR = []: all datavectors are equally probable
% COST    size [C,C+1] Classification cost matrix. COST(I,J) are the costs
%                     of classifying an object from class I as class J.
%                     Column C+1 generates an alternative reject class and
%                     may be omitted, yielding a size of [C,C]. 
%                     An empty cost matrix, COST = [] (default) is interpreted
%                     as COST = ONES(C) - EYE(C) (identical costs of
%                     misclassification).
% LABLIST size [C,N]  class labels corresponding to the unique labels found
%                     in LABELS and thereby to the classes in the dataset.
%                     The order of the items in LABLIST corresponds to the
%                     apriori probablities stored in PRIOR. LABLIST should
%                     only be given explicitely if PRIOR is given and if it
%                     is not equal to 0 and not empty.
% LABTYPE             String defining the label type,
%                     'crisp' for defining classes by integers or strings
%                     'soft' for defining memberships to classes. In this
%                             case LABELS should be a MxC array with numbers
%                             between 0 and 1.
%                     'targets' for defining regression type target values.
%                             Labels should be a MxN numeric array for
%                             defining N targets per object.
% OBJSIZE             number of objects, or vector with its shape. This is
%                     useful if the set of objects can be interpreted as an
%                     image (objects are pixels).
% FEATSIZE            number of features, or vector with its shape. This is
%                     useful if the set of features can be interpreted as an
%                     image (features are pixels).
% IDENT  [M,1]        Cell array, identifier for objects. 
% NAME                String with dataset name
% USER                User definable variable
% VERSION             Date and PRTOOLS version at creation
%
% The fields LABLIST, OBJSIZE, FEATSIZE, IDENT and VERSION are preset by PRTOOLS. 
% The other fields can be set by the user by the below SET commands.
% All fields can be read by GET commands. By STRUCT(A) a dataset A can be
% converted to a structure. By DOUBLE(A) or +A the data can be retrieved.
% HELP DATASETS lists more information.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>) 
% DATASETS, MAPPINGS
