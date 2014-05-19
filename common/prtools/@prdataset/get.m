%GET Get dataset parameter fields 
%
%   [VALUE1,VALUE2,...] = GET(A,FIELD1,FIELD2,...)
%
% INPUT 
%   A       Dataset
%   FIELDx  Field names (strings)
%
% OUTPUT
%   VALUEx  Field values
%
% DESCRIPTION
% Get parameter fields (given as strings in FIELDx) of the dataset A:
%   DATA        datavectors
%   LABELS      labels of the datavectors
%   NLAB        numerical labels, index in lablist
%   FEATDOM     feature domains
%   FEATLAB     feature labels
%   PRIOR       prior probabilities
%   COST        classification cost matrix
%   LABLIST     labels of the classes
%   TARGETS     dataset with soft labels or targets
%   LABTYPE     label type: 'crisp','soft' or 'target'
%   OBJSIZE     number of objects or vector with its shape
%   FEATSIZE    number of features or vector with its shape
%   IDENT       identifier for objects
%   VERSION     version field
%   NAME        string with name of the dataset
%   USER        user field
% These names may be supplied in upper- or lowercase.
%
% EXAMPLES
% [DATA,NLAB] = GET(A,'data','nlab')
%
% SEE ALSO
% DATASETS, DATASET
