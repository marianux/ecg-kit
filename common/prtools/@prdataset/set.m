%SET Set dataset fields 
%
%    A = SET(A,VARARGIN)
%
% Sets dataset fields (given as strings in VARARGIN) of the dataset A.
% E.G.: A = SET(A,'data',DATA,'nlab',NLAB).
% This is not different from using the field specific routines
% (e.g. SETDATA(A,DATA)
%
% List of parameter fields:
%
% DATA        datavectors
% LABELS      labels of the datavectors
% NLAB        nummeric labels, index in lablist
% FEATDOM     feature domains
% FEATLAB     feature labels
% PRIOR       prior probabilities
% COST        classification cost matrix
% LABLIST     labels of the classes
% TARGETS     dataset with soft labels or targets
% LABTYPE     label type: 'crisp','soft' or 'target'
% OBJSIZE     number of objects or vector with its shape
% FEATSIZE    number of features or vector with its shape
% IDENT       identifier for objects
% VERSION     version field
% NAME        string with name of the dataset
% USER        user field
%
% See datasets, dataset for more information
