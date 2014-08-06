%CROSSVAL Error/performance estimation by cross validation (rotation)
% 
%   [ERR,CERR,NLAB_OUT] = CROSSVAL(A,CLASSF,NFOLDS,1,TESTFUN)
%   [ERR,STDS]          = CROSSVAL(A,CLASSF,NFOLDS,NREP,TESTFUN)
%   [ERR,CERR,NLAB_OUT] = CROSSVAL(A,CLASSF,NFOLDS,'DPS',TESTFUN)
%   R                   = CROSSVAL(A,[],NFOLDS,0)
%
% INPUT
%   A          Input dataset
%   CLASSF     The untrained classifier to be tested.
%   NFOLDS     Number of folds 
%              (default: number of samples: leave-one-out)
%   NREP       Number of repetitions (default: 1)
%   TESTFUN    Mapping,evaluation function (default classification error)
%
% OUTPUT
%   ERR        Average test error or performance weighted by class priors.
%   CERR       Unweighted test errors or performances per class
%   NLAB_OUT   Assigned numeric labels
%   STDS       Standard deviation over the repetitions.
%   R          Index array with rotation set
%
% DESCRIPTION
% Cross validation estimation of the error or performance (defined by TESTFUN)
% of the untrained classifier CLASSF using the dataset A. The set is randomly
% permutated and divided in NFOLDS (almost) equally sized parts, using a 
% stratified procedure. The classifier is trained on NFOLDS-1 parts and the 
% remaining part is used for testing. This is rotated over all parts. ERR
% is 
% their weighted avarage over the class priors. CERR are the class error 
% frequencies. The inputs A and/or CLASSF may be cell arrays of datasets and 
% classifiers. In that case ERR is an array with on position ERR(i,j) the 
% error or performance of classifier j for dataset i. In this mode CERR and 
% NLAB_OUT are returned in cell arrays.
%
% For NREP > 1 the mean error(s) over the repetitions is returned in ERR
% and the standard deviations in the observed errors in STDS.
%
% If NREP == 'DPS', crossvalidation is done by density preserving data
% splitting (DPS). In this case NFOLD should be a power of 2.
%
% In case NREP == 0 an index array is returned pointing to a fold for every
% object. No training or testing is done. This is useful for handling
% training and testing outside CROSSVAL.
%
% Note that this routine is identical to the PRCROSSVAL rouitne located in 
% the PRTools main directory. It thereby avoids confusion with the CROSSVAL 
% routine in the Stats toolbox if called by a dataset as a first parameter. 
% Users should preferably call PRCROSSVAL for the PRTools routine and use 
% CROSSVAL for the Stats version.
%
% REFERENCES
% 1. R. Kohavi: A Study of Cross-Validation and Bootstrap for Accuracy
% Estimation and Model Selection. IJCAI 1995: 1137-1145.
% 2. M. Budka, B. Gabrys, Correntropy-based density-preserving data
% sampling as an alternative to standard cross-validation, IJCNN2010, 1-8
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATASETS, MAPPINGS, DPS, CLEVAL, TESTC, PRCROSSVAL
