%DATAFILES Info on the datafile class construction for PRTools
%
% This is not a command, just an information file.
%
% Datafiles in PRTools are in the MATLAB language defined as objects of the
% class PRDATAFILE. They inherit most of their properties of the class PRDATASET.
% They are a generalisation of this class allowing for large datasets 
% distributed over a set of files. Before conversion to a dataset 
% preprocessing can be defined. There are four types of datafiles:
% raw   : Every file is interpreted as a single object in the dataset. These
%         files may, for instance, be images of different size.
% cell  : All files should be mat-files containing just a single variable being
%         a cell array. Its elements are interpreted as objects. The file names 
%         will be used as labels during construction. This may be changed by the
%         user afterwards.
% pre-cooked  : In this case the user should supply a command that reads a
%          file and converts it to a dataset.
% half-baked  : All files should be mat-files, containing a single dataset. 
% mature      : This is a datafile by PRTools, using the SAVEDATAFILE command after
%              execution of all preprocessing defined for the datafile. 
%
% A datafile is, like a dataset, a set consisting of M objects, each described
% by K features. K might be unknown, in which case it is set to zero, K=0.
% Datafiles store an administration about the files or directories in which
% the objects are stored. In addition they can store commands to preprocess
% the files before they are converted to a dataset and postprocessing
% commands, to be executed after conversion to a dataset.
%
% Datafiles are mainly an administration. Operations on datafiles are
% possible as long as they can be stored (e.g. filtering of images for raw
% datafiles, or object selection by GENDAT). Commands that are able to
% process objects sequentially, like NMC and TESTC can be executed on
% datafiles.
%
% Whenever a raw datafile is sufficiently defined by pre- and postprocessing
% it can be converted into a dataset. If this is still a large dataset, not
% suitable for the available memory, it should be stored by the
% SAVEDATAFILE command and is ready for later use. If the dataset is
% sufficiently small it can be directly converted into a dataset by
% PRDATASET.
%
% Intermediate results of datafiles that by the defined preprocessing
% cannot yet be converted into a dataset, can be stored as a new, raw
% datafile by CREATEDATAFILE.
%
% The main commands specific for datafiles are:
% DATAFILE       : constructor. It defines a datafile on a directory.
% ADDPREPROC     : adds preprocessing commands (low level command)
% ADDPOSTPROC    : adds postprocessing commands (low level command)
% FILTM          : user interface to add preprocessing to a datafile.
% CREATEDATAFILE : executes all defined preprocessing and stores the result
%                  as a new, raw datafile.
% SAVEDATAFILE   : executes all defined pre- and postprocessing and stores
%                  the result as a dataset in a set of matfiles.
% DATASET        : conversion to dataset
%
% Datafiles have the following fields, in addition to all dataset fields.
% ROOTPATH       : Absolute path of the datafile
% FILES          : names of directories (for raw datafiles) or mat-files
%                  (for converted datafiles)
% TYPE           : datafile type
% PREPROC        : preprocessing commands in a struct array
% POSTPROC       : postprocessing commands as mappings
% DATASET        : stores all dataset fields. Note that the DATA field as well
%                  as the target field are empty and that the IDENT.FILE_INDEX 
%                  field is used to store for every object a pointer to a file 
%                  or directory in FILES.
%
% Almost all operations defined for datasets are also defined for
% datafiles, with a few exceptions. Also fixed and trained mappings can
% handle datafiles, as they process objects sequentially. The use of
% untrained mappings in combination with datafiles is a problem, as they
% have to be adapted to the sequential use of the objects. Mappings that
% can handle datafiles are indicated in the Contents file.
%
% Subscription of datafiles is only defined for the first arguement, the
% objects, e.g. A(M,:) or even, irregulary, A(M) refer to object number M.
% As the objects in datafiles (e.g. images or time signals) may have different
% lengths, the second subscript, for datasets refering to the feature
% number, is undefined. A(M,N) causes an error of any N. Formally the feature 
% size of a dataset is set to 0. Checking of feature sizes in applying mappings
% to datafiles is disabled.
%
% The possibility to define preprocessing of objects (e.g. images) with
% different sizes makes datafiles useful for handling raw data and
% measurements of features.
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRDATAFILE, ADDPREPROC, ADDPOSTPROC, FILTM, FILTIM, CREATEDATAFILE,
% SAVEDATAFILE

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

