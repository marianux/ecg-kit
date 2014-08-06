%PRDATAFILE Datafile class constructor. This is an extension of PRDATASET.
%
%    A = PRDATAFILE(DIRECTORY,TYPE,READCMD,P1,P2,P3, ...)
%    A = PRDATAFILE(DIRECTORY,READCMD,P1,P2,P3, ...)
%
% INPUT
%   DIRECTORY - Data directory
%   TYPE      - Datafile type (default 'raw')
%   READCMD   - Command (m-file) for reading files in DIRECTORY.
%               Default: IMREAD
%   P1,P2,P3  - Optional parameters of READCMD
%
% OUTPUT
%   A         - Datafile
%
% DESCRIPTION
% Datafiles prepare and enable the handling of datasets distributed over
% multiple files, i.e. all files of DIRECTORY. Datafiles inherit all
% dataset fields. Consequently, most commands defined on datasets also
% operate on datafiles with the exception of a number of trainable
% mappings. There are five types of datafiles defined (TYPE):
% 'raw'          Every file is interpreted as a single object in the
%                dataset. All objects in the same sub-directory of
%                DIRECTORY receive the name of that sub-directory as class
%                label. Files may be preprocessed before conversion to
%                dataset by FILTM. At conversion time they should have the
%                same size (number of features).
% 'cell'         All files in DIRECTORY should be mat-files containing just
%                a single variable being a cell array. Its elements are
%                interpreted as objects. The file names will be used as
%                labels during construction. This may be changed by the
%                user afterwards.
% 'pre-cooked'   It is expected that READCMD outputs for all files a
%                dataset with the same label list and the same feature size.
% 'half-baked'   All files in DIRECTORY should be mat-files,
%                containing a single dataset. All datasets should have the 
%                same label list and the same feature size.
% 'mature'       This is a datafile directory constructed by SAVEDATAFILE. 
%                It executes all processing before creation.
%
%
% For all datafile types holds that execution of mappings (by FILTM or 
% otherwise and conversion to a dataset (by DATASET) is postponed as long as 
% possible. Commands are stored inside one of the datafile fields. 
% Consequently, errors might be detected at a later stage.
%
% The construction by DATAFILE still might be time consuming as for some types
% all files have to be checked. For that reason PRTools attempts to save a 
% mat-file with the DATAFILE definition in DIRECTORY. If it is encountered, it 
% is loaded avoiding a redefinition. 
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% DATAFILES, DATASETS, MAPPINGS, FILTM, SAVEDATAFILE, CREATEDATAFILE
