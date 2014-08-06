%MAPPINGS Info on the mapping class construction of PRTools
% 
% This is not a command, just an information file.
% 
% Mappings in PRTools are in the MATLAB language defined as objects of the
% class PRMAPPING. In the text below, the words 'object' and 'class' are used
% in the pattern recognition sense.
%
% In the Pattern Recognition Toolbox PRTools, there are many commands to
% define, train and use mappings between spaces of different (or equal)
% dimensionalities. Mappings operate mainly on datasets, i.e. variables of
% the type DATASET (see also DATASETS) and generate datasets and/or other
% mappings. For example:
%  
%  if    A  is a M x K dataset (M objects in  a K-dimensional space)
%  and   W  is a K x N mapping (a map from K to N dimensions)
%  then A*W is a M x N dataset (M objects in  a N-dimensional space)
%
% This is enabled by overloading the *-operator for the MAPPING variables.
% A*W is executed by PRMAP(A,W) and may also be called as such.
%
% Mappings can be linear (e.g. a rotation) as well as nonlinear (e.g. a
% neural network). Typically they are used to represent classifiers. In that
% case, a K x C mapping maps a K-feature data vector on the output space of
% a C-class classifier (an exception: some 2-class classifiers, like the
% discriminant functions may be implemented by a mapping onto a 1-dimensional 
% space determined by the distance to the discriminant).
% 
% Mappings are of the data-type MAPPING (CLASS(W) is a MAPPING), have a size
% of K x C if they map from K to C dimensions. Four types of mapping are 
% defined:
%
% untrained,  V = A*W.
%
%   Trains the untrained mapping W, resulting in the trained mapping V. W
%   has to be defined by W = PRMAPPING(MAPPING_FILE,{PAR1, PAR2}), in which
%   MAPPING_FILE is the name of the routine that executes the training and
%   PAR1, and PAR2 are two parameters that have to be included into the call
%   to MAPPING_FILE. Consequently, A*W is executed by PRTools as
%   MAPPING_FILE(A,PAR1,PAR2).
%
%   Example: train the 3-NN classifier on the generated data.
%
% 	  W = knnc([],3);         % untrained classifier
% 	  V = gendatd([50 50])*W; % trained classifier
% 	
% trained,  D = B*V	
%
%   Maps the dataset B on the trained mapping or classifier V, e.g. as
%   trained above. The resulting dataset D has as many objects (rows) as A,
%   but its feature size is now C if V is a K x C mapping. Typically, C is
%   the number of classes in the training set A or a reduced number of
%   features determined by the the training of V. V is defined by
%   V = PRMAPPING(MAPPING_FILE,'trained',DATA,LABELS,SIZE_IN,SIZE_OUT),
%   in which the MAPPING_FILE is the name of the routine that executes the
%   mapping, DATA is a field in which the parameters are stored (e.g.
%   weights) for the mapping execution, LABELS are the feature labels to be
%   assigned to the resulting dataset D = B*V (e.g. the class names) and
%   SIZE_IN and SIZE_OUT are the dimensionalities of the input and output
%   spaces. They are used for error checking only. D = B*V is executed by
%   PRTools as MAPPING_FILE(B,W). Example:
%
%     A = gendatd([50 50],10);	% generate random 10D datasets
%     B = gendatd([50 50],10);
%     W = klm([],0.9);          % untrained mapping, Karhunen-Loeve projection
%     V = A*W;                  % trained mapping V
%     D = B*V;                  % the result of the projection of B onto V
%
% fixed,  D = A*W 
%
%   Maps the dataset A by the fixed mapping W, resulting into a transformed
%   dataset D. Examples are scaling and normalization, e.g. W =
%   PRMAPPING('SIGM','fixed',S) defines a fixed mapping by the sigmoid function 
%   SIGM a scaling parameter S. A*W is executed by PRTools as SIGM(A,S). 
%
%   Example: normalize the distances of all objects in A such that their 
%   city block distances to the origin are one.
%
%     A = gendatb([50 50]);
%     W = normm;
% 	  D = A*W;
%
% combiner,  U = V*W 
%
%   Combines two mappings. The mapping W is able to combine itself with V
%   and produces a single mapping U. A combiner is defined by
%   W = PRMAPPING(MAPPING_FILE,'combiner',{PAR1,PAR2})
%   in which MAPPING_FILE is the name of the routine that executes the
%   combining and PAR1, and PAR2 are the parameters that have to be included
%   into the call to the MAPPING_FILE. Consequently, V*W is executed by
%   PRTools as MAPPING_FILE(V,PAR1,PAR2). In a call as D = A*V*W, first B =
%   A*V is resolved and may result in a dataset B. Consequently, W should be
%   able to handle datasets, and MAPPING_FILE is now called by
%   MAPPING_FILE(B,PAR1,PAR2) Remark: the combiner construction is not
%   necessary, since PRTools stores U = V*W as a SEQUENTIAL mapping (see
%   below) if W is not a combiner. The construction of combiners, however,
%   may increase the transparency for the user and efficiency in
%   computations. Example:
%
%     A = gendatd([50 50],10); % generate random 10D datasets
%     B = gendatd([50 50],10);
%     V = klm([],0.9);         % untrained Karhunen-Loeve (KL) projection
%     W = ldc;                 % untrained linear classifier LDC
%     U = V*W;                 % untrained combiner
%     T = A*U;                 % trained combiner
%     D = B*T;                 % apply the combiner (first KL projection, 
%                              %       then LDC) to B
%
% Differences between the four types of mappings are now summarized for
% a dataset A and a mapping W:
% A*W	  	-  untrained : results in a mapping
% 	  		-  trained   : results in a dataset, size checking
% 		  	-  fixed     : results in a dataset, no size checking
% 			  -  combiner  : treated as fixed      
%
% Suppose V is a fixed mapping, then for the various possibilities of
% the mapping W, the following holds:
% A*(V*W) -  untrained : evaluated as V*(A*V*W), resulting in a mapping
%         -  trained   : evaluated as A*V*W, resulting in a dataset
%         -  fixed     : evaluated as A*V*W, resulting in a dataset
%         -  combiner  : evaluated as A*(V*W), resulting in a dataset
%         
% Suppose V is an untrained mapping, then for the various possibilities of
% the mapping W holds:
% A*(V*W) -  untrained : evaluated as A*V*(A*(A*V)*W), resulting in a mapping
%         -  trained   : evaluated as A*V*W, resulting in a mapping
%         -  fixed     : evaluated as A*V*W, resulting in a mapping
%         -  combiner  : evaluated as A*(V*W), resulting in a mapping     
%
% Suppose V is a trained mapping, then for the various possibilities of
% the mapping W holds:
% A*(V*W) -  untrained : evaluated as V*(A*V*W), resulting in a mapping
%         -  trained   : evaluated as A*V*W, resulting in a dataset
%         -  fixed     : evaluated as A*V*W, resulting in a dataset
%         -  combiner  : evaluated as A*(V*W), resulting in a dataset
%
% The data fields stored in the MAPPING W = A*QDC can be found by
%
% 	STRUCT(W)
%
% which may display:
%
% 	MAPPING_FILE: 'normal_map'
% 	MAPPING_TYPE: 'trained'
%   DATA:          [1x1 struct]
%   LABELS:        [2x1 double]
%   SIZE_IN:       2
%   SIZE_OUT:      2
%   SCALE:         1
%   COST:          []
%   OUT_CONV:      0
%   NAME:          []
%   USER:          []
%   VERSION:       {1x2 cell  }
%
% These fields have the following meaning:
% MAPPING_FILE: Name of the m-file that executes the mapping.
% MAPPING_TYPE: Type of mapping: 'untrained','trained','fixed' or 'combiner'.
% DATA:         Parameters or data for handling or executing the mapping.
% LABELS:       Label list used as FEATLAB for labeling the features of the
%               output PRDATASET.
% SIZE_IN:      Expected input dimensionality of the data to be mapped.
%               If not set, it is neglected, otherwise it is used for the error
%               checking and display of the mapping size on the command line.
% SIZE_OUT:     Dimensionality of the output space. It should correspond to the
%               size of LABLIST. SIZE_OUT may be size vector, e.g. describing
%               the size of an image. See also the FEATSIZE field of PRDATASET.
% SCALE:        Output multiplication factor. If SCALE is a scalar all
%               multiplied by it. SCALE may also be a vector with size as
%               defined by SIZE_OUT to set separate scalings for each output.
% COST:         Classification costs in case the mapping defines a classifier.
% OUT_CONV:     Defines for trained and fixed mappings the output conversion:
%               - 0 : no conversion (to be used for mappings that output 
%               confidences or densities;
%               - 1 : sigmoid (for discriminants that output distances);
%               - 2 : normalisation (for converting densities and confidences
%               into posterior probability estimates;
%               - 3 : for performing sigmoid as well as normalisation.
% NAME:         Name of the mapping, used for informing the user on the 
%               command line, as well as for annotating plots.
% USER:         User field, not used by PRTools.
% VERSION:      Some information related to the version of PRTools used for
%               the mapping definition.
%
%
% The fields can be set by commands like SETMAPPING_FILE, SETDATA, SETLABELS, 
% SETSIZE, and may be retrieved by commands like GETMAPPING_FILE, GETDATA,
% GETLABELS and SETSIZE. Information stored in a mapping can be found
% as follows:
% - By DOUBLE(W) and by +W the content of the W.DATA is returned.
% - DISPLAY(W) writes the size of the mapping, the number of classes and the
% label type on the terminal screen.
% - SIZE(W) returns dimensionalities of input space and output space.
% - SCATTERD(A) makes a scatter-plot of a dataset.
% - SHOW(W) may be used to display images that are stored in mappings with 
% the MAPPING_FILE 'affine'. 
% - Using the dot extension as for structures, e.g. NAME = W.NAME;
% - The routines ISAFFINE, ISCLASSIFIER, ISCOMBINER, ISEMPTY, ISFIXED, 
% ISTRAINED and ISUNTRAINED test on some mapping types and states.
%
% Some standard MATLAB operations have been overloaded for variables of the 
% type PRMAPPING. They are defined as follows:
% W'      Defined for affine mappings only. It returns a transposed mapping.
% [W,V] 	Builds a combined classifier (see STACKED) operating in the same
%         feature space. A * [W V] = [A*W A*V].
% [W;V] 	Builds a combined classifier (see PARALLEL) operating in different
%         feature spaces: [A B] * [W;V] = [A*W B*V]. W and V should be 
%         mappings that correspond to the feature sizes of A and B.
% A*W     Maps a DATASET A by the MAPPING W. This is executed by PRMAP(A,W).
% V*W     Combines the mappings V and W sequentially. This is executed by
%         SEQUENTIAL(V,W).
% W+CON   Adding a constant is defined for affine mappings only. 
% W(:,K)	Output selection. If W is a trained mapping, just the features 
% 		    listed in K are returned.
%
%
% SEE ALSO (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% PRMAPPING, CLASSC, CNORMC, LABELD, 
%
% - Classifiers
%   NMC, KNNC, UDC, LDC, QDC, MOGC, QUADRC, FISHERC, PARZENC, PARZENDC,
%   TREEC, LOGLC, NAIVEBC, SVC, RBSVC, NUSVC, LIBSVC, TREEC, PERLC, BPXNC, 
%   RBNC, LMNC, WEAKC, STUMPC, SUBSC, ADABOOSTC, BAGGINGC, FDSC, VPC, DRBMC
%
% - Classifier Combiners
%   STACKED, PARALLEL, SEQUENTIAL, MEANC, AVERAGEC, PRODC, MEDIANC, MINC,
%   MAXC, VOTEC, WVOTEC, MODSELC, DCSC, RSSCC, MLRC, NAIVEBCC, TRAINCC
%
% - Density Estimation
%   GAUSSM, PARZENM, KNNM, 
%   
% - Dimension Reduction
%   FEATSEL, FEATSELB, FEATSELF, FEATSELI, FEATSELLR, FEATSELM, FEATSELO, 
%   FEATSELP, FEATSELV, BHATM, FISHERM, CHERNOFFM, KLM, KLMS, NLFISHERM,
%   PCAM, REDUCM
%
% - Scaling
%   SCALEM, CMAPM, SIGM, INVSIGM, NORMM
% 
% - Set commands
%   SETBATCH, SETCOST, SETDATA, SETLABELS, SETMAPPING_FILE,SETMAPPING_TYPE,
%   SETNAME, SETOUT_CONV, SETPOSTPROC, SETSCALE, SETSIZE, SETSIZE_IN,
%   SETSIZE_OUT, SETUSER
%
% - Get commands
%   GETBATCH, GETCOST, GETDATA, GETLABELS, GETMAPPING_FILE,, GETMAPPING_TYPE,
%   GETNAME, GETOUT_CONV, GETSCALE, GETSIZE, GETSIZE_IN, GETSIZE_OUT, GETUSER
% 
% - Tests
%   ISAFFINE, ISCLASSIFIER, ISEMPTY, ISPARALLEL, ISSTACKED, ISSEQUENTIAL,
%   ISTRAINED, ISUNTRAINED, ISFIXED, ISCOMBINER, 

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands

      
