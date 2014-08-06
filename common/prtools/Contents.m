% Pattern Recognition Tools (<a href="http://37steps.com/prtools">PRTools Guide</a>)
% Version 5.1.1 14-May-2014
%
%Datasets and Mappings (just most important routines)
%---------------------
%prdataset      - Define dataset from datamatrix and labels
%datasets       - List information on datasets (just help, no command)
%prdatafile     - Define dataset from directory of object files 
%datafiles      - List information on datafiles (just help, no command)
%cat2data       - Create categorical dataset
%classnames     - Retrieve names of classes
%classsizes     - Retrieve sizes of classes
%feat2lab       - Label dataset by one of its features and remove this feature
%gencirc        - Generation of a one-class circular dataset
%genclass       - Generate class frequency distribution
%genlab         - Generate dataset labels
%getlab         - Retrieve object labels from datasets and mappings
%getnlab        - Retrieve nummeric object labels from dataset
%setfeatlab     - Set feature labels in dataset
%getfeatlab     - Get feature labels in dataset
%getfeat        - Retrieve feature labels from datasets and mappings
%setdat         - Change data in dataset for classifier output
%setdata        - Change data in dataset or mapping
%getdata        - Retrieve data from dataset or mapping
%setlabels      - Change labels of dataset or mapping
%getlabels      - Retrieve labels from a dataset
%setprior       - Reset class prior probabilities of dataset
%getprior       - Retrieve class prior probabilities from dataset
%addlabels      - Add additional labelling
%changelablist  - Change current active labeling
%misval         - Fix missing values in a dataset
%multi_labeling - List information on multi-labeling (help only)
%prmapping      - Define and retrieve mapping and classifier from data
%mappings       - List information on mappings (just help, no command)
%renumlab       - Convert labels to numbers
%matchlab       - Match different labelings
%prarff         - Convert ARFF file (WEKA) to PRTools dataset
%remclass       - Remove a class from a dataset
%seldat         - Retrieve a part of a dataset
%selclass       - Retrieve a class from a dataset 
%
%Data Generation (more in prdatasets)
%---------------
%circles3d   - Create a dataset containing 2 circles in 3 dimensions
%lines5d     - Create a dataset containing 3 lines in 5 dimensions
%gendat      - Random sampling of datasets for training and testing
%gensubsets  - Generation of a consistent series of subsets of a dataset
%gendatgauss - Generation of multivariate Gaussian distributed data
%gendatb     - Generation of banana shaped classes
%gendatc     - Generation of circular classes
%gendatd     - Generation of two difficult classes
%gendath     - Generation of Highleyman classes
%gendati     - Generation of random windows from images
%gendatk     - Nearest neighbour data generation
%gendatl     - Generation of Lithuanian classes
%gendatm     - Generation of 8 2d classes
%gendatp     - Parzen density data generation
%gendatr     - Generate regression dataset from data and target values
%gendats     - Generation of two Gaussian distributed classes
%gendatw     - Sample dataset by given weigths
%gendatv     - Generation of a very large dataset
%gentrunk    - Generation of Trunk's example
%prdata      - Read data from file
%seldat      - Select classes / features / objects from dataset
%spirals     - Generation of a two-class spiral dataset
%getwindows  - Get pixel feature vectors around given pixels in image dataset
%prdataset   - Read existing dataset from file
%prdatasets  - Overview and download of standard datasets
%
%Datafiles
%---------
%prdatafile     - Define datafile from set of files in directory
%createdatafile - Save datafile, store intermediate result as raw datafile
%savedatafile   - Save datafile, store intermediate result as mature datafile
%filtm          - Mapping for arbitrary processing of a datafile
%prdatafiles    - Overview and download of standard datafiles
%
%Linear and Quadratic Classifiers (*operate on datasets and datafiles)
%--------------------------------
%fisherc     - Minimum least square linear classifier
%ldc         - Normal densities based linear (muli-class) classifier
%loglc       - Logistic linear classifier
%nmc         - Nearest mean linear classifier
%nmsc        - Scaled nearest mean linear classifier
%quadrc      - Quadratic classifier
%qdc         - Normal densities based quadratic (multi-class) classifier
%udc         - Uncorrelated normal densities based quadratic classifier
%klldc       - Linear classifier based on KL expansion of common cov matrix
%pcldc       - Linear classifier based on PCA expansion on the joint data
%polyc       - Add polynomial features and run arbitrary classifier
%subsc       - Subspace classifier
%statslinc   - Linear classifier from the Stats toolbox
% 
%classc      - Converts a mapping into a classifier
%labeld      - Find labels of objects by classification
%logdens     - Convert density estimates to log-densities for more accuracy
%rejectc     - Creates reject version of exisiting classifier
%testc       - General error estimation routine for trained classifiers
%
%Other Classifiers 
%-----------------
%knnc        - k-nearest neighbour classifier (find k, build classifier)
%testk       - Error estimation for k-nearest neighbour rule
%edicon      - Edit and condense training sets
%statsknnc   - k-nearest neighbour classifier from the Stats toolbox
%
%weakc       - Weak classifier
%stumpc      - Decision stump classifier
%adaboostc   - ADABoost classifier
%
%parzenc     - Parzen classifier
%parzendc    - Parzen density based classifier
%testp       - Error estimation for Parzen classifier
%
%treec       - Construct binary decision tree classifier
%dtc         - Decision tree classifier, rewritten, also for nominal features
%statsdtc    - Decision tree classifier from the Stats toolbox
%randomforestc - Breiman's random forest classifier
%naivebc     - Naive Bayes classifier
%statsnbc    - Naive Bayes classifier from the Stats toolbox
%bpxnc       - Feed forward neural network classifier by backpropagation
%lmnc        - Feed forward neural network by Levenberg-Marquardt rule
%neurc       - Automatic neural network classifier
%perlc       - Linear perceptron 
%rbnc        - Radial basis neural network classifier
%rnnc        - Random neural network classifier
%ffnc        - Feed-forward neural net classifier back-end routine
%bagc        - Feature set classifier, e.g. for multiple-instance learning
%
%fdsc        - Feature based dissimilarity space classifier
%mdsc        - Manhatten distance feature based dissimilarity space classifier
%vpc         - Voted perceptron classifier
%drbmc       - Discriminative restricted Boltzmann machine classifier
%
%libsvc      - Support vector classifier by LIBSVM
%nulibsvc    - Support vector classifier by LIBSVM
%svc         - Support vector classifier
%svo         - Support vector optimizer
%nusvc       - Support vector classifier
%nusvo       - Support vector optimizer
%rbsvc       - Radial basis SV classifier
%kernelc     - General kernel/dissimilarity based classification
%
%Normal Density Based Classification
%-----------------------------------
%distmaha    - Mahalanobis distance
%meancov     - Estimation of means and covariance matrices from multiclass data
%nbayesc     - Bayes classifier for given normal densities
%ldc         - Normal densities based linear (muli-class) classifier
%qdc         - Normal densities based quadratic (multi-class) classifier
%udc         - Uncorrelated normal densities based quadratic classifier
%mogc        - Mixture of gaussians classification
%testn       - Error estimate of discriminant on normal distributions
%
%Feature Selection
%-----------------
%feateval    - Evaluation of a feature set
%featrank    - Ranking of individual feature permormances
%featsel     - Feature Selection
%featselb    - Backward feature selection
%featself    - Forward feature selection
%featsellr   - Plus-l-takeaway-r feature selection
%featseli    - Feature selection on individual performance
%featselm    - Feature selection map, general routine for feature selection
%featselo    - Branch and bound feature selection
%featselp    - Floating forward feature selection
%featselv    - Selection of varying features
%
%Classifiers and tests (general)
%-------------------------------
%bayesc      - Bayes classifier by combining density estimates
%classim     - Classify image using a given classifier
%classc      - Convert mapping to classifier
%labeld      - Find labels of objects by classification
%cleval      - Classifier evaluation (learning curve)
%clevalb     - Classifier evaluation (learning curve), bootstrap version
%clevalf     - Classifier evaluation (feature size curve)
%clevals     - Classifier evaluation (feature /learning curve), bootstrap
%confmat     - Computation of confusion matrix
%costm       - Cost mapping, classification using costs
%prcrossval  - Crossvalidation 
%cnormc      - Normalisation of classifiers
%disperror   - Display error matrix with information on classifiers and datasets
%labelim     - Construct image of labeled pixels
%logdens     - Convert density estimates to log-densities for more accuracy
%loso        - Leave_one_set_out crossvalidation
%mclassc     - Computation of multi-class classifier from 2-class discriminants
%regoptc     - Optimisation of regularisation and complexity parameters
%reject      - Compute error-reject trade-off curve
%prroc       - Receiver-operator curve (ROC)
%shiftop     - Shift operating point of classifier
%testc       - General error estimation routine for trained classifiers
%testd       - Error of dataset applied to given classifier
%testauc     - Estimate error as area under the ROC
%
%Mappings
%--------
%affine      - Construct affine (linear) mapping from parameters
%bhatm       - Two-class Bhattacharryya mapping
%cmapm       - Compute some special maps
%datasetm    - Mapping conversion dataset
%disnorm     - Normalization of a dissimilarity matrix
%featselm    - Feature selection map, general routine for feature selection
%fisherm     - Fisher mapping
%chernoffm   - Chernoff mapping
%invsigm     - Inverse sigmoid map
%filtm       - Arbitrary operation on datafiles/datasets, object by object
%mapm        - Arbitrary mapping operation on doubles and datasets
%gaussm      - Mixture of Gaussians density estimation
%kernelm     - Kernel mapping
%klm         - Decorrelation and Karhunen Loeve mapping (PCA)
%klms        - Scaled version of klm, useful for prewhitening
%knnm        - k-Nearest neighbor density estimation
%mclassm     - Computation of mapping from multi-class dataset
%prmap       - General routine for computing and executing mappings
%mappingtools - Macro defining some mappings
%nlfisherm   - Nonlinear Fisher mapping
%normm       - Object normalization map
%parzenm     - Parzen density estimation
%parzenml    - Optimization of smoothing parameter in Parzen density estimation.
%pcam        - Principal Component Analysis
%pcaklm      - Backend routine for PC and KL mappings
%proxm       - Proximity mapping and kernel construction
%reducm      - Reduce to minimal space mapping
%remoutl     - Remove outliers
%rejectm     - Creates rejecting mapping
%scalem      - Compute scaling data
%sigm        - Simoid mapping
%spatm       - Augment image dataset with spatial label information
%tsnem       - tSNE mapping
%sammonm     - Multi-dimensional scaling by Sammon mapping
%userkernel  - User supplied kernel definition
%
%gtm         - Fit a Generative Topographic Mapping (GTM) by EM
%plotgtm     - Plot a Generative Topographic Mapping in 2D
%som         - Simple routine computing a Self-Organizing Map (SOM)
%prplotsom   - Plot a Self-Organizing Map in 2D
%
%Classifier combiners
%--------------------
%averagec    - Combining linear classifiers by averaging coefficients
%baggingc    - Bootstrapping and aggregation of classifiers
%dcsc        - Dynamic Classifier Selecting Combiner
%modselc     - Model Selection Combiner (Static selection)
%rsscc       - Random subspace combining classifier
%votec       - Voting classifier combiner
%wvotec      - Weighted voting classifier combiner
%maxc        - Maximum classifier combiner
%minc        - Minimum classifier combiner
%meanc       - Mean classifier combiner
%medianc     - Median classifier combiner
%mlrc        - Muli-response linear regression combiner
%naivebcc    - Naive Bayes classifier combiner
%perc        - Percentile combiner
%prodc       - Product classifier combiner
%traincc     - Train combining classifier
%fixedcc     - Fixed combiner construction, back end
%parsc       - Parse classifier or map
%rejectc     - Creates reject version of exisiting classifier
%parallel    - Parallel combining of classifiers
%bagcc       - Feature set combining classifier
%stacked     - Stacked combining of classifiers
%sequential  - Sequential combining of classifiers
%
%
%Regression
%----------
%linearr     - Linear regression
%ridger      - Ridge regression
%lassor      - LASSO
%svmr        - Support vector regression
%ksmoothr    - Kernel smoother
%knnr        - k-nearest neighbor regression
%pinvr       - Pseudo-inverse regression
%plsr        - Partial least squares regression
%plsm        - Partial least squares mapping
%gpr         - Gaussian Process regression
%
%testr       - Mean squared regression error
%rsquared    - R^2-statistic
%
%Handling images in datasets and datafiles
%-----------------------------------------
%data2im     - Convert dataset to image
%getobjsize  - Retrieve image size of feature images in datasets
%getfeatsize - Retrieve image size of object images in datasets
%obj2feat    - Transform object images to feature images in dataset
%feat2obj    - Transform feature images to object images in dataset
%im2feat     - Convert image to feature in dataset
%im2obj      - Convert image to object in dataset
%imsize      - Retrieve size of specific image in datafile
%im_patch    - Find / generate patches in object images
%band2obj    - Convert image bands to objects in dataset
%bandsel     - Select image bands in dataset or datafile
%selectim    - Select image in multi-band object image dataset/datafile
%show        - Display objects in datasets, datafiles and mappings
%im_dbr      - Image Database Retrieval GUI
%
%Operations on images in datasets and datafiles
%----------------------------------------------
%classim     - Classify image using a given classifier
%doublem     - Convert datafile images into double
%filtim      - Image operation on objects in datafiles/datasets
%spatm       - Augment image dataset with spatial label information
%im_box            - Bounding box
%im_center         - Center image
%im_fft            - FFT transform (and more)
%im_gauss          - Gaussian filtering by Matlab
%im_gray           - Multi-band to gray-value conversion
%im_hist_equalize  - Histogram equalization
%im_invert         - Invert image
%im_label          - Labeling binary images
%im_norm           - Normalize images w.r.t. mean and variance
%im_resize         - Resize images
%im_rotate         - Rotate images
%im_scale          - Scale images
%im_select_blob    - Select largest blob
%im_stretch        - Contrast stretching of images
%im_threshold      - Threshold images
%im_unif           - Uniform filtering
%
%Feature extraction from images in datasets and datafiles
%--------------------------------------------------------
%histm         - Convert images to histograms. Trains the bin positions
%im_hist       - Convert images to histograms for fixed bin positions
%im_harris     - Find Harris points in images
%im_moments    - Computes moments as features from object images
%im_mean       - Computes center of gravity
%im_measure    - Computes some measurements
%im_profile    - Computes image profiles
%im_skel_meas  - Skeleton measurements
%im_stat       - Compute some simple statistics
%
%Clustering and distances
%------------------------
%distm       - Distance matrix between two data sets
%emclust     - Expectation - maximization clustering
%proxm       - Proximity mapping and kernel construction
%hclust      - Hierarchical clustering
%kcentres    - k-centres clustering
%prkmeans    - k-means clustering
%modeseek    - Clustering by modeseeking
%
%mds         - Non-linear mapping by multi-dimensional scaling (Sammon)
%mds_cs      - Linear mapping by classical scaling
%mds_init    - Initialisation of multi-dimensional scaling
%mds_stress  - Dissimilarity of distance matrices
%
%Plotting
%--------
%gridsize    - Set gridsize used in the PRTools plot commands
%plotc       - Plot discriminant function for two features
%plote       - Plot error curves
%plotf       - Plot feature distribution
%plotm       - Plot mapping
%ploto       - Plot object functions
%plotr       - Plot regression functions
%plotdg      - Plot dendrgram (see hclust)
%scatterd    - Scatterplot
%scatterdui  - Scatterplot scatterplot with feature selection
%scattern    - Simple, unannotated scatterplot, no axes.
%scatterr    - Scatter regression dataset
%
%Various tests and support routines
%----------------------------------
%cdats              - Support routine for checking datasets
%concatm            - Concatenate cell array of mappings or datasets ({} --> [])
%iscomdset          - Test on compatible datasets
%isdataim           - Test on image dataset
%isdataset          - Test on dataset
%isfeatim           - Test on feature image dataset
%ismapping          - Test on mapping
%isobjim            - Test on object image dataset
%issequential       - Test on sequential mapping
%isstacked          - Test on stacked mapping
%isparallel         - Test on parallel mapping
%issym              - Test on symmetric matrix
%isvaldset          - Test on valid dataset
%isvaldfile         - Test on valid datafile
%matchlablist       - Match entries of label lists
%mapex              - Train and execute mapping on the same dataset
%labcmp             - Compare two label lists and find the differences
%nlabcmp            - Compare two label lists and count the differences
%testdatasize       - Check datasize and convert datafile to dataset
%define_mapping     - Define empty mapping
%mapping_task       - Check mapping task
%trained_mapping    - Defined trained mapping
%trained_classifier - Define trained classifier
%setdefaults        - Substitute defaults
%shiftargin         - Conditional shift of input arguments
%prload             - Load prtools4 mat-files and convert to prtools5
%prtools4to5        - Convert prtools4 directory to prtools5
%
%Examples
%--------
%prex_cleval     - learning curves
%prex_combining  - classifier combining
%prex_confmat    - confusion matrix, scatterplot and gridsize
%prex_datafile   - datafile usage
%prex_datasets   - standard datasets
%prex_density    - Various density plots
%prex_eigenfaces - Use of images and eigenfaces
%prex_matchlab   - K-means clustering and matching labels
%prex_mcplot     - Multi-class classifier plot
%prex_plotc      - Dataset scatter and classifier plot
%prex_mds        - Multi-dimensional scaling and visualisation
%prex_som        - Training a SelfOrganizing Maps
%prex_spatm      - Spatial smoothing of image classification
%prex_cost       - Cost matrices and rejection
%prex_logdens    - Density based classifier improvement
%prex_soft       - Soft label example
%prex_regr       - Regression example
%
%prdownload  - low level routine for retrieving datasets
%prglobal    - set / list all globals and settings
%prversion   - returns version information on PRTools
%prwaitbar   - report PRTools progress by single waitbar
%prwarning   - control PRTools warning level
%prmemory    - controol PRTools large dataset handling
%prtver      - prtools version back end
%typp        - list prtools routine nicely
%
%--- <a href="http://37steps.com/prtools">PRTools Guide</a> ---

% Copyright: R.P.W. Duin, r.p.w.duin@37steps.com
% Faculty EWI, Delft University of Technology
% P.O. Box 5031, 2600 GA Delft, The Netherlands
