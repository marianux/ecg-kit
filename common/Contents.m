% COMMON
%
% Files
%   ADC2realunits                        - Convert adimentional sample values to real units
%   ADC2units                            - Convert adimentional sample values to target voltage units
%   addpath_if_not_added                 - (Internal) Add the path only if not was already added.
%   adjust_string                        - (Internal) Works with strings to center, trim and justify to a certain string width
%   Annotation_process                   - Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
%   AnnotationFilterConvert              - Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
%   arrow                                - (Internal) Creates arrows in grapics
%   autovec_calculation                  - (Internal) Calculate eigenvalues and vectors
%   autovec_calculation_robust           - (Internal) Calculate eigenvalues and vectors using robust covariance estimation
%   bandpass_filter_design               - MATLAB Code
%   BaselineWanderRemovalMedian          - Remove baseline wandering with the median estimation method
%   BaselineWanderRemovalSplines         - Remove baseline wandering with the median estimation method
%   bxb                                  - (Internal) Compares two heartbeat series and produce the confusion matrix as result
%   calc_btime                           - (Internal) Creates a string with the base time
%   calc_co_ocurrences                   - (Internal) Calculate the heartbeats co-ocurrences
%   calc_correlation_gain                - (Internal) Add the path only if not was already added.
%   CalcRRserieQuality                   - (Obsolete) Estimate the quality of QRS complex detections 
%   CalcRRserieRatio                     - Estimate the quality of QRS complex detections 

%   calculateSeriesQuality               - Estimate the quality of QRS complex detections 
%   cluster_data_with_EM_clust           - (Internal) Cluster data with expectation-maximization algorithm
%   colvec                               - (Internal) Reshape the input into a column vector
%   combine_anns                         - (Internal) Create new QRS detections based on other lead/algorithms detections
%   ConcatenateQRSdetectionPayloads      - (Internal) Concatenate two payloads.
%   cprintf                              - displays styled formatted text in the Command Window

%   DelayedCovMat                        - (Internal) Delayed covariance matrix calculation
%   deNaN_dataset                        - (Internal) Replace NaN from PRdatasets
%   design_downsample_filter             - (Internal) Design a filter to downsample signals prior printing to a report

%   disp_option_enumeration              - (Internal) Display an enumeration of options to a string
%   disp_string_framed                   - (Internal) Display a message framed to a string
%   disp_string_title                    - (Internal) Display a message framed to a string
%   DisplayConfusionMatrix               - (Internal) Pretty display in Screen the confusion matrix.
%   DisplayResults                       - (Internal) Pretty-Display results of a classification experiment
%   ds2nfu                               - Convert data space units into normalized figure units. 

%   ECGformat                            - Gets the format of an ECG recording filename
%   ECGtask                              - Defines the class interface for the ECGtask derived classes
%   ECGtask_do_nothing                   - Null ECGtask (for Matlab)
%   ECGtask_ECG_delineation              - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_ECG_delineation_corrector    - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_heartbeat_classifier         - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PCA_proj_basis               - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PPG_ABP_corrector            - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PPG_ABP_detector             - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_corrector                - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_detection                - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_detections_post_process  - ECGtask for ECGwrapper (for Matlab)
%   ECGwrapper                           - Allow acces to ECG recordings of arbitrary format and length.
%   example_worst_ever_ECG_delineator    - (Internal) Example of user-created QRS detector
%   example_worst_ever_QRS_detector      - (Internal) Example of user-created QRS detector
%   exist_distributed_file               - (Internal) Check the existence of a file in a distributed (slow) filesystem
%   get_ECG_idx_from_header              - (Internal) Guess ECG signals indexes in a multimodal recording
%   get_PPG_ABP_idx_from_header          - (Internal) Guess PPG/ABP signals indexes in a multimodal recording
%   getAnnNames                          - (Internal) Get names of annotations from annotation structure
%   GetFunctionInvocation                - (Internal) Create a string with the invocation of a function
%   GTHTMLtable                          - GTHTMLtable - Generate an HTML page with a table of a matrix.
%   init_ghostscript                     - (Internal) Init environment variables for using ghostscript
%   init_WFDB_library                    - (Internal) Init environment variables for using WFDB toolbox
%   isAHAformat                          - (Internal) Check if a recording is in ISHNE format.
%   isHESformat                          - (Internal) Check if a recording is in HES format.
%   isISHNEformat                        - (Internal) Check if a recording is in ISHNE format.
%   isMatlab                             - (Internal) Check if the kit is running on Matlab
%   isOctave                             - (Internal) Check if the kit is running on Octave

%   logit_function                       - (Internal) Logit function
%   matrix2positions                     - (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
%   max_index                            - (Internal) Index of the maximum element in a vector
%   maximize                             - Size a window to fill the entire screen.
%   MedianFilt                           - (Internal) Mean/Median filtering
%   modmax                               - (Internal) Find modulus maxima in a signal
%   my_colormap                          - (Internal) Create a colormap for signal visualization
%   my_ppval                             - PPVAL  Evaluate piecewise polynomial.
%   myzerocros                           - (Internal) Detect zero-crosses in a signal
%   nanmeda                              - (Internal) Calculate the median of absolute deviations from the median (MEDA)
%   pack_signal                          - (Internal) Example of user-created QRS detector
%   parse_pids                           - (Internal) Identify how many PIDs are in total and which is this PID, based on a string formatted this_pid/cant_pids

%   PeakDetection2                       - peaks = PeakDetection2(x,fs,wlen,fp1,fp2,th,flag),
%   PiCA                                 - [y,W,A] = PiCA(x,peaks1,peaks2)
%   plot_auc                             - (Internal) Plot the area under the ROC curve
%   plot_ecg_heartbeat                   - (Internal) obsolete, use plot_ecg_strip
%   plot_ecg_mosaic                      - Plots multidimentional signal in mosaic style
%   plot_ecg_strip                       - Plots and interact with ECG signal

%   plot_roc                             - (Internal) Plot the ROC curve
%   PlotGlobalWaveMarks                  - (Internal) Internal function of plot_ecg_strip
%   PlotWaveMarks                        - (Internal) Internal function of plot_ecg_strip

%   positions2matrix                     - (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
%   PrctileFilt                          - (Internal) Arbitrary percentile filtering
%   progress_bar                         - (Internal) A progress bar class for showing evolution of a process to users
%   progress_bar_ex                      - start of algorithm
%   QRScorrector                         - (Internal) GUI for correcting QRS detections
%   qs_filter_design                     - (Internal) Design the wavelet decomposition filters for wavedet algorithm
%   qs_wt                                - (Internal) Calculates the wavelet transform


%   rand_linespec                        - (Internal) Example of user-created QRS detector
%   read_310_format                      - (Internal) Read the MIT 310 format
%   read_311_format                      - (Internal) Read the MIT 311 format
%   read_AHA_ann                         - Reads ECG annotations in AHA format
%   read_AHA_format                      - Reads ECG recording in AHA format
%   read_AHA_header                      - Reads ECG header in AHA format
%   read_ECG                             - Reads an ECG recording
%   read_HES_ann                         - Reads ECG annotations in HES format
%   read_HES_format                      - Reads ECG recording in HES format
%   read_HES_header                      - Reads ECG header in HES format
%   read_ishne                           - (Internal) Reads ECG recordings in Mortara format
%   read_ishne_ann                       - Reads ECG annotations from ISHNE format
%   read_ishne_header                    - Reads ECG header from ISHNE format
%   read_Mortara                         - (Internal) Reads ECG recordings in Mortara format
%   read_Mortara_header                  - Reads ECG header in Mortara format. 
%   readheader                           - (Internal) function reads the header of signal files
%   reportECG                            - (Internal) function reads the header of signal files
%   rotateticklabel                      - rotates tick labels
%   rowvec                               - (Internal) Reshape a matrix into a row vector
%   Seconds2HMS                          - (Internal) Create a string of hours mins and seconds based on data in seconds
%   set_a_linespec                       - (Internal) Set a series of properties to a line handle
%   set_rand_linespec                    - (Internal) Set random properties a line handle
%   soft_intersect                       - (Internal) Intersection of two sets with tolerance
%   soft_range_conversion                - (Internal) Convert an input range to an output range with a soft function
%   soft_set_difference                  - (Internal) Set difference with tolerance
%   tablas_y_constantes                  - (Internal) Constants for the HES format
%   TaskPartition                        - (Internal) Generate a PIDs work list
%   text_arrow                           - (Internal) Plot an arrow with text in a graphic
%   text_line                            - (Internal) Plot a line with text in a graphic
%   trim_ECG_header                      - (Internal) Trim a header info struct to a subset of signals







%   WFDB_command_prefix                  - System commands to initialize the WFDB toolbox
%   woody_method                         - (Internal) Woody algorithm for heartbeat allignment
%   writeannot                           - (Internal) Write annotation files for biomedical signals in MIT Format. (MEX file)
%   writeheader                          - (Internal) Write an ECG header in MIT format
%   AUC_calc                             - Compute area under the ROC curve (AUC).
%   ECGtask_Delineation_corrector        - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_arbitrary_function           - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_classification_features_calc - ECGtask for ECGwrapper (for Matlab)
%   GetBestQRSdetections                 - Fetch the best QRS detections from an ECGtask_QRSdetections object
%   HasAdminPrivs                        - Checks administrator privileges 
%   default_concatenate_function         - Description:
%   default_finish_function              - Description:
%   list_all_ECGtask                     - List al ECGtask availables
%   matformat_definitions                - (Internal) A definition or header file, for names allowed for signals, header and annotations included in MAT format files
%   ppval                                - Evaluate piecewise polynomial.
%   read_Mortara_format                  - Reads ECG recordings in Mortara format. 
%   read_ishne_format                    - Reads ECG recording in ISHNE format
%   similarity_calculation               - (Internal) Pattern matching function to be used in an arbitrary task
%   sys_cmd_separation_string            - String to issue multiline system commands
%   sys_command_strings                  - Strings to execute typical I/O commands via system calls
