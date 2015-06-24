
Other functions
===============

Several low-level functions that are located in ``your_path\ecg-kit\common\`` 
but are not yet well documented, tested or integrated with other parts of the kit.

General functions
-----------------

- addpath_if_not_added                 - (Internal) Add the path only if not was already added.
- colvec                               - (Internal) Reshape the input into a column vector
- rowvec                               - (Internal) Reshape a matrix into a row vector
- init_ghostscript                     - (Internal) Init environment variables for using ghostscript
- init_WFDB_library                    - (Internal) Init environment variables for using WFDB toolbox
- isMatlab                             - (Internal) Check if the kit is running on Matlab
- isOctave                             - (Internal) Check if the kit is running on Octave
- exist_distributed_file               - (Internal) Check the existence of a file in a distributed (slow) filesystem
- GetFunctionInvocation                - (Internal) Create a string with the invocation of a function
- max_index                            - (Internal) Index of the maximum element in a vector
- modmax                               - (Internal) Find modulus maxima in a signal
- myzerocros                           - (Internal) Detect zero-crosses in a signal
- soft_intersect                       - (Internal) Intersection of two sets with tolerance
- soft_range_conversion                - (Internal) Convert an input range to an output range with a soft function
- soft_set_difference                  - (Internal) Set difference with tolerance
- parse_pids                           - (Internal) Identify how many PIDs are in total and which is this PID, based on a string formatted this_pid/cant_pids
- getAnnNames                          - (Internal) Get names of annotations from annotation structure
- matrix2positions                     - (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
- positions2matrix                     - (Internal) Convert matrix of ECG wave annotations to a struct position format, used in wavedet algorithm
- pack_signal                          - (Internal) Example of user-created QRS detector
- :doc:`progress_bar <progress_bar>`   - (Internal) A progress bar class for showing evolution of a process to users
- progress_bar_ex                      - A progress bar class example
- TaskPartition                        - (Internal) Generate a PIDs work list
- trim_ECG_header                      - (Internal) Trim a header info struct to a subset of signals
- WFDB_command_prefix                  - System commands to initialize the WFDB toolbox
- HasAdminPrivs                        - Checks administrator privileges 
- sys_cmd_separation_string            - String to issue multiline system commands
- sys_command_strings                  - Strings to execute typical I/O commands via system calls


Strings related
----------------

- adjust_string                        - (Internal) Works with strings to center, trim and justify to a certain string width
- calc_btime                           - (Internal) Creates a string with the base time
- disp_string_framed                   - (Internal) Display a message framed to a string
- disp_string_title                    - (Internal) Display a message framed to a string
- disp_option_enumeration              - (Internal) Display an enumeration of options to a string
- DisplayConfusionMatrix               - (Internal) Pretty display in Screen the confusion matrix.
- DisplayResults                       - (Internal) Pretty-Display results of a classification experiment
- Seconds2HMS                          - (Internal) Create a string of hours mins and seconds based on data in seconds


Graphics related
----------------

- arrow                                - (Internal) Creates arrows in grapics
- ds2nfu                               - Convert data space units into normalized figure units. 
- plot_auc                             - (Internal) Plot the area under the ROC curve
- plot_ecg_heartbeat                   - (Internal) obsolete, use plot_ecg_strip
- plot_ecg_mosaic                      - Plots multidimentional signal in mosaic style
- plot_ecg_strip                       - Plots and interact with ECG signal
- plot_roc                             - (Internal) Plot the ROC curve
- PlotGlobalWaveMarks                  - (Internal) Internal function of plot_ecg_strip
- PlotWaveMarks                        - (Internal) Internal function of plot_ecg_strip
- my_colormap                          - (Internal) Create a colormap for signal visualization
- maximize                             - Size a window to fill the entire screen.
- rand_linespec                        - (Internal) Example of user-created QRS detector
- rotateticklabel                      - rotates tick labels
- set_a_linespec                       - (Internal) Set a series of properties to a line handle
- set_rand_linespec                    - (Internal) Set random properties a line handle
- text_arrow                           - (Internal) Plot an arrow with text in a graphic
- text_line                            - (Internal) Plot a line with text in a graphic


Signal processing / statistical methods 
---------------------------------------

- autovec_calculation                  - (Internal) Calculate eigenvalues and vectors
- autovec_calculation_robust           - (Internal) Calculate eigenvalues and vectors using robust covariance estimation
- BaselineWanderRemovalMedian          - Remove baseline wandering with the median estimation method
- BaselineWanderRemovalSplines         - Remove baseline wandering with the median estimation method
- bxb                                  - (Internal) Compares two heartbeat series and produce the confusion matrix as result
- calc_co_ocurrences                   - (Internal) Calculate the heartbeats co-ocurrences
- calc_correlation_gain                - (Internal) Add the path only if not was already added.
- CalcRRserieQuality                   - (Obsolete) Estimate the quality of QRS complex detections 
- CalcRRserieRatio                     - Estimate the quality of QRS complex detections 
- calculateSeriesQuality               - Estimate the quality of QRS complex detections 
- MedianFilt                           - (Internal) Mean/Median filtering
- cluster_data_with_EM_clust           - (Internal) Cluster data with expectation-maximization algorithm
- DelayedCovMat                        - (Internal) Delayed covariance matrix calculation
- deNaN_dataset                        - (Internal) Replace NaN from PRdatasets
- design_downsample_filter             - (Internal) Design a filter to downsample signals prior printing to a report
- bandpass_filter_design               - MATLAB Code
- PeakDetection2                       - peaks = PeakDetection2(x,fs,wlen,fp1,fp2,th,flag),
- PiCA                                 - [y,W,A] = PiCA(x,peaks1,peaks2)
- PrctileFilt                          - (Internal) Arbitrary percentile filtering
- logit_function                       - (Internal) Logit function
- my_ppval                             - PPVAL  Evaluate piecewise polynomial.
- qs_filter_design                     - (Internal) Design the wavelet decomposition filters for wavedet algorithm
- qs_wt                                - (Internal) Calculates the wavelet transform
- nanmeda                              - (Internal) Calculate the median of absolute deviations from the median (MEDA)
- woody_method                         - (Internal) Woody algorithm for heartbeat allignment
- AUC_calc                             - Compute area under the ROC curve (AUC).
- ppval                                - Evaluate piecewise polynomial.
- similarity_calculation               - (Internal) Pattern matching function to be used in an arbitrary task
- combine_anns                         - (Internal) Create new QRS detections based on other lead/algorithms detections



Tasks
-----

- ECGtask                              - Defines the class interface for the ECGtask derived classes
- ECGtask_do_nothing                   - Null ECGtask (for Matlab)
- ECGtask_ECG_delineation              - ECGtask for ECGwrapper (for Matlab)
- ECGtask_ECG_delineation_corrector    - ECGtask for ECGwrapper (for Matlab)
- ECGtask_heartbeat_classifier         - ECGtask for ECGwrapper (for Matlab)
- ECGtask_PCA_proj_basis               - ECGtask for ECGwrapper (for Matlab)
- ECGtask_PPG_ABP_corrector            - ECGtask for ECGwrapper (for Matlab)
- ECGtask_PPG_ABP_detector             - ECGtask for ECGwrapper (for Matlab)
- ECGtask_QRS_corrector                - ECGtask for ECGwrapper (for Matlab)
- ECGtask_QRS_detection                - ECGtask for ECGwrapper (for Matlab)
- ECGtask_QRS_detections_post_process  - ECGtask for ECGwrapper (for Matlab)
- ECGtask_Delineation_corrector        - ECGtask for ECGwrapper (for Matlab)
- ECGtask_arbitrary_function           - ECGtask for ECGwrapper (for Matlab)
- ECGtask_classification_features_calc - ECGtask for ECGwrapper (for Matlab)
- example_worst_ever_ECG_delineator    - (Internal) Example of user-created QRS detector
- example_worst_ever_QRS_detector      - (Internal) Example of user-created QRS detector
- reportECG                            - (Internal) function reads the header of signal files
- QRScorrector                         - (Internal) GUI for correcting QRS detections
- GetBestQRSdetections                 - Fetch the best QRS detections from an ECGtask_QRSdetections object
- list_all_ECGtask                     - List al ECGtask availables


Functions from other projects
-----------------------------

- GTHTMLtable                          - GTHTMLtable - Generate an HTML page with a table of a matrix.
- cprintf                              - displays styled formatted text in the Command Window


I/O signals
-----------

- Annotation_process                   - Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
- AnnotationFilterConvert              - Convert heartbeat type of annotation from valid ECG formats to EC57 AAMI
- ADC2realunits                        - Convert adimentional sample values to real units
- ADC2units                            - Convert adimentional sample values to target voltage units
- ECGwrapper                           - Allow acces to ECG recordings of arbitrary format and length.
- ECGformat                            - Gets the format of an ECG recording filename
- get_ECG_idx_from_header              - (Internal) Guess ECG signals indexes in a multimodal recording
- get_PPG_ABP_idx_from_header          - (Internal) Guess PPG/ABP signals indexes in a multimodal recording
- isAHAformat                          - (Internal) Check if a recording is in ISHNE format.
- isHESformat                          - (Internal) Check if a recording is in HES format.
- isISHNEformat                        - (Internal) Check if a recording is in ISHNE format.
- matformat_definitions                - (Internal) A definition or header file, for names allowed for signals, header and annotations included in MAT format files
- read_310_format                      - (Internal) Read the MIT 310 format
- read_311_format                      - (Internal) Read the MIT 311 format
- read_AHA_ann                         - Reads ECG annotations in AHA format
- read_AHA_format                      - Reads ECG recording in AHA format
- read_AHA_header                      - Reads ECG header in AHA format
- read_ECG                             - Reads an ECG recording
- read_HES_ann                         - Reads ECG annotations in HES format
- read_HES_format                      - Reads ECG recording in HES format
- read_HES_header                      - Reads ECG header in HES format
- read_ishne                           - (Internal) Reads ECG recordings in Mortara format
- read_ishne_ann                       - Reads ECG annotations from ISHNE format
- read_ishne_header                    - Reads ECG header from ISHNE format
- read_Mortara                         - (Internal) Reads ECG recordings in Mortara format
- read_Mortara_header                  - Reads ECG header in Mortara format. 
- read_Mortara_format                  - Reads ECG recordings in Mortara format. 
- read_ishne_format                    - Reads ECG recording in ISHNE format
- readheader                           - (Internal) function reads the header of signal files
- tablas_y_constantes                  - (Internal) Constants for the HES format
- writeannot                           - (Internal) Write annotation files for biomedical signals in MIT Format. (MEX file)
- writeheader                          - (Internal) Write an ECG header in MIT format
- ConcatenateQRSdetectionPayloads      - (Internal) Concatenate two payloads.
- default_concatenate_function         - Description:
- default_finish_function              - Description:


.. toctree::
   :titlesonly:
   :hidden:

   A progress bar <progress_bar>
   List installed ECGtask <list_all_ECGtask>
   More docs soon ... <progress_bar>

