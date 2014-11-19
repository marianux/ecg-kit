% COMMON
%
% Files
%   ADC2realunits                       - 
%   ADC2units                           - 
%   addpath_if_not_added                - 
%   adjust_string                       - 
%   Annotation_process                  - 
%   AnnotationFilterConvert             - 
%   arrow                               - interface
%   autovec_calculation                 - 
%   autovec_calculation_robust          - 
%   bandpass_filter_design              - UNTITLED Returns a discrete-time filter System object.
%   BaselineWanderRemovalMedian         - Script creado para el taller de Matlab sobre procesamiento de ECG
%   BaselineWanderRemovalSplines        - Description:
%   bxb                                 - constants
%   calc_btime                          - 
%   calc_co_ocurrences                  - 
%   calc_correlation_gain               - 
%   CalcRRserieQuality                  - function [q_measure noise_power] = CalcRRserieQuality(signal, heasig, references, noise_power)
%   CalcRRserieRatio                    - Description: 
%   calculate_artificial_QRS_detections - 
%   calculateSeriesQuality              - estimate quality of QRS detections performed
%   cluster_data_with_EM_clust          - 
%   colvec                              - Reshape the input into a column vector
%   combine_anns                        - lreferences = length(time_serie);
%   ConcatenateQRSdetectionPayloads     - 
%   cprintf                             - displays styled formatted text in the Command Window
%   d2p                                 - Identifies appropriate sigma's to get kk NNs up to some tolerance 
%   DelayedCovMat                       - Internal function for calculating covariance matrix in the vicinity of
%   deNaN_dataset                       - 
%   design_downsample_filter            - Returns a discrete-time filter object.
%   detect_QRS_complexes                - 
%   disp_option_enumeration             - This function format a string for fprintf like functions in order to
%   disp_string_framed                  - 
%   disp_string_title                   - 
%   DisplayConfusionMatrix              - 
%   DisplayResults                      - 
%   ds2nfu                              - Convert data space units into normalized figure units. 
%   ECG_avg                             - 
%   ECGformat                           - 
%   ECGtask                             - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_do_nothing                  - Null ECGtask (for Matlab)
%   ECGtask_ECG_delineation             - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_ECG_delineation_corrector   - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_heartbeat_classifier        - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PCA_proj_basis              - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PPG_ABP_corrector           - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_PPG_ABP_detector            - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_corrector               - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_detection               - ECGtask for ECGwrapper (for Matlab)
%   ECGtask_QRS_detections_post_process - ECGtask for ECGwrapper (for Matlab)
%   ECGwrapper                          - ECGwrapper for Matlab
%   example_worst_ever_ECG_delineator   - Example of user-created QRS detector
%   example_worst_ever_QRS_detector     - Example of user-created QRS detector
%   exist_distributed_file              - 
%   get_ECG_idx_from_header             - 
%   get_PPG_ABP_idx_from_header         - 
%   getAnnNames                         - 
%   GetFunctionInvocation               - 
%   GTHTMLtable                         - GTHTMLtable - Generate an HTML page with a table of a matrix.
%   init_ghostscript                    - 
%   init_WFDB_library                   - 
%   isAHAformat                         - Check if a recording is in ISHNE format.
%   isHESformat                         - Check if a recording is in HES format.
%   isISHNEformat                       - Check if a recording is in ISHNE format.
%   isMatlab                            - 
%   isOctave                            - 
%   IterPartition                       - no es recomendable hacerlo mas grande de cant_recs la particion.
%   logit_function                      - 
%   matrix2positions                    - 
%   max_index                           - 
%   maximize                            - Size a window to fill the entire screen.
%   MedianFilt                          - 
%   modmax                              - Function which returns the indexes of vector x in which there are
%   my_colormap                         - 
%   my_ppval                            - PPVAL  Evaluate piecewise polynomial.
%   myzerocros                          - Returns the index of the input vector in which the first zero crossing is located.
%   nanmeda                             - 
%   pack_signal                         - Description:
%   parse_pids                          - 
%   PCA_window                          - [y,W,A] = PCA_window(x, peaks, win_size)
%   PeakDetection2                      - peaks = PeakDetection2(x,fs,wlen,fp1,fp2,th,flag),
%   PiCA                                - [y,W,A] = PiCA(x,peaks1,peaks2)
%   plot_auc                            - 
%   plot_ecg_heartbeat                  - obsolete, use plot_ecg_strip.m
%   plot_ecg_mosaic                     - Description: 
%   plot_ecg_strip                      - This function plot ECG signals, eventually with annotation marks such as
%   plot_n_QRS                          - 
%   plot_roc                            - 
%   PlotGlobalWaveMarks                 - 
%   PlotWaveMarks                       - 
%   pos2mit                             - output = pos2mit(ann_samples,ann_label,ann_lead,ann_path,ann_filename)
%   positions2matrix                    - 
%   PrctileFilt                         - 
%   progress_bar                        - Description:
%   progress_bar_ex                     - start of algorithm
%   QRScorrector                        - Description: 
%   qs_filter_design                    - Prototype:
%   qs_wt                               - Prototype:
%   queue2read                          - 
%   queue2write                         - 
%   rand_linespec                       - 
%   read_310_format                     - Reads ECG recordings in 310 format. Implements the documentation
%   read_311_format                     - Reads ECG recordings in 311 format. Implements the documentation
%   read_AHA_ann                        - Reads ECG annotations in AHA format. Implements the documentation
%   read_AHA_format                     - Reads ECG recordings in AHA format. Implements the documentation
%   read_AHA_header                     - Reads ECG headers in AHA format. Implements the documentation
%   read_ECG                            - 
%   read_HES_ann                        - La traduccion se termina haciendo en LoadDatabaseCustomization.m
%   read_HES_format                     - Reads ECG recordings in HES (Biosigna) format. Implements the documentation
%   read_HES_header                     - Reads ECG headers in HES (Biosigna's) format. Implements the documentation
%   read_ishne                          - Reads ECG recordings in ISHNE format from THEW databases. Implements the documentation available in:
%   read_ishne_ann                      - Reads ECG annotations in ISHNE format from THEW databases. Implements the documentation available in:
%   read_ishne_header                   - Reads ECG headers in ISHNE format from THEW databases. Implements the documentation available in:
%   read_Mortara                        - Reads ECG recordings in Mortara format. 
%   read_Mortara_header                 - Reads ECG header in Mortara format. 
%   readheader                          - function reads the header of DB signal files
%   reportECG                           - This function creates a report of the signals handled by the ECG wrapper
%   rotateticklabel                     - rotates tick labels
%   rowvec                              - Reshape a matrix into a row vector
%   Seconds2HMS                         - 
%   set_a_linespec                      - 
%   set_rand_linespec                   - 
%   soft_intersect                      - Find the intersection between index sequence val1-2 within a windows
%   soft_range_conversion               - Sigmoid function mapping an input range to an output range:
%   soft_set_difference                 - Find the set difference between index sequence idx1-2 within a windows
%   tablas_y_constantes                 - 
%   TaskPartition                       - no es recomendable hacerlo mas grande de cant_recs la particion.
%   text_arrow                          - interface
%   text_line                           - interface
%   trim_ECG_header                     - UNTITLED Summary of this function goes here
%   tsne                                - Performs symmetric t-SNE on dataset this_X
%   tsne_d                              - Performs symmetric t-SNE on the pairwise Euclidean distance matrix D
%   tsne_p                              - Performs symmetric t-SNE on affinity matrix P
%   unqueue2read                        - 
%   unqueue2write                       - 
%   wavedet_interface_heartbeats        - 
%   wavedet_QRS_detection_mix           - attemp to build a better detection from single-lead detections.
%   WFDB_command_prefix                 - string to change the working directory.
%   woody_method                        - 
%   writeannot                          - Matlab mex-file
%   writeheader                         - WRITHEAD function writes header file for signal data struct in directory header_path
