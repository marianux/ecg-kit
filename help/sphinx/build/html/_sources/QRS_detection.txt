
QRS detection
=============

This document describes how to perform QRS detection.

Description
-----------

Heartbeat detection is probably one of the first and most important
tasks when you process cardiovascular recordings. The ECGkit has several
algorithms implemented:

-  `Wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__
-  `gqrs <http://www.physionet.org/physiotools/wag/gqrs-1.htm>`__
-  `wqrs <http://www.physionet.org/physiotools/wag/wqrs-1.htm>`__
-  `Pan and
   Tompkins <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&arnumber=4122029>`__ 
-  `ECGpuwave <http://www.physionet.org/physiotools/ecgpuwave/>`__
-  `EP limited <http://www.eplimited.com/confirmation.htm>`__ 
-  Aristotle (not distributed with the kit)

You can use any or all the algorithms as you will see below or you can
even add your own algorithms if you follow an easy interface, as
described :ref:`below <Adding_a_custom_detection_algorithm>`. The results
are stored in a single file, that you can use to perform other
subsequent tasks, such as :doc:`ECG delineation <ECGdelineation>` or
even the visual :doc:`inspection <plot_ecg_strip>` or
:doc:`correction <QRScorrector>` of the algorithms results. For a quick
reference about heartbeat detection you may want to check this
:ref:`example <QRS_automatic_detection>`

 

Input Arguments
---------------

The properties that this task uses are the following:

``progress_handle`` — Used to track the progress within your function. ``[] (default)``

	progress\_handle, is a handle to a :doc:`progress\_bar <progress_bar>`
	object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data. ``tempdir() (default)``

	Full path to a directory with write privileges.

``detectors`` — The QRS detection algorithms to use. ``'all-detectors' (default)``

	This property controls which algorithms are used. A cell string or char array with any of the following names

	- all-detectors
	- `wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__
	- `pantom <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&arnumber=4122029>`__
	- aristotle
	- `gqrs <http://www.physionet.org/physiotools/wag/gqrs-1.htm>`__
	- `sqrs <http://www.physionet.org/physiotools/wag/sqrs-1.htm>`__
	- `wqrs <http://www.physionet.org/physiotools/wag/wqrs-1.htm>`__
	- `ecgpuwave <http://www.physionet.org/physiotools/ecgpuwave/>`__
	- `epltdqrs1 or epltdqrs2 <http://www.eplimited.com/confirmation.htm>`__

``only_ECG_leads`` — Process only ECG signals. ``true (default)``

	Boolean value. Find out which signals are ECG based on their header
	description.

``gqrs_config_filename`` — A configuration filename for the gqrs algorithm. ``[] (default)``

	A full filename with the configuration for the gqrs algorithm. See the
	algorithm `web <http://www.physionet.org/physiotools/wag/gqrs-1.htm>`__
	page for details.

``detection_threshold`` — A threshold to control the sensitivity of the detector. ``1 (default)``
	
	Use higher values to reduce false detections, or lower values to reduce the number of missed beats.

.. _payload_prop:

``payload`` — An arbitrary format variable to be passed to your user-defined algorithm. ``[] (default)``

	This variable can be useful for passing data to your own function, in addition to the interface described
	:ref:`below <Adding_a_custom_detection_algorithm>`.

``CalculatePerformance`` — Calculate algorithm performances based on gold standard reference detections. ``false (default)``

	Boolean value. Calculate the algorithm performance based on the
	reference annotations found by the ECGwrapper object in the same folder
	where the signals are. This reference annotations are loaded, if
	detected, in the ECG\_annotations property.
 
``bRecalcQualityAndPerformance`` — Recalculate algorithm performances without performing heartbeat detection. ``false (default)``

	Boolean value. Recalculate algorithm performances without performing 
	heartbeat detection. This feature is useful for BIG jobs where only 
	a small change in performance calculation methodology was changed.
	. The cached result is passed to the task via the payload property 
	as in the example:
 
.. code::

	% If only recaclulate performance is needed.
	ECGw.ECGtaskHandle.bRecalcQualityAndPerformance = true;    
	% in order to avoid skiping the task
	ECGw.cacheResults = false; 
	%payload has the current detections
	cached_filenames = ECGw.GetCahchedFileName('QRS_detection');
	ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
 
``bRecalculateNewDetections`` — Calculate only heartbeat detections not performed before. ``false (default)``

	Boolean value. The heartbeat detection is only performed in the 
	algorithms not executed in a previous cached result. The cached 
	result is passed to the task via the payload property as in the 
	example:
 
.. code::

    % this is needed in order to recalculate tasks.
	% Useful if a new detector is added 
	ECGw.ECGtaskHandle.bRecalculateNewDetections = true;    
	% in order to avoid skiping the task
	ECGw.cacheResults = false; 
	%payload has the current detections
	cached_filenames = ECGw.GetCahchedFileName('QRS_detection');
	ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
			
			
.. _Adding_a_custom_detection_algorithm:

Adding a custom detection algorithm
-----------------------------------

Adding your own QRS detectors to the kit is very simple. Ensure that
your function implements this interface:

.. code::

    function [positions_single_lead, position_multilead] = 
			
			your_QRS_detector( ECG_matrix, ECG_header, progress_handle, payload_in) 
                            

where the arguments are:

	**ECG\_matrix**, is a matrix size ``[ECG\_header.nsamp ECG\_header.nsig]``

	**ECG\_header**, is a struct with info about the ECG signal, see :ref:`ECG header <ECG_header_description>` 
	for details.
	
	**progress\_handle**, is a handle to a `progress\_bar <progress_bar.htm>`__
	object, that can be used to track the progress within your function.

	**payload\_in**, is a user variable, of arbitrary format, allowed to be sent
	to your function. It is sent via the :ref:`payload property <payload_prop>` 
	of this class, for example:

.. code::
	
	% One variable
	this_ECG_wrapper.ECGtaskHandle.payload = your_variable;
	
	% Several variables with a cell container
	this_ECG_wrapper.ECGtaskHandle.payload = {your_var1 your_var2};
	
	% Or the result of a previous task, in this case QRS manual correction (if available)
	% or the automatic detection if not.
	cached_filenames = this_ECG_wrapper.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
	this_ECG_wrapper.ECGtaskHandle.payload = load(cached_filenames);

and the output of your function must be:

	**positions\_single\_lead**, a cell array size ECG\_header.nsig with the QRS
	sample locations found in each lead.

	**position\_multilead**, a numeric vector with the QRS locations calculated
	using multilead rules.
 

Examples
--------

Create the ECGtask\_QRS\_detection object.

.. code::

    % with the task name
    ECG_w.ECGtaskHandle = 'QRS_detection';
    
	% or create an specific handle to have more control
    ECGt_QRSd = ECGtask_QRS_detection();

and then you are ready to set the algorithms to use. In the following
example you have several possible set-ups.

.. code::

	% select an specific algorithm. Default: Run all detectors
	ECGt_QRSd.detectors = 'wavedet'; % Wavedet algorithm based on
	ECGt_QRSd.detectors = 'pantom';  % Pan-Tompkins alg.
	ECGt_QRSd.detectors = 'gqrs';    % WFDB gqrs algorithm.
	% Example of how you can add your own QRS detector.
	ECGt_QRSd.detectors = 'user:example_worst_ever_QRS_detector';    
	% "your_QRS_detector_func_name" can be your own detector.
	ECGt_QRSd.detectors = 'user:your_QRS_detector_func_name';    
	ECGt_QRSd.detectors = {'wavedet' 'gqrs' 'user:example_worst_ever_QRS_detector'};
                            

Finally set the task to the wrapper object, and execute the task.

.. code::

    ECG_w.ECGtaskHandle= ECGt_QRSd; % set the ECG task
    ECG_w.Run();

You can check the result of this task, with either the :doc:`detection
corrector <QRScorrector>` or the :doc:`visualization
functions <plot_ecg_strip>`.

Also check this :ref:`example <QRS_automatic_detection>` for
further information.

.. _QRS_det_result_format:

Results format
--------------
 
The result file will have ``ECG_header.nsig x algorithms_used`` variables, which can later be recovered 
as a ``struct`` variable, with fields named according to ``[ 'algorithm_name' '_' 'lead_name' ]``. Each
of this fields is a ``struct`` itself with a single field called ``time``, where the actual QRS detections are.
In addition, another ``struct`` variable called ``series_quality`` is stored in order to provide a quality metric of 
the detections created. This metric is found in the ``ratios`` field, a higher ratio means better detections.
Each ratio corresponds with a name in the ``AnnNames`` field.


More About
----------

Here are some external references about heartbeat detection:

-  `PhysioNet/Computing in Cardiology Challenge
   2014 <http://physionet.org/challenge/2014/>`__
-  `Physionet <http://www.physionet.org/>`__
-  A video demo in `Youtube <https://www.youtube.com/watch?v=QrM-aYANUns&index=2&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7>`__

See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`ECG delineation <ECGdelineation>` \| :doc:`examples <examples>`

