
Arbitrary tasks
===============

This document describes how to use arbitrary tasks with the ECGkit.


Description
-----------

Sometimes the task you need to perform on ECG signals is too simple to
develop a new ECGtask, like computing some statistics, or apply a linear
filter, or any type of transformation you may need to perform to the
signal. For those cases you may found arbitrary tasks useful.

 

Input Arguments
---------------

The properties that this task uses are the following:

``progress_handle`` — Used to track the progress within your function. ``[] (default)``

	progress\_handle, is a handle to a :doc:`progress\_bar <progress_bar>`
	object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data. ``tempdir() (default)``

	Full path to a directory with write privileges.

``only_ECG_leads`` — Process only ECG signals ``true (default)`` 

	Boolean value. Find out which signals are ECG based on their ``ECG_header.desc`` description.

``payload`` — An arbitrary format variable to be passed to your function. ``[] (default)`` 

	This variable can be useful for passing data to your own function.

``signal_payload`` — Consider the result of your arbitrary function as a signal. ``false (default)`` 

	Boolean value that indicates the ECGwrapper to produce a signal instead of a result payload.

``lead_idx`` — The signal indexes that your function will affect. ``[] (default)`` 

	A positive integer array with values from 1 to ``ECG_header.nsig``.

``function_pointer`` — The pointer to your arbitrary function. ``[] (default)`` 

Your function must follow this prototype:

.. code::

    function result = your_function( ECG_matrix, ECG_header, progress_handle, payload_in)  

							
where the arguments are:

	**ECG\_matrix**, is a matrix size ``[ECG\_header.nsamp ECG\_header.nsig]``

	**ECG\_header**, is a struct with info about the ECG signal, such as:

		- *freq*, is the sampling frequency of ECG\_matrix signal.

		- *desc*, description strings about each of the leads/signals.

		- *nsamp* is the number of samples of ECG\_matrix.

		- *nsig* is the amount of leads or signals of ECG\_matrix.

		- *gain* is a vector of [nsig × 1] with the gain of each lead ( ADCsamples / μV ).

		- *adczero* is a vector of [nsig × 1] with the offset of each lead in ADC samples.
		
		and others described in the `Physionet header <http://www.physionet.org/physiotools/wag/header-5.htm>`__.

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


and the output of your function must be a struct variable **result**, or
if it is a signal, ensure to make ``true`` the ``signal_payload`` property.


Examples
--------

This example is used in the QRScorrector function to perform
template-matching on an ECGwrapper (arbitrary big recording) object.

.. code::

	aux_w = ECGwrapper('recording_name', 'your_path/recname');
	aux_w.ECGtaskHandle = 'arbitrary_function';

	% This is in case you want always to recalculate results, no caching
	aux_w.cacheResults = false;

	% Use first and third columns-signals
	aux_w.ECGtaskHandle.lead_idx = [1 3];

	% Produce a signal as a result
	aux_w.ECGtaskHandle.signal_payload = true;

	% Add a user-string to identify the run
	aux_w.ECGtaskHandle.user_string = ['similarity_calc_for_lead_' num2str(sort(lead_idx)) ];

	% add your function pointer
	aux_w.ECGtaskHandle.function_pointer = @similarity_calculation;

	% and any data your function may need.
	aux_w.ECGtaskHandle.payload = pattern2detect;
	% and you are ready to go !
	aux_w.Run
                            

.. _arbitrary_result_format:

Results format
--------------
 
The result file will have ``ECG_header.nsig x algorithms_used`` variables, which can later be recovered 
as a ``struct`` variable, with fields named according to ``[ 'algorithm_name' '_' 'lead_name' ]``. Each
of this fields is a ``struct`` itself with a single field called ``time``, where the actual QRS detections are.
In addition, another ``struct`` variable called ``series_quality`` is stored in order to provide a quality metric of 
the detections created. This metric is found in the ``ratios`` field, a higher ratio means better detections.
Each ratio corresponds with a name in the ``AnnNames`` field.
							
							
							
See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`QRS detection <QRS_detection>` \| :doc:`ECG delineation <ECGdelineation>` \| :doc:`examples <examples>`
