
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

	**progress\_handle**, is a handle to a :doc:`progress\_bar <progress_bar>`
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

``finish_func_pointer`` — A pointer to your arbitrary finish function. ``@default_finish_function (default)`` 

A function that will operate over the whole result of your arbitrary function, after the payloads resulting of each 
iteration were concatenated. This is only used when the result of your ``function_pointer`` is **not** a signal 
(``signal_payload = false``). Your function must follow this prototype:

.. code::

	payload = your_finish_function(payload, ECG_header)

	
where the arguments are:

	**payload**, is the complete payload.

	**ECG\_header**, is a struct with info about the ECG signal, see above for reference.
	
and this function will change the payload variable as according to your needs and return it to the ECGwrapper object.
	
``concate_func_pointer`` — The pointer to your arbitrary concatenate function. ``@default_concatenate_function (default)`` 

A function that will concatenate or integrate the information produced in each part of your recording, when the result of 
your ``function_pointer`` is **not** a signal (``signal_payload = false``). Your function must follow this prototype:

.. code::

	payload = your_concatenate_function(plA, plB)

	
where the arguments are:

	**plA** and **plB** are the two payloads to concatenate
	
and this function will integrate or concatenate both payloads into the resulting payload. This resulting payload, will be
plA in the next iteration of concatenation. The ``default_concatenate_function`` just concatenate payloads:

.. code::

	% The default behavior of the concatenate function is to concatenate
	% payloads vertically or row-wise.
	if( isempty(plA) )
		payload = plB;
	else
		payload = [plA; plB];
	end


Examples
--------

1. Arbitrary task producing a **signal** as a result
 
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
 	aux_w.user_string = ['similarity_calc_for_lead_' num2str(sort(lead_idx)) ];
 
 	% add your function pointer
 	aux_w.ECGtaskHandle.function_pointer = @similarity_calculation;
 
 	% and any data your function may need.
 	aux_w.ECGtaskHandle.payload = pattern2detect;
 	% and you are ready to go !
 	aux_w.Run
 							
 	
 
2. Arbitrary task producing an **arbitrary result**
 	
 This is achieved by defining 3 properties (function handles) that perform:
 
 - The arbitrary task, which produces an arbitrary result ``function_pointer``
 
 - The concatenation of these results ``concate_func_pointer``
 
 - The final result calculation, when all results are concatenated. ``finish_func_pointer``
 
 The configuration of the ECGwrapper object is quite simple:
 
  	.. code::
		
		cd your_path\ecg-kit\examples
		ECGw = ECGwrapper( 'recording_name', 'your_path\ecg-kit\recordings\208')
		% no overlapp needed between signal partitions
		ECGw.partition_mode = 'ECG_contiguous';
		ECGw.ECGtaskHandle = 'arbitrary_function';
		ECGw.ECGtaskHandle.function_pointer = @my_mean;
		ECGw.ECGtaskHandle.concate_func_pointer = @my_concatenate_mean;
		ECGw.ECGtaskHandle.finish_func_pointer = @my_finish_mean;
		ECGw.Run
 
 
 
  The result is stored in a ``mat`` file.
 	
 	
 .. code-block:: none
 	
 	Description of the process:
 	 + Recording: d:\mariano\misc\ecg-kit\recordings\208.dat
 	 + Task name: arbitrary_function                             
 
 
 	##############
 	# Work done! #
 	##############
 
 
 	Results saved in
 	 + your_path\ecg-kit\recordings\208_arbitrary_function.mat	

 
 The arbitrary functions used to calculate the mean in an arbitrary large recording are:
 
 	- ``\ecg-kit\examples\my_mean.m`` In this function we only accumulate and count 
 	  the size of the accumulation.
 	  
 	.. code::
 
 		function result = my_mean(x)
 
 		result.the_sum = sum(x);
 		result.the_size = size(x,1);
 
 
 	- ``\ecg-kit\examples\my_concatenate_mean.m`` This function calculate the final 
 	  accumulation and counting.
 	  
 	.. code::
 
 		function payload = my_concatenate_mean(plA, plB)
 
 		if( isempty(plA) )
 			payload = plB;
 		else
 			payload.the_sum = plA.the_sum + plB.the_sum;
 			payload.the_size = plA.the_size + plB.the_size;
 		end
 
 
 	- ``\ecg-kit\examples\my_finish_mean.m`` In this function the mean calculation 
 	  is performed.
 	  
 	.. code::
 
 		function result_payload = my_finish_mean(payload, ECG_header)
 
 		result_payload.mean = payload.the_sum ./ payload.the_size;
 
 	
.. _arbitrary_result_format:

Results format
--------------
 
The format of the results depends on the ``signal_payload`` property, if it is a signal it will be in `MIT format`_.
Otherwise, the results depends on the user-defined output of 
							
							
							
See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`QRS detection <QRS_detection>` \| :doc:`ECG delineation <ECGdelineation>` \| :doc:`examples <examples>`

 
.. _`MIT format`: http://www.physionet.org/physiotools/wag/signal-5.htm
