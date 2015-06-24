
ECG delineation
===============

This document describes how to perform automatic delineation or wave
segmentation on ECG signals.

Description
-----------

Automatic wave segmentation or delineation is exclusively performed by 
`wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__ algorithm.

Input Arguments
---------------

The properties that this task uses are the following:

``progress_handle`` — Used to track the progress within your function. ``[] (default)``

	progress\_handle, is a handle to a :doc:`progress\_bar <progress_bar>`
	object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data. ``tempdir() (default)``

	Full path to a directory with write privileges.

``delineators`` — The ECG delineation algorithms to use ``'all-delineators' (default)``

	This property controls which algorithms are used. A cell string or char with any of the following names

	- *'all-delineators'*

	- `'wavedet' <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__

 

``only_ECG_leads`` — Process only ECG signals ``true (default)`` 

	Boolean value. Find out which signals are ECG based on their ``ECG_header.desc`` 
	description.

``wavedet_config`` — A structure for configuring wavedet algorithm. ``[] (default)`` 

	Undocumented yet, use it only if you know what you are doing.

``payload`` — An arbitrary format variable. ``[] (default)`` 

	This variable can be useful for passing data to your own delineation function
	(described :ref:`below <Adding_a_custom_delineation_algorithm>`) or to
	provide visually audited QRS detections to the delineation algorithm.

.. _Adding_a_custom_delineation_algorithm:

Adding a custom delineation algorithm
-------------------------------------

Adding your own delineator to the kit is very simple. Ensure that your
function implements this interface:

.. code::

    function [positions_single_lead, position_multilead] = 
	
		your_ECG_delineation( ECG_matrix, ECG_header, progress_handle, payload_in)  
                            

where the arguments are:

	**ECG\_matrix**, is a matrix size [ECG\_header.nsamp ECG\_header.nsig]

	.. _ECG_header_description:
	
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

and the output of your function must be:

.. _delineation_struct:

	**positions\_single\_lead**, is an **structure array** of ``ECG_header.nsig`` elements with *at least* the following wave fiducial points as fields:
	
	- ``'Pon'`` P wave onset
	- ``'P'`` P wave peak
	- ``'Poff'`` P wave offset
	- ``'QRSon'`` QRS complex onset
	- ``'qrs'`` QRS fiducial point, obtained from QRS detection.
	- ``'Q'`` Q wave peak
	- ``'R'`` R wave peak
	- ``'S'`` S wave peak
	- ``'QRSoff'`` QRS complex offset
	- ``'Ton'`` T wave onset
	- ``'T'`` T wave peak
	- ``'Toff'`` T wave offset

	**position\_multilead**, is a single structure with *at least* the wave fiducial points described above.
	This delineation is commonly calculated from the single lead delineations, in order to obtain a unique wave 
	fiducial point per heartbeat.
	

Examples
--------

Create the *ECGtask\_ECG\_delineation* object.

.. code::

    % with the task name
        ECG_w.ECGtaskHandle = 'ECG_delineation';
    % or create an specific handle to have more control
        ECGt = ECGtask_ECG_delineation();

and then you are ready to set the algorithms to use. In the following
example you have several possible set-ups.

.. code::

    % select an specific algorithm. Default: Run all detectors
            ECGt.delineators = 'wavedet'; % Wavedet algorithm based on
            % "your_delineator_func_name" can be your own delineator.
			ECGt.delineators = 'user:your_delineator_func_name';    
            ECGt.delineators = {'wavedet' 'user:your_delineator_func_name'};
                            

Finally set the task to the wrapper object, and execute the task.

.. code::

            ECG_w.ECGtaskHandle= ECGt; % set the ECG task
            ECG_w.Run();

You can check the result of this task, with either the :doc:`delineator
corrector <ECG_delineation_corrector>` or the :doc:`visualization
functions <plot_ecg_strip>`.

Also check this :ref:`example <ECG_automatic_delineation>`
for further information.

.. _Delineation_result_format:

Results format
--------------
 
The result file will have a ``struct`` variable with the name of the algorithm (only *wavedet* at the time of 
writing this). Inside this, it will contain one :ref:`delineation struct <delineation_struct>` per ECG lead 
in the ``ECG_header.desc`` field, plus another called ``multilead`` which is a delineation accounting with the 
information present in all leads.


More About
----------

This publication describes the
`wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__
algorithm.

See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`QRS detection <QRS_detection>` \| :doc:`examples <examples>`

