
ABP/PPG peak detection
======================

This document describes how to perform automatic peak detection in
pulsatile signals.


Description
-----------

This task perform peak detection in pulsatile signals such as arterial
blood pressure (ABP) or plethysmographic (PPG). The task uses two
algorithms to achieve pulse detection, and, as in QRS detection task you
can choose to use any of them.

 

Input Arguments
---------------

The properties that the ECGtask\_PPG\_ABP\_detector class accepts are
described below. The usage of these properties is restricted to
low-level programming, you can use this task through the ECGwrapper as
is shown in the example below.

``progress_handle`` — Used to track the progress within your function. ``[] (default)``

	progress\_handle, is a handle to a :doc:`progress\_bar <progress_bar>`
	object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data. ``tempdir() (default)``

	A folder to store temporary data. Full path to a directory with write privileges.

``lead_config`` — Select which signals to process. ``'PPG-ABP-only' (default)`` 

	This property control on which signals the pulse detection algorithms will be applied. 
	A cell string or char with any of the following names:

	- '*all-leads*'. Process all leads.

	- '*PPG-ABP-only*'. Detect pulsatile signals based on their ``ECG_header.desc``
	  description variable.

	- '*User-defined-leads*'. Tell the algorithm (with ``PPG_ABP_idx property``)
	  which signal indexes to process, from 1 to ``ECG_header.nsig``

``PPG_ABP_idx`` — The indexes corresponding to pulsatile signals ``[] (default)`` 

	A value from 1 to ``ECG_header.nsig`` indicating the column indexes (of the
	signal matrix ``ECG``) where pulsatile signals are located. By default this task
	process all signals.
	

``detectors`` — The PPG/ABP detection algorithms to use ``'all-detectors' (default)`` 

	Select which algorithm to use. A cell string or char with any of the following names:

	- 'all-detectors'

	- `'wavePPG' <http://dx.doi.org/10.1109/JBHI.2013.2267096>`__

	- `'wabp' <http://www.physionet.org/physiotools/wag/wabp-1.htm>`__


Examples
--------

Create the *ECGtask\_PPG\_ABP\_detector* object.

.. code::

    % with the task name
        ECG_w.ECGtaskHandle = 'PPG_ABP_detector';
    % or create an specific handle to have more control
        ECGt_PPG = ECGtask_PPG_ABP_detector();

and then you are ready to set the algorithms to use. In the following
example you have several possible set-ups.

.. code::

    % select an specific algorithm. Default: Run all detectors
        ECGt_PPG.detectors = 'wavePPG'; % A J. Lazaro algorithm for peak detection
        ECGt_PPG.detectors = 'wabp';  % Another algorithm from Physionet
                            

Finally set the task to the wrapper object, and execute the task.

.. code::

            ECG_w.ECGtaskHandle= ECGt_PPG; % set the ECG task
            ECG_w.Run();

You can check the result of this task, with either the :doc:`detection
corrector <ABP_PPG_peak_correction>` or the :doc:`visualization
functions <plot_ecg_strip>`.

Also check this :ref:`example <PPG_ABP_pulse_detection>` for
further information.

.. _pulse_det_result_format:

Results format
--------------
 
The result file have the same format than :ref:`QRS detection task <QRS_det_result_format>`.


More About
----------

Here are some external references about pulse detection:

-  **?? Add some**

See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`QRS detection <QRS_detection>` \| :doc:`examples <examples>`
	
