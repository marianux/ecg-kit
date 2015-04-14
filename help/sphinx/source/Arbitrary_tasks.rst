
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

```only_ECG_leads`` — Process only ECG signals `__\ ``true`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boolean value. Find out which signals are ECG based on their header
description.

```payload`` — An arbitrary format variable to be passed to your user-defined algorithm. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This variable can be useful for passing data to your own function, not
covered in the interface described
`below <#Adding_a_custom_detection_algorithm>`__.

```signal_payload`` — Treat the result of your arbitrary function as a signal. `__\ ``false`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boolean value that indicates the ECGwrapper to produce a signal or
result payload.

```lead_idx`` — The signal indexes that your function will affect. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A positive integer array with values from 1 to ECG\_header.nsig.

```function_pointer`` — The pointer to your arbitrary function. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Your function must follow this prototype:

.. code::
    function result = your_function( ECG_matrix, ECG_header, progress_handle, payload_in)  
                            

where the arguments are:

ECG\_matrix, is a matrix size [nsamp length(lead\_idx)], where nsamp is
handled internally by the ECG wrapper and `lead\_idx <#lead_idx_prop>`__
property is user-defined.

ECG\_header, is a struct with info about the ECG signal, such as:

freq, the sampling frequency

desc, description about the signals.

and others described `here <Copy_of_ECGtask.htm>`__

progress\_handle, is a handle to a `progress\_bar <progress_bar.htm>`__
object, that can be used to track the progress within your function.

payload\_in, is a user variable, of arbitrary format, allowed to be sent
to your function. It is sent, via the `payload
property <#payload_prop>`__ of this class, for example:

.. code::
        % One variable
        this_ECG_wrapper.ECGtaskHandle.payload = your_variable;
        
        % Several variables with a cell container
        this_ECG_wrapper.ECGtaskHandle.payload = {your_var1 your_var2}; 

and the output of your function must be a result (struct) variable, or
can be handled as a signal with signal\_payload property.

Examples
--------

This example is used in the QRScorrector function to perform
template-matching on an ECGwrapper (arbitrary big recording) object.

.. code::
    aux_w = ECGwrapper('recording_name', 'your_path/recname');
    aux_w.ECGtaskHandle = 'arbitrary_function';
    % This is in case you want always to recalculate results, no cachingaux_w.cacheResults = false;
    % Use first and third columns-signalsaux_w.ECGtaskHandle.lead_idx = [1 3];
    % Produce a signal as a resultaux_w.ECGtaskHandle.signal_payload = true;
    % Add a user-string to identify the runaux_w.ECGtaskHandle.user_string = ['similarity_calc_for_lead_' num2str(sort(lead_idx)) ];
    % add your function pointeraux_w.ECGtaskHandle.function_pointer = @similarity_calculation;
    % and any data your function may need.aux_w.ECGtaskHandle.payload = pattern2detect;% and you are ready to go !
    aux_w.Run
                            

See Also
--------

```ECGtask`` <ECGtask.html>`__ \| ``ECGwrapper`` \|
```examples`` <examples.html>`__

