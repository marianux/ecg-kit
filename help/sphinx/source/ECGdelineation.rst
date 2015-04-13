
ECG delineation
===============

This document describes how to perform automatic delineation or wave
segmentation on ECG signals.

`expand all in page `__

 

Description
-----------

Automatic wave segmentation or delineation is exclusively performed by
wavedet algorithm.

 

Input Arguments
---------------

The properties that this task uses are the following:

```progress_handle`` — used to track the progress within your function. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

progress\_handle, is a handle to a `progress\_bar <progress_bar.htm>`__
object, that can be used to track the progress within your function.

```tmp_path`` — The path to store temporary data `__\ ``'tempdir()'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges.

```delineators`` — The ECG delineation algorithms to use `__\ ``'all-delineators'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A cell string or char with any of the following names

all-delineators

`wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__

 

```only_ECG_leads`` — Process only ECG signals `__\ ``true`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Boolean value. Find out which signals are ECG based on their header
description.

```wavedet_config`` — A structure for customizing wavedet algorithm. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Undocumented yet, use only if you know what you are doing.

```payload`` — An arbitrary format variable. `__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

This variable can be useful for passing data to your own function
(described `below <#Adding_a_custom_detection_algorithm>`__) or to
provide delineation algorithms visually audited QRS detections.

 

Adding a custom delineation algorithm
-------------------------------------

Adding your own delineator to the kit is very simple. Ensure that your
function implements this interface:

.. code::
    function [positions_single_lead, position_multilead] = your_ECG_delineation( ECG_matrix, ECG_header, progress_handle, payload_in)  
                            

where the arguments are:

ECG\_matrix, is a matrix size [ECG\_header.nsamp ECG\_header.nsig]

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
        
        % Or the result of a previous task, in this case QRS manual correction (if available)
        % or the automatic detection if not.
        cached_filenames = this_ECG_wrapper.GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
        this_ECG_wrapper.ECGtaskHandle.payload = load(cached_filenames);

and the output of your function must be:

positions\_single\_lead, a cell array size ECG\_header.nsig with the QRS
sample locations found in each lead.

position\_multilead, a numeric vector with the QRS locations calculated
using multilead rules.

Examples
--------

Create the ECGtask\_ECG\_delineation object.

.. code::
    % with the task name
        ECG_w.ECGtaskHandle = 'ECG_delineation';
    % or create an specific handle to have more control
        ECGt = ECGtask_ECG_delineation();

and then you are ready to set the algorithms to use. In the following
example you have several possible setups.

.. code::
    % select an specific algorithm. Default: Run all detectors
            ECGt.delineators = 'wavedet'; % Wavedet algorithm based on
            ECGt.delineators = 'user:your_delineator_func_name';    % "your_delineator_func_name" can be your own delineator.
            ECGt.delineators = {'wavedet' 'user:your_delineator_func_name'};
                            

Finally set the task to the wrapper object, and execute the task.

.. code::
            ECG_w.ECGtaskHandle= ECGt; % set the ECG task
            ECG_w.Run();

You can check the result of this task, with either the `delineation
corrector <ECG_delineation_corrector.htm>`__ or the `visualization
functions <plot_ecg_strip.htm>`__.

Also check this `example <examples.html#ECG_automatic_delineation>`__
for further information.

 

More About
----------

This publication describes the
`wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__
algorithm:

See Also
--------

```ECGtask`` <ECGtask.html>`__ \| ``QRS                   detection`` \|
```examples`` <examples.html>`__

