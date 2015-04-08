
ABP/PPG peak detection task
===========================

Allow the access to ECG recordings of arbitrary format and length.

`expand all in page `__

Syntax
------

-  ``ECGw = ECGwrapper()`` `example <ECGwrapper.html#ecgw_ex_noarg>`__
-  ``ECGw = ECGwrapper(Name,Value)``
   `example <ECGwrapper.html#ecgw_ex_namevalue>`__

 

Description
-----------

ECG wrapper is a class to allow the access to cardiovascular signal
recordings of several formats (MIT, ISHNE, AHA, HES, MAT) and lengths,
from minutes to days. Also an ECGtask object can be plugged to it to
perform several tasks, such as QRS detection and ECG delineation among
others.

`example <ECGwrapper.html#ecgw_ex_noarg>`__

``ECGwrapper()`` creates an ECGwrapper object with its default values.
You can set later through its properties the recording filename to work
with, and the task to perform on it.

`example <ECGwrapper.html#ecgw_ex_namevalue>`__

``ECGwrapper(Name,Value``) creates an ECGwrapper object with options
specified by one or more name-value pair arguments.

 

Input Arguments
---------------

Name-Value Pair Arguments
~~~~~~~~~~~~~~~~~~~~~~~~~

Specify optional comma-separated pairs of ``Name,Value`` arguments.
``Name`` is the argument name and ``Value`` is the corresponding value.
``Name`` must appear inside single quotes (``' '``). You can specify
several name and value pair arguments in any order as
``Name1,Value1,...,NameN,ValueN``.

**Example:**
``'recording_name','/ecg_recordings/rec1.dat',                                       'ECGtaskHandle', 'QRS_detection' ``\ specifies
to create an ECGwrapper to detect heartbeats in recording
'/ecg\_recordings/rec1.dat' LaTeX output file format and excludes the
code from the output.

```'recording_name'`` — ECG Recording full name `__\ ``''`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The full path filename of the ECG recording.

```'recording_format'`` — ECG recording format `__\ ``''`` (default auto-detect) \| ``'MIT'`` \| ``'ISHNE'`` \| ``'AHA'`` \| ``'HES'`` \| ``'MAT'``\ \| ``'Mortara'``
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The format of the ECG recording. By default or if not specified, the
wrapper will attemp to auto-detect the format.

+--------------------------------------------------------------------+----------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------+----------------------------------------------------------------------------+---------------------------------------------------------------------------------+
| Output Format                                                      |
| String Value                                                       |
+====================================================================+======================================================================+===================================================================================================================================================================================================================+=========================================================================+============================================================================+=================================================================================+
| 'MIT'                                                              | 'ISHNE'                                                              | 'AHA'                                                                                                                                                                                                             | 'HES'                                                                   | 'MAT'                                                                      | 'Mortara'                                                                       |
| ``MIT                                                   format``   | ``ISHNE                                                   format``   | ``American                                                   Heart Association ECG                                                   Database or Physionet                                                   ``   | ``Biosigna                                                   format``   | ``Matlab                                                   file format``   | ``Mortara                                                   SuperECG format``   |
+--------------------------------------------------------------------+----------------------------------------------------------------------+-------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------------+-------------------------------------------------------------------------+----------------------------------------------------------------------------+---------------------------------------------------------------------------------+

```'this_pid'`` — Process identification `__\ ``'1/1'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case working in a multiprocess environment, this value will identify
the current process. Can be a numeric value, or a string of the form
'N/M'. This pid is N and the total amount of pid's to divide the whole
work is M.

```'tmp_path'`` — The path to store temporary data `__\ ``'tempdir()'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges.

```'output_path'`` — The output path to store results `__\ ``'fileparts(recording_name)'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges. By default will be the
same path of the recordings.

```'ECGtaskHandle'`` — The task to perform. `__\ ``''`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The task to perform, can be the name of the task, or an ECGtask object.
Available ECGtasks can be listed with
`list\_all\_ECGtask() <matlab:doc('list_all_ECGtask')>`__ command.

```'partition_mode'`` — The way that this object will partition lengthy signals `__\ ``'ECG_overlapped'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The way that this object will partition lengthy signals:

-  'ECG\_contiguous' no overlapp between segments.

-  'ECG\_overlapped' overlapp 'overlapping\_time' between segments.

-  'QRS' do the partition based on the annotations in
   ECG\_annotations.time property. Typically but not necessary are QRS
   annotations.

```'overlapping_time'`` — Time in seconds of overlapp among consequtive segments `__\ ``30`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time in seconds of overlapp among consequtive segments. This segment is
useful for ensuring transitory responses of systems to be finished.

```'cacheResults'`` — Save intermediate results to recover in case of failure `__\ ``true`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Save intermediate results to recover in case of failure.

```'syncSlavesWithMaster'`` — Time in seconds of overlapp among consequtive segments `__\ ``false`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In multiprocess environments sometimes it is useful to terminate all
pid's together in order to start subsequent tasks synchronously.

```'repetitions'`` — Times to repeat the ECGtask `__\ ``1`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case the ECGtask is not deterministic, the repetition property allows
to repeat the task several times.

 

Examples
--------

`collapse all `__

`Create the simplest ECG wrapper object `__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create the ECGwrapper object.

.. code::

    >> ECG_w = ECGwrapper()
    ECG_w = 
    ############################
    # ECGwrapper object config #
    ############################
    +ECG recording: None selected
    +PID: 1/1
    +Repetitions: 1
    +Partition mode: ECG_overlapped
    +Function name: Null task
    +Processed: false
                    

Then, in your script or in the command window you can type:

.. code::

    >> ECG_w.recording_name = 'some_path\100';
    >> ECG_w.ECGtaskHandle = 'QRS_detection'
    ECG_w = 
    ############################
    # ECGwrapper object config #
    ############################
    +ECG recording: some_path\100 (auto)
    +PID: 1/1
    +Repetitions: 1
    +Partition mode: ECG_overlapped
    +Function name: QRS_detection
    +Processed: false
                    

Now, you just want to run the task by executing:

.. code::

    >> ECG_w.Run();
                    

`Create an ECGwrapper object for an specific recording and task `__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, we create the same object of the previous example but
using the name-value .

.. code::

    >> ECG_w = ECGwrapper( ...
            'recording_name', 'some_path\100', ...
            'recording_format', 'MIT', ...
            'ECGtaskHandle', 'QRS_detection', ...
            )
    ECG_w = 
    ############################
    # ECGwrapper object config #
    ############################
    +ECG recording: some_path\100 (auto)
    +PID: 1/1
    +Repetitions: 1
    +Partition mode: ECG_overlapped
    +Function name: QRS_detection
    +Processed: false
                        
    >> ECG_w.Run();
                    

 

More About
----------

`expand all `__

Other resources
---------------

-  `Physionet.org <http://physionet.org/>`__
-  `Telemetric and Holter ECG Warehouse
   (THEW) <http://thew-project.org/>`__
-  `Pablo Laguna research group at University of
   Zaragoza <http://diec.unizar.es/~laguna/personal/publicaciones/publicaciones.htm>`__
-  `Computing in Cardiology <http://cinc.org/>`__

See Also
--------

```ECGtask`` <ECGtask.html>`__ \| ```examples`` <examples.html>`__

