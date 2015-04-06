 

+--------------------------------------+--------------------------------------+
| |image2|                             |
| |image3|                             |
| `Show <javascript:onShowHideClick()> |
| `__\ `Hide <javascript:onShowHideCli |
| ck()>`__                             |
+--------------------------------------+--------------------------------------+

-  Contents
-  Index
-  Glossary

| 

plot\_ecg\_mosaic
=================

Low level function to produce mosaic charts of signals and annotations.

`expand all in page <javascript:void(0);>`__

Syntax
======

-  ``ECGw = ECGwrapper()`` `example <ECGwrapper.html#ecgw_ex_noarg>`__
-  ``ECGw = ECGwrapper(Name,Value)``
   `example <ECGwrapper.html#ecgw_ex_namevalue>`__

 

Description
===========

This is the main class of the toolbox since it allows access to
cardiovascular signal recordings of several formats
(`MIT <http://www.physionet.org/physiotools/wag/signal-5.htm>`__,
`ISHNE <http://thew-project.org/THEWFileFormat.htm>`__,
`AHA <https://www.ecri.org/Products/Pages/AHA_ECG_DVD.aspx>`__, HES,
`MAT <Matlab_format.htm>`__) and lengths, from seconds to days. The
objective of this class is to provide `ECGtask <ECGtask.htm>`__  class,
a common interface to access data and perform specific tasks. Briefly,
this class sequentially reads data and passes to the Process method of
the `ECGtask <ECGtask.htm>`__ plugged in the
`ECGtaskHandle <#inputarg_ECGtask>`__ property. Some common tasks, such
as `QRS detection <examples.html#QRS_automatic_detection>`__ and `ECG
delineation <examples.html#ECG_automatic_delineation>`__, can be easily
invoked. Also other predefined tasks or your own code can be adapted as
is shown in the `examples <examples.html>`__.

A more detailed description of this class, together with an explanation
of how you can easily hook your algorithms to this class is
`here <extensions.htm>`__.

Finally the results produced by the `ECGtask <ECGtask.htm>`__ are stored
in order to ease reproducibility and backup of your experiments, or to
be used of subsequent tasks as shown in the
`examples <examples.html>`__.

 

Input Arguments
===============

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
'/ecg\_recordings/rec1.dat'.

```'recording_name'`` — ECG Recording full name <javascript:void(0);>`__\ ``''`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The full path filename of the ECG recording.

```'recording_format'`` — ECG recording format <javascript:void(0);>`__\ ``''`` (default auto-detect) \| ``'MIT'`` \| ``'ISHNE'`` \| ``'AHA'`` \| ``'HES'`` \| ``'MAT'``\ \| ``'Mortara'``
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

```'this_pid'`` — Process identification <javascript:void(0);>`__\ ``'1/1'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case working in a multiprocess environment, this value will identify
the current process. Can be a numeric value, or a string of the form
'N/M'. This pid is N and the total amount of pid's to divide the whole
work is M.

```'tmp_path'`` — The path to store temporary data <javascript:void(0);>`__\ ``'tempdir()'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges.

```'output_path'`` — The output path to store results <javascript:void(0);>`__\ ``'fileparts(recording_name)'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges. By default will be the
same path of the recordings.

```'ECGtaskHandle'`` — The task to perform. <javascript:void(0);>`__\ ``''`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The task to perform, can be the name of the task, or an ECGtask object.
Available ECGtasks can be listed with
`list\_all\_ECGtask() <matlab:doc('list_all_ECGtask')>`__ command.

````

```'partition_mode'`` — The way that this object will partition lengthy signals <javascript:void(0);>`__\ ``'ECG_overlapped'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

The way to do batch partition in lengthy signals:

-  'ECG\_contiguous' no overlapp between segments.

-  'ECG\_overlapped' overlapp of 'overlapping\_time' among segments.
   This can be useful if your task have a transient period to avoid.

-  'QRS' do the partition based on the annotations provided in
   ECG\_annotations.time property. This option is useful if your task
   works in the boundaries of a fiducial point (commonly a heartbeat),
   and not in the whole signal. This partition mode ignores those parts
   of the recording without annotations.

```'overlapping_time'`` — Time in seconds of overlapp among consequtive segments <javascript:void(0);>`__\ ``30`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Time in seconds of overlapp among consequtive segments. This segment is
useful for ensuring the end of all transients within a task.

```'cacheResults'`` — Save intermediate results to recover in case of failure <javascript:void(0);>`__\ ``true`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Save intermediate results to recover in case of errors. Useful for long
jobs or recordings.

```'syncSlavesWithMaster'`` — Time in seconds of overlapp among consequtive segments <javascript:void(0);>`__\ ``false`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In multiprocess environments sometimes it is useful to terminate all
pid's together in order to start subsequent tasks synchronously. This
value forces all parts of a multipart process to wait until all other
parts finish.

```'repetitions'`` — Times to repeat the ECGtask <javascript:void(0);>`__\ ``1`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In case the ECGtask is not deterministic, the repetition property allows
to repeat the task several times.

 

Examples
========

`collapse all <javascript:void(0);>`__

`Create the simplest ECG wrapper object <javascript:void(0);>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Create the ECGwrapper object.

.. code:: programlisting

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

.. code:: programlisting

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

.. code:: programlisting

    >> ECG_w.Run();
                    

`Create an ECGwrapper object for an specific recording and task <javascript:void(0);>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

In this case, we create the same object of the previous example but
using the name-value .

.. code:: programlisting

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
==========

`expand all <javascript:void(0);>`__

 

`Other resources <javascript:void(0);>`__
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

-  `Physionet.org <http://physionet.org/>`__
-  `Telemetric and Holter ECG Warehouse
   (THEW) <http://thew-project.org/>`__
-  `Pablo Laguna research group at University of
   Zaragoza <http://diec.unizar.es/~laguna/personal/publicaciones/publicaciones.htm>`__
-  `Computing in Cardiology <http://cinc.org/>`__

See Also
========

```ECGtask`` <ECGtask.html>`__ \| ```examples`` <examples.html>`__

 

.. |image0| image:: template/my_layout/Search.png
   :target: #
.. |image1| image:: template/my_layout/Print.png
   :target: javascript:window.print()
.. |image2| image:: template/my_layout/Search.png
   :target: #
.. |image3| image:: template/my_layout/Print.png
   :target: javascript:window.print()
