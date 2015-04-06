 

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

ABP/PPG peak detection
======================

This document describes how to perform automatic peak detection in
pulsatile signals.

`expand all in page <javascript:void(0);>`__

 

Description
===========

This task perform peak detection in pulsatile signals such as arterial
blood pressure (ABP) or plethysmographic (PPG). The task uses two
algorithms to achieve pulse detection, and, as in QRS detection task you
can choose to use any of them.

 

Input Arguments
===============

The properties that the ECGtask\_PPG\_ABP\_detector class accepts are
described below. The usage of these properties is restricted to
low-level programming, you can use this task through the ECGwrapper as
is shown in the example below.

```progress_handle`` — used to track the progress within your function. <javascript:void(0);>`__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

progress\_handle, is a handle to a `progress\_bar <progress_bar.htm>`__
object, that can be used to track the progress within your function.

```tmp_path`` — The path to store temporary data <javascript:void(0);>`__\ ``'tempdir()'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

Full path to a directory with write privileges.

```lead_config`` — Select which signals to process <javascript:void(0);>`__\ ``'PPG-ABP-only'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A cell string or char with any of the following names

'all-leads'. Process all leads

'PPG-ABP-only'. Detect pulsatile signals based on their ECG\_header.desc
description variable.

'User-defined-leads'. Tell the algorithm (with PPG\_ABP\_idx property)
which signal indexes to process, from 1 to ECG\_header.nsig

```PPG_ABP_idx`` — The indexes corresponding to pulsatile signals <javascript:void(0);>`__\ ``[]`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A value from 1 to ECG\_header.nsig indicating the column indexes (of the
signal matrix of nsamp x nsig) where pulsatile signals are located.

```detectors`` — The PPG/ABP detection algorithms to use <javascript:void(0);>`__\ ``'all-detectors'`` (default)
~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~~

A cell string or char with any of the following names

all-detectors

`wavePPG <http://dx.doi.org/10.1109/JBHI.2013.2267096>`__

`wabp <http://www.physionet.org/physiotools/wag/wabp-1.htm>`__

 

Examples
========

Create the ECGtask\_PPG\_ABP\_detector object.

.. code:: codeinput

    % with the task name
        ECG_w.ECGtaskHandle = 'PPG_ABP_detector';
    % or create an specific handle to have more control
        ECGt_PPG = ECGtask_PPG_ABP_detector();

and then you are ready to set the algorithms to use. In the following
example you have several possible setups.

.. code:: codeinput

    % select an specific algorithm. Default: Run all detectors
        ECGt_PPG.detectors = 'wavePPG'; % A J. Lazaro algorithm for peak detection
        ECGt_PPG.detectors = 'wabp';  % Another algorithm from Physionet
                            

Finally set the task to the wrapper object, and execute the task.

.. code:: codeinput

            ECG_w.ECGtaskHandle= ECGt_PPG; % set the ECG task
            ECG_w.Run();

You can check the result of this task, with either the `detection
corrector <ABP_PPG_peak_correction.htm>`__ or the `visualization
functions <plot_ecg_strip.htm>`__.

Also check this `example <examples.html#PPG_ABP_pulse_detection>`__ for
further information.

 

More About
==========

Here are some external references about pulse detection:

-  `Physionet Challenge 2014 <http://physionet.org/challenge/2014/>`__

See Also
========

```ECGtask`` <ECGtask.html>`__ \| ``QRS               detection`` \|
```ECG delineation`` <ECGdelineation.htm>`__ \|
```examples`` <examples.html>`__

 

.. |image0| image:: template/my_layout/Search.png
   :target: #
.. |image1| image:: template/my_layout/Print.png
   :target: javascript:window.print()
.. |image2| image:: template/my_layout/Search.png
   :target: #
.. |image3| image:: template/my_layout/Print.png
   :target: javascript:window.print()
