
QRS correction
==============

This document describes how to perform inspection and correction of
automatic heartbeat detection.

Description
-----------

Automatic heartbeat detection is commonly well performed in those
recordings with stable heart rhythms and QRS morphologies. In those
cases where these situations are not met, many problems arise and
automatic detection is not easy performed. This task provides a
graphical user interface (GUI) to ease verification, correction and even
manual detection.

 

Input Arguments
---------------

The properties that the ECGtask\_QRS\_corrector class handle are
described below. The usage of these properties is restricted to
low-level programming, you can use this task through the ECGwrapper as
is shown in the example below.

``progress_handle`` — Used to track the progress within your function. ``[] (default)``

	progress\_handle, is a handle to a :doc:`progress\_bar <progress_bar>`
	object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data. ``tempdir() (default)``

	Full path to a directory with write privileges.

.. _payload_prop:

``payload`` — An arbitrary format variable to be passed to your user-defined algorithm. ``[] (default)``

	This propery is typically used to pass the automatic detection results. See the example below.

``caller_variable`` — An arbitrary variable name in the caller workspace. ``'payload' (default)``

	This variable in the caller workspace, will be assigned after user-interaction ends. ECGtask\_QRS\_corrector uses 'payload' as default variable in
	order to save the result of edition/verification with the GUI.

Examples
--------

Create the ECGtask\_QRS\_corrector object.

.. code::

    ECGw.ECGtaskHandle = 'QRS_corrector';
    % this is to use previous saved results as starting point, if any available
    cached_filenames = ECG_all_wrappers(ii).GetCahchedFileName({'QRS_corrector' 'QRS_detection'});
    ECGw.ECGtaskHandle.payload = load(cached_filenames{1});
    ECGw.Run();

Then the following GUI appears

|image4|

and the command window shows the following message:

.. code::

    #############################
    # User interaction required #
    #############################
    This ECGtask allow user interaction. Press [CTRL + G] in figure 1 to save results and press F5 (Run) to continue.
    K>>

see the videos in
`YouTube <https://www.youtube.com/watch?v=qgWjvsvafVg&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=3>`__
for a more detailed demo about things you can do with the GUI. The demo shows how to
add/remove heartbeats annotation, browse the detections through the whole recording and
cut/copy/paste heartbeat detections between different annotations. Other features not 
described in this video were added in `this other <https://www.youtube.com/watch?v=qgWjvsvafVg&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=3>`__. 
After edition/verification of the automatic delineation, press CTRL+G to save
results in the 'payload' variable of the caller workspace. Then press F5
to save the results to disk.

Results format
--------------
 
The result file have the same format than :ref:`QRS detection task <QRS_det_result_format>`.


More About
----------

Here are some external references about heartbeat detection:

-  A video demo in `Youtube <https://www.youtube.com/watch?v=qgWjvsvafVg&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=3>`__


See Also
--------

 :doc:`ECGtask <ECGtask>` \| :doc:`ECG delineation <ECGdelineation>` \| :doc:`examples <examples>`

.. |image4| image:: QRS_corrector.PNG
