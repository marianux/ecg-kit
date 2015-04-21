
reportECG
=========

Description
-----------

This function creates a report of a signal handled by an ECGwrapper object.
The report includes several views of the signals at different time scales. 
In addition, you have the possibility to overprint information from other 
ECGtask results, such as QRS detections, wave delineation, and heartbeat types. 
Some aspects of the report can be configured as the detail degree, the length 
of each time scale and the report format.


.. toctree::
   :titlesonly:
   :hidden:

   A signal visualization tool <plot_ecg_strip>
   A mosaic visualization tool <plot_ecg_mosaic>

   
Syntax
------

The function prototype is


.. code::

    function reportECG(ECG_w, detailLevel, report_mode, win_lengths, report_format, filename)
	

where the arguments are:

 - ``ECG_w`` An ECGwrapper object as the signal handler.
 
 - ``detailLevel`` The report detail level:
   
   - 'HighDetail'
   
   - 'MediumDetail'
   
   - 'LowDetail'
   
   
   A higher detail level means report the whole recording at every time resolution defined 
   in "win_lengths". High resolution also means larger reports. ``LowDetail (default)``.
 
 - ``report_mode`` Information from other tasks like QRS detection/delineation/classification 
   added to the signals in case available mode. Possible values are:
               
   'full' 
   
   ``ECG only (default)`` 
   
   'QRS detection' 
   
   'Wave delineation'
               
   'Heartbeat classification'
 
 - ``win_lengths`` The amount and size (in seconds) of each scale length present in the report. 
   ``[60*60 30*60 60 7] (default)``. It means 1 hour - 30 min - 1 min and 7 seconds.
 
 - ``report_format`` The report format of the document. ``PDF (default)``.
 
 - ``filename`` The report filename. ``rec_folder\rec_name.report_format (default)``.


Examples
--------

The example folder has some examples of the use of the reporting functions.

.. code::

    reportECG(ECGw, 'LowDetail', 'full');

This is an example of an ECG overview
	
.. image:: 208_full_03.png

And this with more information overprint

.. image:: 208_full_14.png
	
	
See Also
--------

 :doc:`Plot ECG strip <plot_ecg_strip>` \| :doc:`Plot ECG mosaic <plot_ecg_mosaic>`
