
Accessing results
=================

Results are stored in a ``mat`` file for compatibility reasons. The format depends on the 
task that generated the results, but a typical procedure to grab data from experiments 
is:

.. code::

    ECGw = ECGwrapper('recording_name', 'your_rec_filename');
    result_filename = ECGw.GetCahchedFileName('QRS_detection');
    results = load(cached_filenames{1});


In this example, the results from the previous QRS detection experiment is loaded in the ``results`` 
variable. The format for the specific tasks was described in the following links:

- :ref:`QRS detection <QRS_det_result_format>` 

- ABP/PPG pulse detection tasks have the same format of :ref:`QRS detection <QRS_det_result_format>`

- :ref:`ECG delineation <Delineation_result_format>`

- :ref:`Heartbeat classifier <Classifier_det_result_format>`

- :ref:`Arbitrary tasks <arbitrary_result_format>`

