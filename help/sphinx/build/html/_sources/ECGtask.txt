
ECGtask
=======

Perform an specific task to an ECG signal. This document defines the standard interface that tasks must implement.

.. toctree::
   :maxdepth: 4
   :titlesonly:
   :hidden:
   
   QRS detection <QRS_detection>
   QRS correction <QRScorrector>
   ABP/PPG peak detection <ABP_PPG_peak_detection>
   ABP/PPG peak correction <ABP_PPG_peak_correction>
   ECG delineation <ECGdelineation>
   ECG delineation correction <ECG_delineation_corrector>
   Heartbeat classification <ECG_heartbeat_classifier>
   Arbitrary tasks <Arbitrary_tasks>



Description
-----------

The ECGtask is an abstract class definition where the minimum interface
requirements are specified, in order that your own tasks can be safely
plugged into `ECGwrapper <ECGwrapper>` objects. As an example of
how to use this interface, see the derived classes for :doc:`QRS 
detection <QRS_detection>` and :doc:`ECG delineation <ECGdelineation>`, 
among others that can be listed
with the :doc:`list\_all\_ECGtask <list_all_ECGtask>` function:

-  :doc:`QRS detection <QRS_detection>`         
-  :doc:`QRS correction <QRScorrector>`         
-  :doc:`ECG delineation <ECGdelineation>`         
-  :doc:`ECG delineation correction <ECG_delineation_corrector>`         
-  :doc:`ABP/PPG peak detection <ABP_PPG_peak_detection>`         
-  :doc:`ABP/PPG peak correction <ABP_PPG_peak_correction>`         
-  :doc:`Heartbeat classification <ECG_heartbeat_classifier>`         
-  :doc:`Arbitrary tasks <Arbitrary_tasks>`         

These tasks are the core of this kit and you will probably refer to them
before you extend the functionality with your own tasks.

Properties
----------

All tasks must implement the following properties with its attributes:

.. code::

	properties(GetAccess = public, Constant)


``name`` — The name of the task.

``target_units`` — The signal units required by the task. Possible values are: 
 
 ADCu, raw ADC samples.
 
 nV, uV, mV, V, voltage.

``doPayload`` — Boolean. Does this task generates a payload to be stored ? 


.. code::

	properties(GetAccess = public, SetAccess = private)


``memory_constant`` — A coefficient to indicate the ECGwrapper how big should 
be a batch processing part. The size of each part is calculated as 

.. code::

	user = memory;
	batch_size = memory_constant * user.MaxPossibleArrayBytes;

	
``started`` — Boolean. Did the task executed the Start method ?


.. code::

	properties(GetAccess = public, SetAccess = public)

	
	
``progress_handle`` — is a handle to a :doc:`progress\_bar <progress_bar>`
object, that can be used to track the progress within your function.

``tmp_path`` — The path to store temporary data.


Methods
------- 

All tasks must implement the following methods:


``Start`` — The task initialization method.

This task initialize specific aspects of the task.

.. code::

	Start(obj, ECG_header, ECG_annotations)
        
							
where the arguments are:

	**ECG\_header**, is a struct with info about the ECG signal, See :ref:`here <ECG_header_description>` 
	for a description.

	**ECG_annotations**, Commonly QRS detections, signal quality annotations or other type of measurements
	included with the recordings. Some documentation about annotations in `Physionet <http://www.physionet.org/physiobank/annotations.shtml>`__.

	
``Process`` — The task core processing function.

This task is the responsible of do the actual work of the ECGtask. This mehtod is called by an ECGwrapper
all the times needed to process the whole recording.
	
.. code::

        payload = Process(ECG, 
			  ECG_start_offset, 
			  ECG_sample_start_end_idx, 
			  ECG_header, 
			  ECG_annotations, 
			  ECG_annotations_start_end_idx )

		
where the arguments are:


	**ECG**, is a matrix size ``[ECG_header.nsamp ECG_header.nsig]``
		
	**ECG\_start\_offset**, is the location of ECG(1,:) within the whole signal.
		
	**ECG\_header**, is a struct with info about the ECG signal, See :ref:`here <ECG_header_description>` 
	for a description.
		
	**ECG_annotations**, Commonly QRS detections, signal quality annotations or other type of measurements
	included with the recordings. Some documentation about annotations in `Physionet <http://www.physionet.org/physiobank/annotations.shtml>`__.
		
	**ECG\_annotations\_start\_end\_idx**, are the start and end indexes corresponding
	to the first and last element of ECG_annotations in the current iteration.

as a result, this method must produce a **payload** variable, that will be handled by the ECGwrapper object.
		
``Concatenate`` — This method is responsible of the payload union after all the processing.

After the execution of all *Process* steps, each payload must be put together with this method. The ECGwrapper
object will call this method once for each payload created, building a final payload.
	
.. code::

        payload = Concatenate(plA, plB)

		
where the arguments are:

	**plA** and **plB**, are two payloads created with the *Process* method.

and as a result, this method creates **payload**, the union of **plA** and **plB**.	

	
``Finish`` — This task perform the last calculation over the whole payload.

After the concatenation of payloads, the whole payload is sent to this method to 
perform any final calculation.
	
.. code::

        payload = Finish(obj, payload, ECG_header)
		
where the arguments are:

	**payload**, is the payload created with all the *Concatenate* method invocation.
				
	**ECG\_header**, is a struct with info about the ECG signal, See :ref:`here <ECG_header_description>` 
	for a description.
				
As a result, the final payload is generated, which later will be stored by the ECGwrapper object.



More About
----------

-  `Physionet.org <http://physionet.org/>`__
-  `Telemetric and Holter ECG Warehouse
   (THEW) <http://thew-project.org/>`__
-  `Pablo Laguna research group at University of
   Zaragoza <http://diec.unizar.es/~laguna/personal/publicaciones/publicaciones.htm>`__
-  `Computing in Cardiology <http://cinc.org/>`__

 

See Also
--------

:doc:`ECGwrapper <ECGwrapper>` \|
:doc:`ECG\_delineation <ECGdelineation>` \|
:doc:`list\_all\_ECGtask <list_all_ECGtask>`

