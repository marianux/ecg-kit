
ECGtask
-------

Perform an specific task to an ECG signal.

`expand all in page `__

 

Description
-----------

The ECGtask is an abstract class definition where the minimum interface
requirements are specified, in order that your own tasks can be safely
plugged into `ECGwrapper <ECGwrapper.html>`__ objects. As an example of
how to use this interface, see the derived classes for `QRS
detection <QRS_detection.htm>`__ and `ECG
delineation <ECGdelineation.htm>`__, among others that can be listed
with the `list\_all\_ECGtask <list_all_ECGtask.htm>`__ function:

-  `ECG\_delineation <ECGdelineation.htm>`__            
-  `ECG\_delineation\_corrector <ECG_delineation_corrector.htm>`__  
-  `PPG\_ABP\_corrector <ABP_PPG_peak_correction.htm>`__          
-  `PPG\_ABP\_detector <ABP_PPG_peak_detection.htm>`__           
-  `QRS\_corrector <QRScorrector.htm>`__              
-  `QRS\_detection <QRS_detection.htm>`__              
-  `ECG\_heartbeat\_classifier <ECG_heartbeat_classifier.htm>`__   
-  `arbitrary\_function <Arbitrary_tasks.htm>`__

These tasks are the core of this kit and you will probably refer to them
before you extend the functionality with your own tasks.

 

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

`ECGwrapper <ECGwrapper.html>`__ \|
`ECG\_delineation <ECGdelineation.htm>`__ \|
`list\_all\_ECGtask <list_all_ECGtask.htm>`__

