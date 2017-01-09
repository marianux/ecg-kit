
Matlab file format
==================

The ecg-kit allows users to use Matlab native format to store and access recordings. For using it,
ensure that your mat file follows these naming convention:
 
 
 1. Your signal variable must be named ``'sig', 'signal' or 'ECG'``. Remember that the signals or leads 
    must be placed column-wise, that is ``signal = [ lead1 lead2 lead3 ... lead_nsig]``.
 2. The header or signal information must be stored in a struct named ``'header', 'heasig' or 'hea'``. The header variable is a struct with the fields defined :ref:`here <header_format>`.
 3. (Optional) In case of including annotations or QRS detections, be sure to be a struct named ``'ann', 
    'annotations' or 'qrs'`` and which includes the fields described for the MIT format in `Physionet <http://www.physionet.org/physiobank/annotations.shtml>`__.

    * time: the time within the recording (recorded in the annotation file as the sample number of the sample to which the annotation "points")
    * anntyp [sic]: a numeric annotation code (see ecgcodes.h for definitions)
    * subtyp [sic], chan, num: three small integers (between -128 to 127) that specify context-dependent attributes (see the documentation for each database for details)
    * aux: a free text string


See the ``ecg-kit\common\matformat_definitions.m`` for more details.
