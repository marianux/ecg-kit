.. ECGkit documentation master file, created by
   sphinx-quickstart on Thu Apr  2 12:01:03 2015.
   You can adapt this file completely to your liking, but it should at least
   contain the root `toctree` directive.

Welcome to ECGkit's documentation!
==================================
This toolbox is a collection of Matlab tools that I used, adapted or developed during my PhD and post-doc work with the `Besicos group at University of Zaragoza <http://diec.unizar.es/~laguna/personal/>`__, Spain and at the `National Technological University <http://www.electron.frba.utn.edu.ar/>`__ of Buenos Aires, Argentina. The ECG-kit has tools for reading, processing and presenting results, as you can see in the `documentation <http://ecg-kit.readthedocs.org/en/master/>`__ or in these demos on `Youtube <https://www.youtube.com/watch?v=8lJtkGhrqFw&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=1>`__.

The main feature of the this toolbox is the possibility to use several popular algorithms for ECG processing, such as:

* Algorithms from Physionet's `WFDB software package <http://physionet.org/physiotools/wfdb.shtml>`__
* QRS detectors, such as `gqrs <http://www.physionet.org/physiotools/wag/gqrs-1.htm>`__, `wqrs <http://www.physionet.org/physiotools/wag/gqrs-1.htm>`__, `wavedet <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__, `ecgpuwave <http://www.physionet.org/physiotools/ecgpuwave/>`__, `Pan & Tompkins <http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&arnumber=4122029>`__, `EP limited <http://www.eplimited.com/confirmation.htm>`__
* `Wavedet ECG delineator <http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf>`__
* Pulse wave detectors as `wabp <http://www.physionet.org/physiotools/wag/wabp-1.htm>`__ and `wavePPG <http://dx.doi.org/10.1109/JBHI.2013.2267096>`__
* `a2hbc <https://code.google.com/p/a2hbc/>`__ and `EP limited <http://www.eplimited.com/confirmation.htm>`__ heartbeat classifiers.
* And other scritps for inspecting, correcting and reporting all these results.

with the same application programmer interface (API) directly in Matlab, under Windows, Linux or Mac. The kit also implements a recording interface which allows processing several ECG formats, such as HL7-aECG, MIT, ISHNE, HES, Mortara, and AHA, of arbitrary recording size (the record so far is a 1 week recording of 3 leads, sampled at 500 Hz).

.. image:: ex_ABP_PPG_Registro_01M_full_Pagina_05.png

.. image:: QRS_corrector.PNG

.. image:: 208_full_14.png

.. image:: 1.png

.. image:: 2.png

This kit also includes many open-source projects such as `WFDB Toolbox for MATLAB and Octave <http://physionet.org/physiotools/matlab/wfdb-app-matlab/>`__ from `Physionet <http://physionet.org/>`__, `PRtools <http://prtools.org/>`__, `Libra <https://wis.kuleuven.be/stat/robust/LIBRA>`__, `export_fig <http://undocumentedmatlab.com/blog/export_fig>`__ from `undocumented Matlab <http://undocumentedmatlab.com/>`__, and other open-source scripts that have their proper references to the original projects or authors.

Voluntary contributions
-----------------------
Many thanks to Andrés Demski from UTN who helped to this project before he learned how to use it. To **all** the friends in Zaragoza, Porto and Lund, but in special to the ones closest to the project:

* Pablo Laguna, Juan Pablo Martínez, Rute Almeida and Juan Bolea, for the wavedet ECG delineator and many parts of the Biosig browser project that were adapted to this project.
* Jesús Lázaro and Eduardo Gil for the PPG / ABP pulse detection code.

Involuntary contributions
-------------------------
The acknowledgements also goes to all these people, important in many ways to the fulfilment of this project

* George Moody, Wei Zong, Ikaro Silva, for all the software of `Physionet <http://physionet.org/>`__.
* Reza Sameni, for his `Open-Source ECG Toolbox (OSET) <http://www.oset.ir>`__
* Bob Duin and all the team behind `PRtools <http://prtools.org/>`__
* Yair Altman from `undocumented Matlab <http://undocumentedmatlab.com/>`__
* Diego Armando Maradona for `this <https://github.com/marianux/ecg-kit/blob/master/common/genio_inspirador.jpeg?raw=true>`__.

.. toctree::
   :maxdepth: 4
   :titlesonly:
   :includehidden:
   :hidden:

   Getting started <getting_started>
   Examples <first_example>
   Accessing signals  <ECGwrapper>
   Performing tasks  <ECGtask>
   Result format <results_format>
   Plotting and reporting <reportECG>
   Other functions <other_functions>
   Extending the ECGkit <extensions>
