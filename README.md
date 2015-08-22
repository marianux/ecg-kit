Check the web page of this project at 

# http://marianux.github.io/ecg-kit/

# Welcome to the ecg-kit !

This toolbox is a collection of Matlab tools that I used, adapted or developed during my PhD and post-doc work with the [Besicos group at University of Zaragoza](http://diec.unizar.es/~laguna/personal/), Spain and at the [National Technological University](http://www.electron.frba.utn.edu.ar/) of Buenos Aires, Argentina. The ECG-kit has tools for reading, processing and presenting results, as you can see in the [documentation](http://ecg-kit.readthedocs.org/en/master/) or in these demos on [Youtube](https://www.youtube.com/watch?v=8lJtkGhrqFw&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=1).

The main feature of the this toolbox is the possibility to use several popular algorithms for ECG processing, such as:

* Algorithms from Physionet's [WFDB software package](http://physionet.org/physiotools/wfdb.shtml)
* QRS detectors, such as [gqrs](http://www.physionet.org/physiotools/wag/gqrs-1.htm), [wqrs](http://www.physionet.org/physiotools/wag/gqrs-1.htm), [wavedet](http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf), [ecgpuwave](http://www.physionet.org/physiotools/ecgpuwave/), [Pan & Tompkins](http://ieeexplore.ieee.org/xpl/articleDetails.jsp?reload=true&arnumber=4122029), [EP limited](http://www.eplimited.com/confirmation.htm)
* [Wavedet ECG delineator](http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf)
* Pulse wave detectors as [wabp](http://www.physionet.org/physiotools/wag/wabp-1.htm) and [wavePPG](http://dx.doi.org/10.1109/JBHI.2013.2267096)
* [a2hbc](https://code.google.com/p/a2hbc/) and [EP limited](http://www.eplimited.com/confirmation.htm) heartbeat classifiers.
* And other scritps for inspecting, correcting and reporting all these results. 

with the same application programmer interface (API) directly in Matlab, under Windows, Linux or Mac. The kit also implements a recording interface which allows processing several ECG formats, such as MIT, ISHNE, HES, Mortara, and AHA, of arbitrary recording size (the record so far is a 1 week recording of 3 leads, sampled at 500 Hz).

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/ex_ABP_PPG_Registro_01M_full_Pagina_05.png)

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/QRS_corrector.PNG)

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/208_full_14.png)

This kit also includes many open-source projects such as [WFDB Toolbox for MATLAB and Octave](http://physionet.org/physiotools/matlab/wfdb-app-matlab/) from [Physionet](http://physionet.org/), [PRtools](http://prtools.org/), [Libra](https://wis.kuleuven.be/stat/robust/LIBRA), [export_fig](http://undocumentedmatlab.com/blog/export_fig) from [undocumented Matlab](http://undocumentedmatlab.com/), and other open-source scripts that have their proper references to the original projects or authors.

## Voluntary contributions
Many thanks to Andrés Demski from UTN who helped to this project before he learned how to use it. To **all** the friends in Zaragoza, Porto and Lund, but in special to the ones closest to the project:

* Pablo Laguna, Juan Pablo Martínez, Rute Almeida and Juan Bolea, for the wavedet ECG delineator and many parts of the Biosig browser project that were adapted to this project. 
* Jesús Lázaro and Eduardo Gil for the PPG / ABP pulse detection code.
* Li-wei Lehman from Physionet/MIT helped a lot in testing the first versions of the kit.

## Involuntary contributions
The acknowledgements also goes to all these people, important in many ways to the fulfilment of this project

* George Moody, Wei Zong, Ikaro Silva, for all the software of [Physionet](http://physionet.org/).
* Reza Sameni, for his [Open-Source ECG Toolbox (OSET)](http://www.oset.ir)
* Bob Duin and all the team behind [PRtools](http://prtools.org/)
* Yair Altman from [undocumented Matlab](http://undocumentedmatlab.com/)
* Diego Armando Maradona for [this](https://github.com/marianux/ecg-kit/blob/master/common/genio_inspirador.jpeg?raw=true).
