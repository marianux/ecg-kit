# Welcome to the ecg-kit !

This toolbox is a collection of Matlab tools used, adapted or developed by me during my PhD and post-doc work in the [Besicos group at University of Zaragoza](http://diec.unizar.es/~laguna/personal/), Spain and the [National Technological University](http://www.electron.frba.utn.edu.ar/) of Buenos Aires, Argentina. The ECG-kit has tools for reading, processing and presenting results, as you can see in the [documentation](http://ecg-kit.readthedocs.org/en/latest/index.html) or in these demos on [Youtube](https://www.youtube.com/watch?v=8lJtkGhrqFw&list=PLlD2eDv5CIe9sA2atmnb-DX48FIRG46z7&index=1).

The main feature of the this toolbox is the possibility to use several popular algorithms for ECG processing, such as:

* Algorithms from Physionet's [WFDB software package](http://physionet.org/physiotools/wfdb.shtml)
* QRS detectors, such as gqrs, wqrs, ecgpuwave, Pan & Tompkins, [EP limited](http://www.eplimited.com/confirmation.htm)
* [Wavedet ECG delinator](http://diec.unizar.es/~laguna/personal/publicaciones/wavedet_tbme04.pdf)
* [a2hbc heartbeat classifier](https://code.google.com/p/a2hbc/)
* And other scritps for inspecting, correcting and reporting all these results. 

with the same application programmer interface (API) directly in Matlab, under Windows, Linux or Mac. The kit also implements a recording interface which allows processing several ECG formats, such as MIT, ISHNE, HES, Mortara, and AHA, of arbitrary recording size (the record so far is a 1 week recording of 3 leads, sampled at 500 Hz).

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/ex_ABP_PPG_Registro_01M_full_Pagina_05.png)

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/QRS_corrector.PNG)

![Image not found](http://ecg-kit.readthedocs.org/en/latest/_images/208_full_14.png)

This kit also includes many open-source projects such as [WFDB Toolbox for MATLAB and Octave](http://physionet.org/physiotools/matlab/wfdb-app-matlab/) from [Physionet](http://physionet.org/), [PRtools](http://prtools.org/), [Libra](https://wis.kuleuven.be/stat/robust/LIBRA), [export_fig](http://undocumentedmatlab.com/blog/export_fig) from [undocumented Matlab](http://undocumentedmatlab.com/), and other open-source scripts that have their proper references to the original projects or authors.
