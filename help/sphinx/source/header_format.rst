
Header Format
=============

.. _header_format:

The header struct variable includes the technical information that describes an ECG recording. The ECG recording is referred as the ECG\_matrix in this document, and is a integer matrix with the raw ADC samples. This structure has a set of mandatory fields that must be included:

	- *freq*, is the sampling frequency of ECG\_matrix signal, or in other words, the reciprocal of the sampling interval or time (in seconds) elapsed between ADC samples.

	- *desc*, description strings about each of the leads/signals. It is a char matrix of [nsig × description_length], where description_length is the length of the longer description string.

	- *nsamp* is the number of samples of ECG\_matrix, or the size of the first matrix dimension (rows).

	- *nsig* is the amount of leads or signals of ECG\_matrix, or the size of the second matrix dimension (columns).

	- *gain* is a double precision vector of [nsig × 1] with the gain of each lead ( ADCsamples / μV ).

	- *adczero* is a double precision vector of [nsig × 1] with the offset of each lead in ADC samples.
	
	
Also other fields described in the `Physionet header <http://www.physionet.org/physiotools/wag/header-5.htm>`__ are accepted by the kit, check the source code for more details.
