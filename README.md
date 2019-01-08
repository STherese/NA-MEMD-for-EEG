# NA-MEMD-for-EEG
Material related to "Unmixing oscillatory brain activity by EEG source localization and empirical mode decomposition". 

Implementation of the noise-assisted multivariate empiral mode decomposition (NA-MEMD) method, builds on code from Rehman & Mandic see:
http://www.commsp.ee.ic.ac.uk/~mandic/research/emd.htm, N. Rehman and D. P. Mandic, "Multivariate Empirical Mode 
Decomposition," Proceedings of the Royal Society A, vol. 466, no. 2117, pp. 1291-1302, 2010. http://doi.org/10.1098/rspa.2009.0502
and 
N. Rehman, D. P. Mandic, Filter Bank Property of Multivariate Empirical Mode Decomposition. IEEE Transactions on Signal Processing, 59(5), 
pp. 2421–2426, 2011. http://doi.org/10.1109/TSP.2011.2106779.

The following is required to run the code:
- The MEMD toolbox from http://www.commsp.ee.ic.ac.uk/~mandic/research/emd.htm
- A forward model/lead field matrix must be provided (dimensions: channels x vertices)
