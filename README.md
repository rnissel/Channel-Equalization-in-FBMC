# Channel Equalization in FBMC


This repository simulates an FBMC transmission over a doubly-selective channel and shows the Bit Error Ratio (BER) for different MMSE equalization and interference cancellation methods. 
All figures from R. Nissel et al. [“FBMC-OQAM in Doubly-Selective Channels:
A New Perspective on MMSE Equalization”](https://publik.tuwien.ac.at/files/publik_259771.pdf), IEEE SPAWC, 2017, can be reproduced. 

For more information about FBMC, see [https://github.com/rnissel/FBMC](https://github.com/rnissel/FBMC). Note that my measurements indicate that one-tap equalizers are sufficient in many practical wireless communication scenarios, see my [PhD thesis](http://publik.tuwien.ac.at/files/publik_265168.pdf)  (perfectly time and frequency synchronized). Only in some rare cases, enhanced equalization methods, as presented here, might be useful.


## Requirements
We used Windows 7 (64bit) and Matlab R2013b/2016a, but newer versions (and some older) should also work.

## Reproducible Figures
The figure numbers are the same as in  [“FBMC-OQAM in Doubly-Selective Channels:
A New Perspective on MMSE Equalization”](https://publik.tuwien.ac.at/files/publik_259771.pdf):

* **Figure  1**: 
Please run [`PlotInterferencePower2D.m`](PlotInterferencePower2D.m) after uncommenting line 26.

* **Figure  2**: 
Please run [`PlotInterferencePower2D.m`](PlotInterferencePower2D.m).

* **Figure  3**: 
Just an illustration.


* **Figure  4**: 
Please run [`BER_SISO.m`](BER_SISO.m) after uncommenting line 34. Moreover, to truly reproduce Figure 4, line 31 must be uncommented (increases the simulation time). 

* **Figure  5**: 
Please run [`BER_SISO.m`](BER_SISO.m). To truly reproduce Figure 5, line 31 must be uncommented (increases the simulation time).

* **Figure  6**: 
Please run [`BER_SISO.m`](BER_SISO.m). To truly reproduce Figure 6, line 31 must be uncommented (increases the simulation time).


* **Figure  7**: 
Please run [`BER_MIMO.m`](BER_SISO.m). To truly reproduce Figure 7, the lines 31-32 must be uncommented (increases the simulation time).



## Please Cite Our Paper

    @inproceedings{Nissel2017spawcA,
		author    = {R. Nissel and M. Rupp and R. Marsalek},
		booktitle = {IEEE International Workshop on Signal Processing Advances in Wireless Communications (SPAWC)},
		title     = {{FBMC-OQAM} in doubly-selective channels: A new perspective on {MMSE} equalization},
		year 	  = {2017},
		pages 	  = {1-5}, 
		doi 	  = {10.1109/SPAWC.2017.8227806},
		month 	  = {July},
	}


## References
- R. Nissel, M. Rupp, and R. Marsalek, [“FBMC-OQAM in Doubly-Selective Channels:
A New Perspective on MMSE Equalization”](https://publik.tuwien.ac.at/files/publik_260162.pdf), IEEE International Workshop on Signal Processing Advances in Wireless Communications (SPAWC), Sapporo, Japan, July 2017.
- R. Nissel, S. Schwarz, and M. Rupp, [“Filter bank multicarrier modulation schemes for future mobile communications”](https://publik.tuwien.ac.at/files/publik_260162.pdf), IEEE Journal on Selected Areas in Communications, vol. 35, no. 8, pp. 1768–1782, 2017.
- R. Nissel, [“Filter bank multicarrier modulation for future wireless systems”](http://publik.tuwien.ac.at/files/publik_265168.pdf), Dissertation, TU Wien, 2017.



