# Glamor
This program allows calculation of static and dynamic REDOR curves in solid-state NMR. It uses analytical formula and thus
allows to calculate the interaction of one observed spin to thousands of recoupled spins and can also calculate the average
of different observed spins (neglecting homonuclear dipolar coupling).
Details and theoretical background of the program can be found here:
http://dx.doi.org/10.1016/j.ssnmr.2012.10.001
https://www.chemie-biologie.uni-siegen.de/ac/jsadg/publications/files/053-manuscript.pdf

It is compiled with: cd source ; make glamor
Copy the binary "glamor" wherever you need it. You may have to edit the Makefile in the source directory.

The required inputfiles can be found in the folder "inputfiles". Required files:
time-sequence.dat
control-sequence.dat
input.xyz

They have to be modified according to the system that has to be calculated
and then copied into a separate folder. Start the program for example by typing:
cd inputfiles
../source/glamor

Citations are appreciated.

Authors: Janin Glänzer, Vinicius Ribeiro Celinski, Johannes Weber, Jörn Schmedt auf der Günne
