Mean square flux noise
====

This repository contains the programs needed to calculate mean square flux noise (MSFN) of superconducting SQUID or qubit loops. Binaries of freely available 64-bit superconducting version of FastHenry and InductEx are included. My Matlab code implements the class `iewasher`, which defines properties and methods relevant for a single loop. The program `batchIEW` creates an array of `iewasher` object, where we vary a specific geometric parameter, and `runBatchIEW` controls the execution of the computation. The program `plotIEWarray` nicely plots the results.

I'm just about finished a draft describing the computation, and I'll post it as soon as it's submitted.

Steven M. Anton<br>
Department of Physics<br>
UC Berkeley
