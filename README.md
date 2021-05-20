# observation-bias-gw

This repository has example scripts for calculating the observation bias in gravitational-wave detections for non-spinning black holes. 

**If you use them please cite as ...**

Feel free to use and change the parameters in the code for increasing/decreasing precision or for changing the sensitivity used for the detectors. Scripts have detailed explanation in them with in line comments.

Required python packages: numpy, scipy, h5py, gwsurrogate, lal (for antenna patterns), time (for timing the runs). Install these first. Also make sure you load the NRHybSur3dq8 waveforms by ```gwsurrogate.catalog.pull('NRHybSur3dq8')``` before your first use of the gwsurrogate package.


calculateSNR.py calculates power signal-to-noise ratios generated in LIGO Hanford, LIGO Livingston and Virgo for a face-on oriented binary (inclination angle=0) ignoring the antenna factor which is 1 Gpc away and in the absence of redshift,on a grid of black hole masses. The current parameters for the code is 5<=m2<=99 and m2<=m1<=min(99,10\*m2) for integer m1 and m2. The code uses the PSD of the GW190412 event which is NOT provided in this repository. To use it download GW190412.tar from https://dcc.ligo.org/LIGO-P2000223/public/. Only the GW190412.h5 file is necessary. The output of this script is a .npy file for each detector. The computation takes about 20 mins for each detector with 4 CPU cores at 2.8 GHz (not all the cores are fully used, so this is a conservative estimate). The outputs for the configuration here are also given which are SNRo3H1.npy, SNRo3L1.npy, SNRo3V1.npy. These have the values for the E_0 in the paper (see Eq. 4 or 9).

calculate_pdet takes the outputs of calculateSNR.py and calculates the detection probability for the mass grid with power SNR threshold 64 (=amplitude SNR 8) assuming uniform distribution in sky, uniform orientation for binary's plane and a r^2/(1+z)^4 distribution in distance which corresponds to uniform distribution in comoving volume with a constant local event rate density in the presence of redshift. The marginalization is done over 20 points in each angular dimension and 300 points in distance (total 48\*10^6) for each mass combination. Overall, with the evaluation in the mass space it considers order of 10^11 points. The computation takes about 1 hour with 4 CPU cores at 2.8 GHz (for this estimate CPU usage was around 20% so it the calculation was probably done with only one core, so this is probably a very conservative estimate for 4 cores). The output for the configuration here is also given in the file hlvo3.npy. It has the values for P(D|m1,m2) (like in Eq 2. in the paper but without marginalization over m2).

The two scripts complete their runs in about 2 hours for a 3 detector network considering about 10^11 points in the 7 dimensional parameter space with 4 CPU cores at 2.8 GHz.

Questions to be sent to Doğa Veske dv2397@columbia.edu.
