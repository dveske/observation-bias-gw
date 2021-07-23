# -*- coding: utf-8 -*-
"""
@author: Doğa Veske
"""
import numpy
import gwsurrogate
import h5py
import time
import xml.etree.ElementTree as ET

tic=time.time()

sur=gwsurrogate.LoadSurrogate('NRHybSur3dq8')#name of the waveform families used

snr=numpy.zeros((100,100)) #Consider integer masses until 99 solar masses

d=5*10**-4 #sampling period in seconds, appropriate for signal below 1 kHz

f1=h5py.File('GW190412.h5','r') #name of the event that we use its PSD
for ifo in ('H1','L1','V1'): #compute SNR generated in each detector separately
    xp=f1['C01:SEOBNRv4P']['psds'][ifo][:,0] #Note that these are ONE SIDED PSDs.
    yp=f1['C01:SEOBNRv4P']['psds'][ifo][:,1]

    x=numpy.arange(1000,dtype=int) 
    y=numpy.zeros(1000)
    for i in x:#reduce the PSD's sampling intervals to 1 Hz as we don't need too much detail for this calculation
        y[i]=yp[numpy.argsort(abs(i-xp))[0]]

    for m2 in range (5,100): #Consider small mass between [5,99] solar masses
	
         for m1 in range (m2,min(100,10*m2+1)): #Consider heavy mass between [m2,min(99,10*m2)] solar masses. Mass ratios >10 are not applicable
            t,h,dyn=sur(m1/m2,[0,0,0],[0,0,0],dt=d,f_low=15,units='mks',dist_mpc=1000,M=m1+m2) #the gw waveform starting at 15 Hz for an event with specified masses, located 1 Gpc away, face on orientation, no effect of redshift, no spin - all modes separately
            t1,h1,dyn1=sur(m1/m2,[0,0,0],[0,0,0],dt=d,f_low=15,units='mks',dist_mpc=1000,M=m1+m2,inclination=0,phi_ref=0) #same gw waveform, but with combined modes for inclination=0
            newft=numpy.zeros(len(t))
            count=0
            for i in numpy.fft.fftfreq(t.shape[-1],d=d): #align the sampled frequencies of the waveform and the PSD, for multiplication in frequency domain
                i=abs(i)
                newft[count]=(y[numpy.argsort(abs(i-x))[0]])
                count+=1
            snr[m1][m2]=(2*numpy.abs(numpy.sum(numpy.fft.fft(h1.real)*numpy.fft.fft(h[2,2].real).conj()/newft))/len(t)*d*0.5*(5/numpy.pi)**0.5)**2/(2*numpy.abs(numpy.sum(numpy.fft.fft(h[2,2].real)*numpy.fft.fft(h[2,2].real).conj()/newft))/len(t)*d*0.5*(5/numpy.pi)**0.5*0.5*(5/numpy.pi)**0.5) #compute SNR for plus polarization. Factor of 2 is for correcting the one sided PSD. For 0 inclination, SNR is same for + and x polarizations.
    print(ifo,time.time()-tic) #print the elapsed time
    numpy.save('SNRo3'+ifo+'.npy',snr)#name of the file containing the SNRs for each mass combination
