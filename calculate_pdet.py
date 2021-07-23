# -*- coding: utf-8 -*-
"""
@author: Doğa Veske
"""

import numpy
from scipy.integrate import quad
from numpy import pi
from scipy.optimize import least_squares
import time
from lal import antenna

t=time.time()

H0=67.6 #Hubble constant in km/s/Mpc
d1=lambda x: ((1+x)**3*0.31+0.69)**-0.5 #local energy density paramters for matter=0.31, for cosmological constant=0.69
def fr(rs): #returns luminosity distance for a redshift
    return 3*10**5*(1+rs)*quad(d1,0,rs)[0]/H0
fr=numpy.vectorize(fr)

n=20 #the number of points in each dimension of the 4 dimensional angle grid (2 sky positions, 2 orbital orientation of the binary)

alpha=numpy.linspace(0,2*pi,num=n) #right ascension - uniform between 0,2pi
psi=numpy.linspace(0,2*pi,num=n) #polarization angle - uniform between 0,2pi
delta=numpy.linspace(-1,1,num=n) #cos(declination) - uniform between -1,+1
delta=numpy.arccos(delta) #declination
io=numpy.linspace(-1,1,num=n) #cos(inclination) - uniform between -1,+1

#values of the antenna patterns of the three detectors on the 4D angle grid 
#Livingston
f=antenna.AntennaResponse('L1',alpha,delta,psi=psi,times=1000000000).plus
fx=antenna.AntennaResponse('L1',alpha,delta,psi=psi,times=1000000000).cross
anl=numpy.outer(f**2,(0.5*(1+io**2))**2)+numpy.outer(fx**2,io**2)
anl=numpy.reshape(anl,(n,n,n,n))
#Hanford
fh=antenna.AntennaResponse('H1',alpha,delta,psi=psi,times=1000000000).plus
fxh=antenna.AntennaResponse('H1',alpha,delta,psi=psi,times=1000000000).cross
anh=numpy.outer(fh**2,(0.5*(1+io**2))**2)+numpy.outer(fxh**2,io**2)
anh=numpy.reshape(anh,(n,n,n,n))
#Virgo
fv=antenna.AntennaResponse('V1',alpha,delta,psi=psi,times=1000000000).plus
fxv=antenna.AntennaResponse('V1',alpha,delta,psi=psi,times=1000000000).cross
anv=numpy.outer(fv**2,(0.5*(1+io**2))**2)+numpy.outer(fxv**2,io**2)
anv=numpy.reshape(anv,(n,n,n,n))

#Load SNRs for each detector computed beforehand
o3l=numpy.load('SNRo3L1.npy')
o3h=numpy.load('SNRo3H1.npy')
o3v=numpy.load('SNRo3V1.npy')

al=numpy.zeros(100) #compute the alpha(m2/m1) function that characterizes the change in SNR due to redshift for each mass ratio
for j in range (10,100):
    a=numpy.zeros(j)
    for i in range (10,j):
        a[i]=o3l[int(100/j*i)][i]
    def ff(x):
        return (numpy.linspace(10,j,len(a)-10)**(x[0])*x[1]-a[10:])
    res_1 = least_squares(ff, numpy.array([1,a[-1]/j**1]))
    al[j]=res_1.x[0]

pdet=numpy.zeros((100,100))
zz=numpy.linspace(0,3.5,num=100000)
rr=fr(zz) #calculate luminosity distances corresponding to redshifts up to 3.5

pr=rr**2/(1+zz)**4
r=numpy.random.choice(rr,size=300,p=pr/numpy.sum(pr)) #sample the distance points according to r^2/(1+z)^4 distribution
r=r[numpy.argsort(r)]
z=numpy.array([zz[numpy.argsort(abs(ra-rr))[0]] for ra in r]) #find corresponding redshifts

def c(i,j):#compute the detection probability with the HLV network at O3 sensitivity with detection threshold POWER SNR=64 (=amplitude SNR=8). Marginalized over the (20^4 x 300) points in the angle and distance grid
    return numpy.sum((1000**2*numpy.outer((1+z)**al[int(numpy.round(100*j/i))-1]/r**2,(o3l[i,j]*anl+o3h[i,j]*anh+o3v[i,j]*anv)))>64)/(300*n**4) 


for j in range (5,100): #Consider small mass between [5,99] solar masses, in accordance with the SNR code
    for i in range (j,min(100,10*j+1)): #Consider heavy mass between [m2,min(99,10*m2)] solar masses, in accordance with the SNR code
        pdet[i][j]=c(i,j)

print(time.time()-t)#print elapsed time
numpy.save('hlvo3.npy',pdet) #name of the file containing the detection probability given m1,m2
