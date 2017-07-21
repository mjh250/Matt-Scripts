# -*- coding: utf-8 -*-
"""
Created on Fri Jul 18 14:31:36 2014

@author: rbowman
"""

import numpy as np
import numpy.linalg
import matplotlib.pyplot as plt
import scipy.signal

try:
    spectra
    peak_z
    correlation_coefficient
    standard_deviation
    wavelengths
    reference
    print "Not reloading data"
except NameError:
    import h5py
    data = h5py.File("summary.hdf5","r")
    scannumber = 1
    spectra = data['particleScanSummaries/scan%d/spectra' % scannumber]
    wavelengths = np.array(spectra.attrs.get('wavelengths'))
    reference = np.array(spectra.attrs.get('reference'))
    spectra = np.array(spectra)
    peak_z = np.array(data['particleScanSummaries/scan%d/peak_z' % scannumber],copy=True)
    correlation_coefficient = np.array(data['particleScanSummaries/scan%d/correlation_coefficient' % scannumber])
    standard_deviation = np.array(data['particleScanSummaries/scan%d/standard_deviation' % scannumber])
    data.close()

def wl_slice(start_wl=-np.Inf, end_wl=np.Inf):
    return slice(np.min(np.where(start_wl < wavelengths)), np.max(np.where(wavelengths < end_wl)))
    
good_wavelengths =  wl_slice(450,850)#this can be used to extract just the good bits

def find_mode(data,axis=None,**kwargs):
    if axis is None:
        hist, bins = np.histogram(data[~np.isnan(data)], **kwargs)
        return (bins[:-1]/2+bins[1:]/2)[hist.argmax()] #return the centre position of the best-filled histogram bin
    else:
        mode = np.zeros(np.delete(data.shape, axis)) #to take the mode along one axis, first delete that axis to get the output shape
        it = np.nditer(mode,flags=['multi_index'],op_flags=['writeonly'])
        while not it.finished:
            it[0] = find_mode(data[it.multi_index[0:axis]+(Ellipsis,)+it.multi_index[axis:]], axis=None, **kwargs)
            it.iternext()
        return mode
#Allow for slight wobble in Z by fixing the mean position for 550-600nm as zero, then find the RMS deviation from the mean 
ca_smoothed = scipy.signal.convolve(peak_z,scipy.signal.gaussian(13,3)[np.newaxis,:])
ca_distance = np.std(ca_smoothed[:,wl_slice(500,750)] - find_mode(ca_smoothed[:,wl_slice(550,600),np.newaxis],axis=1,bins=100), axis=1)

def find_mode_and_width(scores, fraction=0.9, **kwargs):
    hist, binedges = np.histogram(scores, **kwargs)
    bincentres = (binedges[:-1]+binedges[1:])/2
    mode_index = hist.argmax()
    width = 0
    while np.sum(hist[max(mode_index-width,0):min(mode_index+width,len(hist))]) < fraction * len(scores):
        width +=1
    return bincentres[mode_index], max(binedges[min(mode_index+width+1,len(binedges)-1)]-bincentres[mode_index], 
                                       bincentres[mode_index]-binedges[max(mode_index-width,0)])

mean_r = np.mean(correlation_coefficient[:,ca_filter_wavelengths],axis=1)
"""Method 1: draw an ellipse"""
mean_r_mode, mean_r_width = find_mode_and_width(mean_r,0.75,bins=200)
ca_distance_mode, ca_distance_width = find_mode_and_width(ca_distance,0.7,bins=200)
combined_distance = np.sqrt(((mean_r-mean_r_mode)/mean_r_width)**2 + ((ca_distance-ca_distance_mode)/ca_distance_width)**2)
good_spectra = np.where(combined_distance<1)[0] #pick points outside the threshold

"""
plt.plot(mean_r, ca_distance, 'k.')
angles = np.arange(0,np.pi*2,np.pi/100)
plt.plot([mean_r_mode + mean_r_width*np.cos(a) for a in angles],[ca_distance_mode + ca_distance_width*np.sin(a) for a in angles],'r-')
"""


#Principal Component Analsyis - with much help from 
#http://stackoverflow.com/questions/1730600/principal-component-analysis-in-python
"""
plt.figure()
plt.ioff()
for i in good_spectra:
    plt.plot(wavelengths[good_wavelengths], peak_z[i,good_wavelengths] - np.mean(peak_z[i,wl_slice(550,600)]),'+')
plt.draw()
plt.ion()
"""
#First, we should condition the data
#When working with chromatic aberration, let's bin everything outside the good
#spectral range, so

A = np.array(spectra[good_spectra,good_wavelengths]) #limit to the wavelengths where there's signal
mask = np.isnan(A)                   #fill in NaNs by interpolating nearest non-NaN values
A[mask] = np.interp(np.flatnonzero(mask), np.flatnonzero(~mask), A[~mask])
A -= np.nanmean(A, axis=0) #subtract the means

print "Running SVD (may take some time!)"
U, s, V = np.linalg.svd(A) 

#To get the data in principal component space, we do:
#for the first measurement
principal_components = s * np.dot(V,A.T).T

"""
f, ax= plt.subplots(10,1,sharex=True) #plot the principal components
for i, a in enumerate(ax):
     a.plot(wavelengths[good_wavelengths],V[i,:])
"""