# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:39:28 2017

@author: mjh250
"""

import numpy as np
from scipy.interpolate import interp1d
from scipy.optimize import curve_fit
from scipy import asarray as ar, exp
import matplotlib.pyplot as plt


class FringeVisibilityAnalysis:
    def __init__(self):
        return

    def envelope(self, data):
        q_u = np.zeros(data.shape)
        q_l = np.zeros(data.shape)

        # Prepend the first value of (data) to the interpolating values. This
        # forces the model to use the same starting point for both the upper
        # and lower envelope models.
        u_x = [0, ]
        u_y = [data[0], ]

        l_x = [0, ]
        l_y = [data[0], ]

        # Detect peaks and troughs and mark their location in u_x,u_y,l_x,l_y
        # respectively.
        for k in xrange(1, len(data)-1):
            if ((np.sign(data[k]-data[k-1]) == 1) and
               (np.sign(data[k]-data[k+1]) == 1)):
                u_x.append(k)
                u_y.append(data[k])

            if ((np.sign(data[k]-data[k-1]) == -1) and
               ((np.sign(data[k]-data[k+1])) == -1)):
                l_x.append(k)
                l_y.append(data[k])

        # Append the last value of (data) to the interpolating values. This
        # forces the model to use the same ending point for both the upper
        # and lower envelope models.
        u_x.append(len(data)-1)
        u_y.append(data[-1])

        l_x.append(len(data)-1)
        l_y.append(data[-1])

        # Fit suitable models to the data. Here I am using cubic splines,
        # similarly to MATLAB.
        u_p = interp1d(u_x, u_y, kind='cubic',
                                      bounds_error=False, fill_value=0.0)
        l_p = interp1d(l_x, l_y, kind='cubic',
                                      bounds_error=False, fill_value=0.0)

        # Evaluate each model over the domain of (data)
        for k in xrange(0, len(data)):
            q_u[k] = u_p(k)
            q_l[k] = l_p(k)

        # Return upper and lower envelopes
        return q_u, q_l

if __name__ == '__main__':
    scan_analyzer = FringeVisibilityAnalysis()
    filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\"
                "2017_10_20\\Supercontinuum650to650nm\\CrossectionsAt70\\")
    dataa = []
    datab = []
    datac = []
    data = []
    micrometers = [1700, 1800, 1900, 2000, 2100, 2200]
    trimmedData = []
    dataCenterPixel = 800
    dataPixelWidth = 75
    upperEnvelopes = []
    lowerEnvelopes = []
    visibilityArrays = []
    visibility = []
    j = 1

    for i in micrometers:
        dataa.append(np.loadtxt(filepath+"Micrometer"+str(i)+"microns_a.txt"))
        datab.append(np.loadtxt(filepath+"Micrometer"+str(i)+"microns_b.txt"))
        datac.append(np.loadtxt(filepath+"Micrometer"+str(i)+"microns_c.txt"))
        data.append(np.mean([dataa[-1], datab[-1], datac[-1]], axis=0))
        upperEnvelopes.append(scan_analyzer.envelope(
                data[-1][dataCenterPixel-dataPixelWidth:
                         dataCenterPixel+dataPixelWidth])[0])
        lowerEnvelopes.append(scan_analyzer.envelope(
                data[-1][dataCenterPixel-dataPixelWidth:
                         dataCenterPixel+dataPixelWidth])[1])
        trimmedData.append(data[-1][dataCenterPixel-dataPixelWidth:
                           dataCenterPixel+dataPixelWidth])
        visibilityArrays.append((upperEnvelopes[-1]-lowerEnvelopes[-1]) /
                                (upperEnvelopes[-1]+lowerEnvelopes[-1]))
        visibility.append(max(visibilityArrays[-1]))

        plt.figure(j)
        plt.plot(trimmedData[-1])
        plt.plot(upperEnvelopes[-1])
        plt.plot(lowerEnvelopes[-1])
        j = j + 1

    x = ar(micrometers)
    y = ar(visibility)
    yweight = ar(visibility/sum(visibility))
    n = len(x)                               # the number of data
    mean = sum(x*yweight)                    # note this correction
    sigma = sum(yweight*(x-mean)**2)/n       # note this correction

    def gaus(x, a, x0, sigma):
        return a*exp(-(x-x0)**2/(2*sigma**2))

    popt, pcov = curve_fit(gaus, x, y, p0=[0.5, mean, sigma])

    gaussX = np.arange(x[0], x[-1], ((x[-1]-x[0])/100.0))
    plt.figure(j+1)
    plt.plot(x, y, 'b+:', label='Visibility')
    plt.plot(gaussX, gaus(gaussX, *popt), 'r', label=('Gaussian fit, mean = ' +
             str(mean) + ', sigma = ' + str(sigma)))
    plt.legend()
    plt.xlabel("Micrometer position (microns)")
    plt.ylabel("Peak visibility")

    plt.show()

    print("Finished running script!")
