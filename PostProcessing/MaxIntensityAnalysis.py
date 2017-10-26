# -*- coding: utf-8 -*-
"""
Created on Thu Oct 12 11:39:28 2017

@author: mjh250
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt
import bottleneck
# import datetime


class MaxIntensityAnalysis:
    def __init__(self):
        return

    def getScanDataFile(self, filepath):
        datafile = h5py.File(filepath, 'a')
        return datafile

    def getScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/scan%s/particle%s'
                           % (scan_number, particle_number)]
        return dataset

if __name__ == '__main__':
    scan_analyzer = MaxIntensityAnalysis()
    filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\"
                "2017_09_19\\processedData.hdf5")
    processed_filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\"
                          "mjh250\\2017_09_19\\")
    data = scan_analyzer.getScanDataFile(filepath)
    scan_number_list = []
    particle_number_list = []
    infinity3_averaged_maxima_list = []
    white0order_averaged_maxima_list = []
    laser0order_averaged_maxima_list = []
    for name in data:
        if name.startswith('scan'):
            scan_number_list.append(len(scan_number_list))
        for particle_name in data[name]:
            if particle_name.startswith('particle'):
                particle_number_list.append(len(particle_number_list))

    # --- Iterate through particle scans and process files
    for scan_number in scan_number_list:
        for particle_number in particle_number_list:
            # --- Print name of file being analyzed to track progress
            print('Processing scan%s/particle%s'
                  % (scan_number, particle_number))

            # --- GET SCAN DATA & INITIALIZE PROCESSED FILE ---
            datafile = scan_analyzer.getScanDataSet(data, scan_number,
                                                    particle_number)
            processedData = h5py.File(
                                  processed_filepath +
                                  "IntensityAnalysisData.hdf5", "a")
            pdata = processedData.create_group(
                          'scan%s/particle%s' % (scan_number, particle_number))

            # --- INFINITY 3 FIRST IMAGE ---
            infinity3_maxima_to_average = 10
            infinity3_image = np.array(
                                   datafile['Infinity3_First_Processed_Image'])
            infinity3_z = -bottleneck.partition(
                    -infinity3_image.flatten(),
                    infinity3_maxima_to_average)[:infinity3_maxima_to_average]
            infinity3_averaged_maxima_list.append(np.mean(infinity3_z))

            # --- RAMAN WHITE LIGHT IMAGE ---
            white0order_maxima_to_average = 10
            white0order_image = np.array(
                          datafile['Raman_White_Light_0Order_Processed_Image'])
            white0order_z = -bottleneck.partition(
                -white0order_image.flatten(),
                white0order_maxima_to_average)[:white0order_maxima_to_average]
            white0order_averaged_maxima_list.append(np.mean(white0order_z))
            
            # --- RAMAN LASER LIGHT IMAGE ---
            laser0order_maxima_to_average = 10
            laser0order_image = np.array(
                          datafile['Raman_Laser_0Order_Processed_Image'])
            laser0order_z = -bottleneck.partition(
                -laser0order_image.flatten(),
                laser0order_maxima_to_average)[:laser0order_maxima_to_average]
            laser0order_averaged_maxima_list.append(np.mean(laser0order_z))

    # CREATE HISTOGRAM INFINITY3
    fig1 = plt.figure()
    hist, bin_edges, patches = plt.hist(infinity3_averaged_maxima_list, 10)
    plt.xlabel('Average of the 10 brightest pixels')
    plt.ylabel('Number of occurences')
    plt.title('Intensities of the brightest features in the studied particles'
              ' (Infinity3 Images)')
    plt.grid(True)
    plt.show()
    print(hist)
    print(bin_edges)

    # CREATE HISTOGRAM WHITE LIGHT 0ORDER
    fig2 = plt.figure()
    hist, bin_edges, patches = plt.hist(white0order_averaged_maxima_list, 10)
    plt.xlabel('Average of the 10 brightest pixels')
    plt.ylabel('Number of occurences')
    plt.title('Intensities of the brightest features in the studied particles'
              ' (0 Order Raman White Light Images)')
    plt.grid(True)
    plt.show()
    print(hist)
    print(bin_edges)
    
    # CREATE HISTOGRAM LASER LIGHT 0ORDER
    fig3 = plt.figure()
    hist, bin_edges, patches = plt.hist(laser0order_averaged_maxima_list, 10)
    plt.xlabel('Average of the 10 brightest pixels')
    plt.ylabel('Number of occurences')
    plt.title('Intensities of the brightest features in the studied particles'
              ' (0 Order Raman Laser Light Images)')
    plt.grid(True)
    plt.show()
    print(hist)
    print(bin_edges)

    processedData.close()
    data.close()
    print("Finished running script!")
