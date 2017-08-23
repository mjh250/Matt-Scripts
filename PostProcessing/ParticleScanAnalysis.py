# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:33:43 2017

@author: mjh250
"""

import h5py
import numpy as np
import matplotlib.pyplot as plt


class ParticleScanAnalysis:
    def __init__(self):
        return

    def getGrayscaleImage(self, image_name, data):
        image = np.array(data[image_name])
        r, g, b = image[:, :, 0], image[:, :, 1], image[:, :, 2]
        # Convert to gray using MATLAB's (NTSC/PAL) implementation of rgb2gray:
        gray = 0.2989 * r + 0.5870 * g + 0.1140 * b
        return gray

    def processImage(self, data_image, bias_image=None, background_image=None,
                     white_reference_image=None):
        if bias_image is None:
            bias_image = np.full_like(data_image, 0)
        if background_image is None:
            background_image = np.full_like(data_image, 0)
        if white_reference_image is None:
            white_reference_image = np.full_like(data_image, 255)
        print(np.shape(data_image))
        print(np.shape(bias_image))
        print(np.shape(background_image))
        processed_image = ((data_image - bias_image - background_image) /
                           (white_reference_image - background_image))
        return processed_image

    def getScanData(self, filepath, scan_number, particle_number):
        file = h5py.File(filepath, 'r')
        dataset = file['/particleScans/scan%s/scan_%s'
                       % (scan_number, particle_number)]
        return dataset

if __name__ == '__main__':
    filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\"
                "2017_08_21\\test_scans_2.hdf5")
    scan_analyzer = ParticleScanAnalysis()
    data = scan_analyzer.getScanData(filepath, 0, 0)
    data_image = scan_analyzer.getGrayscaleImage(
                                    'Infinity3_FirstWhiteLight_Image', data)
    bias_image = scan_analyzer.getGrayscaleImage('Infinity3_Bias_Image', data)
    background_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_FirstBkgndWhiteLight_Image', data)
    img = scan_analyzer.processImage(data_image, bias_image, background_image)
    # Plot
    fig1 = plt.figure(figsize=(20, 20))
    plt.imshow(img, cmap='gray')
    plt.show()
