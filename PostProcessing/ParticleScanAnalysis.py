# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:33:43 2017

@author: mjh250
"""

import h5py
import numpy as np
# import matplotlib.pyplot as plt
import datetime


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
        processed_image = ((data_image - bias_image - background_image) /
                           (white_reference_image - background_image))
        return processed_image

    def processAndor(self, data_image, bias_image=None, background_image=None):
        if bias_image is None:
            bias_image = np.full_like(data_image, 0)
        if background_image is None:
            background_image = np.full_like(data_image, 0)
        processed_image = data_image - bias_image - background_image
        return processed_image

    def getScanDataFile(self, filepath):
        datafile = h5py.File(filepath, 'a')
        return datafile

    def getScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/particleScans/scan%s/scan_%s'
                           % (scan_number, particle_number)]
        return dataset

    def reZeroImage(self, image):
        reZeroedImage = image - min(image.flatten())
        return reZeroedImage

if __name__ == '__main__':
    scan_analyzer = ParticleScanAnalysis()
    filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\"
                "2017_10_05\\test_scans - Copy.hdf5")
    processed_filepath = ("C:\\Users\\mjh250\\Documents\\Local mjh250\\"
                          "mjh250\\2017_10_05\\")
    data = scan_analyzer.getScanDataFile(filepath)
    scan_number_list = []
    particle_number_list = []
    for name in data['/particleScans']:
        if name.startswith('scan'):
            scan_number_list.append(len(scan_number_list))
    for name in data['/particleScans/scan0']:
        if name.startswith('scan_'):
            particle_number_list.append(len(particle_number_list))

    # --- Iterate through particle scans and process files
    for scan_number in scan_number_list:
        for particle_number in particle_number_list:
            # --- Print name of file being analyzed to track progress
            print('Processing scan%s/scan_%s' % (scan_number, particle_number))

            # --- GET SCAN DATA & INITIALIZE PROCESSED FILE ---
            datafile = scan_analyzer.getScanDataSet(data, scan_number,
                                                    particle_number)
            processedData = h5py.File(
                                  processed_filepath+"processedData.hdf5", "a")
            pdata = processedData.create_group(
                          'scan%s/particle%s' % (scan_number, particle_number))

            # --- FIRST WHITE LIGHT IMAGE ---
            data_image = scan_analyzer.getGrayscaleImage(
                                   'Infinity3_FirstWhiteLight_Image', datafile)
            bias_image = scan_analyzer.getGrayscaleImage(
                                           'Infinity3_Bias_Image', datafile)
            background_image = scan_analyzer.getGrayscaleImage(
                              'Infinity3_FirstBkgndWhiteLight_Image', datafile)
            img = scan_analyzer.processImage(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                                  "Infinity3_First_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- SECOND WHITE LIGHT IMAGE ---
            data_image = scan_analyzer.getGrayscaleImage(
                                  'Infinity3_SecondWhiteLight_Image', datafile)
            bias_image = scan_analyzer.getGrayscaleImage(
                                              'Infinity3_Bias_Image', datafile)
            background_image = scan_analyzer.getGrayscaleImage(
                       'Infinity3_SecondWhiteLight_atBkgndLoc_Image', datafile)
            img = scan_analyzer.processImage(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                                 "Infinity3_Second_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- LASER BEAM IMAGE ---
            data_image = scan_analyzer.getGrayscaleImage(
                                        'Infinity3_Laser_Beam_Image', datafile)
            bias_image = scan_analyzer.getGrayscaleImage(
                                              'Infinity3_Bias_Image', datafile)
            img = scan_analyzer.processImage(data_image, bias_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset("Infinity3_Laser_Beam_Processed_Image",
                                        data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- LASER BEAM IMAGE AT BACKGROUND LOCATION ---
            data_image = scan_analyzer.getGrayscaleImage(
                             'Infinity3_Laser_Beam_Image_atBkgndLoc', datafile)
            bias_image = scan_analyzer.getGrayscaleImage(
                                              'Infinity3_Bias_Image', datafile)
            img = scan_analyzer.processImage(data_image, bias_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                  "Infinity3_Laser_Beam_Processed_Image_atBkgndLoc", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- RAMAN LASER ZERO ORDER ---
            data_image = np.array(datafile['Raman_Laser_0Order_int'])
            bias_image = np.array(datafile['Raman_Bias_0Order_int'])
            background_image = np.array(
                                 datafile['Raman_Laser_0Order_atBkgndLoc_int'])
            img = scan_analyzer.processAndor(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                               "Raman_Laser_0Order_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- RAMAN LASER SPECTRUM ---
            data_image = np.array(datafile['Raman_Laser_Spectrum_int'])
            bias_image = np.array(datafile['Raman_Bias_Spectrum_int'])
            background_image = np.array(
                               datafile['Raman_Laser_Spectrum_atBkgndLoc_int'])
            img = scan_analyzer.processAndor(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                             "Raman_Laser_Spectrum_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- RAMAN WHITE LIGHT ZERO ORDER ---
            data_image = np.array(datafile['Raman_White_Light_0Order_int'])
            bias_image = np.array(datafile['Raman_Bias_0Order_int'])
            background_image = np.array(
                                datafile['Raman_White_Light_Bkgnd_0Order_int'])
            img = scan_analyzer.processAndor(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                         "Raman_White_Light_0Order_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # --- RAMAN WHITE LIGHT SPECTRUM ---
            data_image = np.array(datafile['Raman_White_Light_Spectrum_int'])
            bias_image = np.array(datafile['Raman_Bias_Spectrum_int'])
            background_image = np.array(
                              datafile['Raman_White_Light_Bkgnd_Spectrum_int'])
            img = scan_analyzer.processAndor(
                                      data_image, bias_image, background_image)
            zimg = scan_analyzer.reZeroImage(img)
            pimg = pdata.create_dataset(
                       "Raman_White_Light_Spectrum_Processed_Image", data=zimg)
            pimg.attrs.create("timestamp", datetime.datetime.now().isoformat())

    # Plot
    # fig1 = plt.figure(figsize=(20, 20))
    # plt.imshow(img, cmap='gray')
    # plt.show()

    processedData.close()
    data.close()
    print("Finished running script!")
