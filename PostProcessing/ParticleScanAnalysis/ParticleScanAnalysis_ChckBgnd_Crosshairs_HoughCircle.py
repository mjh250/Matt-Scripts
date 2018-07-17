# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:33:43 2017

@author: mjh250
""" 

from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog
import matplotlib.pyplot as plt
from matplotlib.ticker import MultipleLocator

import h5py
import numpy as np
import datetime
import sys
import window
import shutil
import time
import operator
from scipy import ndimage
import cv2


sys.path.insert(0,"C:/Users/mjh250/Documents/Local mjh250/mjh250/Matt Scripts/PostProcessing/From Ilya")
from ShapeAnalysis import ShapeAnalysis


class MainDialog(QMainWindow, window.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self)

        self.btnRun.clicked.connect(self.btnRun_clicked)
        self.btnScanDataSelect.clicked.connect(
                                    self.btnScanDataSelect_clicked)
        self.btnOutputSelect.clicked.connect(self.btnOutputSelect_clicked)

        self.btnRunSWA.clicked.connect(self.btnRunSWA_clicked)
        self.btnOutputSelectSWA.clicked.connect(
                                    self.btnOutputSelectSWA_clicked)

        self.btnRunRadialProfile.clicked.connect(
                                    self.btnRunRadialProfile_clicked)
        self.btnOutputSelectRadial.clicked.connect(
                                    self.btnOutputSelectRadial_clicked)

        self.btnClose.clicked.connect(self.btnClose_clicked)

    def btnScanDataSelect_clicked(self):
        filepaths = QFileDialog.getOpenFileNames()[0]
        if filepaths:
            self.txtScanDataFilepath.setText(', '.join(filepaths))
        else:
            return

    def btnOutputSelect_clicked(self):
        filepath = QFileDialog.getSaveFileName(
                self,
                "Save File",
                "",
                "Data (*.hdf5 *.h5)")[0]
        if filepath:
            self.txtOutputFilepath.setText(filepath)
        else:
            return

    def btnRun_clicked(self):
        print("Initializing...")
        badBgndThreshold = 1.493e+06
        yBotCrop = 65
        yTopCrop = 68
        xLeftCrop = 780
        xRightCrop = 770
        filepaths = self.txtScanDataFilepath.text().split(', ')
        processed_filepath = self.txtOutputFilepath.text()
        if len(filepaths) > 1:
            processed_filepaths = [processed_filepath.split('.')[0] + str(i) +
                                   '.' + processed_filepath.split('.')[1]
                                   for i in range(len(filepaths))]
        else:
            processed_filepaths = [processed_filepath]

        for fileNum in range(len(filepaths)):
            filepath = filepaths[fileNum]
            processed_filepath = processed_filepaths[fileNum]
            try:
                shutil.os.remove(processed_filepath)
            except OSError:
                pass

            scan_analyzer = ParticleScanAnalysis()
            data = scan_analyzer.getScanDataFile(filepath)

            print("Preparing list of particles...")
            scan_number_list = []
            particle_number_list = []
            for sname in data:
                if sname.startswith('ParticleScannerScan_'):
                    scan_number_list.append(len(scan_number_list))
                    particle_num_sublist = []
                    for pname in data[sname]:
                        if pname.startswith('Particle_'):
                            particle_num_sublist.append(
                                                    len(particle_num_sublist))
                    particle_number_list.append(particle_num_sublist)
            
            if len(scan_number_list) == 0: # Handle case of old scan naming convention
                for sname in data['particleScans']:
                    if sname.startswith('scan'):
                        scan_number_list.append(len(scan_number_list))
                        particle_num_sublist = []
                        for pname in data['particleScans/'+sname]:
                            if pname.startswith('scan'):
                                particle_num_sublist.append(
                                                        len(particle_num_sublist))
                        particle_number_list.append(particle_num_sublist)
                # del particle_number_list[0][-1] # Prevent last particle crash
                                      
            for sublist in particle_number_list:
                sublist.sort()

            # --- Iterate through particle scans and process files
            old_background_image = None
            current_bgnd_particle_number = None
            rmn0_bad_bgnd_particle_list = []
            for scan_number in scan_number_list:
                for particle_number in particle_number_list[scan_number]:
                    # --- Print name of file being analyzed to track progress
                    print('Processing ParticleScannerScan_%s/Particle_%s'
                          % (scan_number, particle_number))

                    # --- GET SCAN DATA & INITIALIZE PROCESSED FILE ---
                    datafile = scan_analyzer.getScanDataSetGeneral(
                                            data, scan_number, particle_number)
                    processedData = h5py.File(processed_filepath, "a")
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
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())

                    # --- SECOND WHITE LIGHT IMAGE ---
                    data_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_SecondWhiteLight_Image', datafile)
                    bias_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_Bias_Image', datafile)
                    background_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_SecondWhiteLight_atBkgndLoc_Image',
                            datafile)
                    img = scan_analyzer.processImage(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                            "Infinity3_Second_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())

#                    # --- LASER BEAM IMAGE ---
#                    data_image = scan_analyzer.getGrayscaleImage(
#                            'Infinity3_Laser_Beam_Image', datafile)
#                    bias_image = scan_analyzer.getGrayscaleImage(
#                            'Infinity3_Bias_Image', datafile)
#                    img = scan_analyzer.processImage(data_image, bias_image)
#                    zimg = scan_analyzer.reZeroImage(img)
#                    pimg = pdata.create_dataset(
#                            "Infinity3_Laser_Beam_Processed_Image", data=zimg)
#                    pimg.attrs.create(
#                            "timestamp", datetime.datetime.now().isoformat())
#
#                    # --- LASER BEAM IMAGE AT BACKGROUND LOCATION ---
#                    data_image = scan_analyzer.getGrayscaleImage(
#                            'Infinity3_Laser_Beam_Image_atBkgndLoc', datafile)
#                    bias_image = scan_analyzer.getGrayscaleImage(
#                            'Infinity3_Bias_Image', datafile)
#                    img = scan_analyzer.processImage(data_image, bias_image)
#                    zimg = scan_analyzer.reZeroImage(img)
#                    pimg = pdata.create_dataset(
#                          "Infinity3_Laser_Beam_Processed_Image_atBkgndLoc",
#                          data=zimg)
#                    pimg.attrs.create(
#                            "timestamp", datetime.datetime.now().isoformat())

                    # --- RAMAN LASER ZERO ORDER ---
                    data_image = np.array(datafile['Raman_Laser_0Order_int'])
                    bias_image = np.array(datafile['Raman_Bias_0Order_int'])
                    
                    # Find a good background image
                    if "2017_12_14" in filepath:
                        tmp_datafile = scan_analyzer.getScanDataSetGeneral(
                                            data, scan_number, 1) # Prticle 1 has a good background
                        background_image = np.array(
                            tmp_datafile['Raman_Laser_0Order_atBkgndLoc_int'])
                        current_bgnd_particle_number = 1
                    elif "2018_02_16" in filepath:
                        tmp_datafile = scan_analyzer.getScanDataSetGeneral(
                                            data, scan_number, 32) # Prticle 32 has a good background
                        background_image = np.array(
                            tmp_datafile['Raman_Laser_0Order_atBkgndLoc_int'])
                        current_bgnd_particle_number = 32
                    else:
                        background_image = np.array(
                                datafile['Raman_Laser_0Order_atBkgndLoc_int'])
                        if old_background_image is None:
                            old_background_image = background_image
                            current_bgnd_particle_number = particle_number
                            
                        bgndCrop = (background_image[yBotCrop:(data_image.shape[0]-yTopCrop),
                                         xLeftCrop:(data_image.shape[1]-xRightCrop)])
                        oldBgndCrop = (old_background_image[yBotCrop:(data_image.shape[0]-yTopCrop),
                                         xLeftCrop:(data_image.shape[1]-xRightCrop)])
                        bgndSum = bgndCrop.sum()
                        oldbgndSum = oldBgndCrop.sum()
                        
                        if ((bgndSum > badBgndThreshold) and not (oldbgndSum > badBgndThreshold)):
                            print("Bad background image (Count sum = " + str(bgndSum) + "), using a previous one.")
                            background_image = old_background_image
                            rmn0_bad_bgnd_particle_list.append(particle_number)
                        elif (bgndSum > badBgndThreshold):
                            print("Bad background image (Count sum = " + str(bgndSum) + "), searching for a good one.")
                            found = False
                            for pn in range(1,len(particle_number_list[scan_number])-particle_number):
                                print("Looking forward by " + str(pn) + " particles for a good background.")
                                tmp_datafile = scan_analyzer.getScanDataSetGeneral(
                                                data, scan_number, particle_number+pn)
                                tmp_background_image = np.array(tmp_datafile['Raman_Laser_0Order_atBkgndLoc_int'])
                                tmpBgndCrop = (tmp_background_image[yBotCrop:(data_image.shape[0]-yTopCrop),
                                         xLeftCrop:(data_image.shape[1]-xRightCrop)])
                                if (tmpBgndCrop.sum() < badBgndThreshold):
                                    found = True
                                    background_image = tmp_background_image
                                    current_bgnd_particle_number = particle_number+pn
                                    break
                            rmn0_bad_bgnd_particle_list.append(particle_number)
                            if found is False:
                                print("WARNING: Could not find a suitable background image!")
                                
                        if oldbgndSum > badBgndThreshold:
                            old_background_image = background_image
                            current_bgnd_particle_number = particle_number
                        
                    # Process image
                    img = scan_analyzer.processAndor(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    imgMin = zimg.min()
                    imgMax = zimg.max()
                    zimgCrop = (zimg[yBotCrop:(data_image.shape[0]-yTopCrop),
                                     xLeftCrop:(data_image.shape[1]-xRightCrop)])
                    imgThumbIntegral = zimgCrop.sum()
                    # Find 10 max values
                    imgAvgMax = scan_analyzer.findMax10Average(zimgCrop)
                    # Save calculated values
                    pimg = pdata.create_dataset(
                            "Raman_Laser_0Order_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    pimg.attrs.create("minimum", imgMin)
                    pimg.attrs.create("maximum", imgMax)
                    pimg.attrs.create("integral", imgThumbIntegral)
                    pimg.attrs.create("average of 10 maxima", imgAvgMax)

                    # --- RAMAN LASER SPECTRUM ---
                    bgnd_datafile = scan_analyzer.getScanDataSetGeneral(
                                            data, scan_number, current_bgnd_particle_number)
                    data_image = np.array(datafile['Raman_Laser_Spectrum_int'])
                    bias_image = np.array(datafile['Raman_Bias_Spectrum_int'])
                    background_image = np.array(
                            bgnd_datafile['Raman_Laser_Spectrum_atBkgndLoc_int'])
                    img = scan_analyzer.processAndor(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                            "Raman_Laser_Spectrum_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    pimg.attrs.create(
                            "wavelengths", datafile['Raman_Laser_Spectrum_wl'])

                    # --- RAMAN WHITE LIGHT ZERO ORDER ---
                    data_image = np.array(
                            datafile['Raman_White_Light_0Order_int'])
                    bias_image = np.array(datafile['Raman_Bias_0Order_int'])
                    background_image = np.array(
                            datafile['Raman_White_Light_Bkgnd_0Order_int'])
                    img = scan_analyzer.processAndor(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    imgMin = zimg.min()
                    imgMax = zimg.max()
                    zimgCrop = (zimg[yBotCrop:(data_image.shape[0]-yTopCrop),
                                     xLeftCrop:(data_image.shape[1]-xRightCrop)])
                    imgThumbIntegral = zimgCrop.sum()
                    # Find 10 max values
                    imgThumbMax = []
                    zimgCropMaxRemoved = [item for sublist in zimgCrop for item in sublist]
                    for i in range(0, 10):
                        index, value = max(enumerate(zimgCropMaxRemoved),
                                           key=operator.itemgetter(1))
                        imgThumbMax.append(value)
                        del zimgCropMaxRemoved[index]
                    imgAvgMax = np.mean(imgThumbMax)
                    # Save calculated values
                    pimg = pdata.create_dataset(
                        "Raman_White_Light_0Order_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    pimg.attrs.create("minimum", imgMin)
                    pimg.attrs.create("maximum", imgMax)
                    pimg.attrs.create("integral", imgThumbIntegral)
                    pimg.attrs.create("average of 10 maxima", imgAvgMax)

                    # --- RAMAN WHITE LIGHT SPECTRUM ---
                    data_image = np.array(
                            datafile['Raman_White_Light_Spectrum_int'])
                    bias_image = np.array(
                            datafile['Raman_Bias_Spectrum_int'])
                    background_image = np.array(
                            datafile['Raman_White_Light_Bkgnd_Spectrum_int'])
                    img = scan_analyzer.processAndor(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                            "Raman_White_Light_Spectrum_Processed_Image",
                            data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    pimg.attrs.create(
                      "wavelengths", datafile['Raman_White_Light_Spectrum_wl'])

            processedData.close()
            data.close()
            print("Finished processing " + filepath.split('/')[-1] + ' !')
            print("The following particles were found to have nonzero "
                  "background in Raman 0 order, so a different background "
                  "image was used:")
            print(rmn0_bad_bgnd_particle_list)

    def btnOutputSelectSWA_clicked(self):
        filepath = QFileDialog.getExistingDirectory(
                self,
                "Choose output directory")
        if filepath:
            self.txtOutputFilepathSWA.setText(filepath)
        else:
            return

    def btnRunSWA_clicked(self):
        print("Initializing...")
        filepaths = self.txtScanDataFilepath.text().split(', ')
        processed_filepath = self.txtOutputFilepathSWA.text()
        if len(filepaths) > 1:
            processed_filepaths = [processed_filepath + str(i)
                                   for i in range(len(filepaths))]
        else:
            processed_filepaths = [processed_filepath]

        for fileNum in range(len(filepaths)):
            filepath = filepaths[fileNum]
            processed_filepath = processed_filepaths[fileNum]
            try:
                shutil.rmtree(processed_filepath)
                time.sleep(1)  # Wait for permissions to be freed up
                shutil.os.mkdir(processed_filepath)
                shutil.os.mkdir(processed_filepath+'/spots')
                shutil.os.mkdir(processed_filepath+'/rings')
                shutil.os.mkdir(processed_filepath+'/junks')
                shutil.os.mkdir(processed_filepath+'/asymmetrics')
                shutil.os.mkdir(processed_filepath+'/dims')
                shutil.os.mkdir(processed_filepath+'/UNCATEGORIZED')
            except OSError:
                shutil.os.mkdir(processed_filepath)
                shutil.os.mkdir(processed_filepath+'/spots')
                shutil.os.mkdir(processed_filepath+'/rings')
                shutil.os.mkdir(processed_filepath+'/junks')
                shutil.os.mkdir(processed_filepath+'/asymmetrics')
                shutil.os.mkdir(processed_filepath+'/dims')
                shutil.os.mkdir(processed_filepath+'/UNCATEGORIZED')

            scan_analyzer = ParticleScanAnalysis()
            data = scan_analyzer.getScanDataFile(filepath)

            print("Preparing list of particles...")
            scan_number_list = []
            particle_number_list = []
            for sname in data:
                if sname.startswith('scan'):
                    scan_number_list.append(int(sname[4:]))
                    particle_num_sublist = []
                    for pname in data[sname]:
                        if pname.startswith('particle'):
                            particle_num_sublist.append(int(pname[8:]))
                    particle_number_list.append(particle_num_sublist)
            
            for sublist in particle_number_list:
                sublist.sort()
            
            # Find range of counts to which images shall be scaled:
            pixel_max_list_Raman0Order = []
            pixel_max_list_White0Order = []
            print("Finding colourbar scaling parameters...")
            for n, scan_number in enumerate(scan_number_list):
                pixel_max_sublist_Raman0Order = []
                pixel_max_sublist_White0Order = []
                for particle_number in particle_number_list[n]:
                    datafile = scan_analyzer.getProcessedScanDataSet(
                                            data, scan_number, particle_number)
                    pixel_max_sublist_Raman0Order.append(
                            datafile["Raman_Laser_0Order_Processed_Image"].attrs['average of 10 maxima'])
                    pixel_max_sublist_White0Order.append(
                            datafile["Raman_White_Light_0Order_Processed_Image"].attrs['average of 10 maxima'])
                pixel_max_list_Raman0Order.append(pixel_max_sublist_Raman0Order)
                pixel_max_list_White0Order.append(pixel_max_sublist_White0Order)
                
            pixel_max_list_Raman0Order_r = list(pixel_max_list_Raman0Order)
            pixel_max_list_White0Order_r = list(pixel_max_list_White0Order)
            
            Raman0OrderMax = []
            White0OrderMax = []
            for n, scan_number in enumerate(scan_number_list):
                for i in range(0, min(10, len(pixel_max_list_Raman0Order_r[n]))):
                    index, value = max(enumerate(pixel_max_list_Raman0Order_r[n]),
                                       key=operator.itemgetter(1))
                    Raman0OrderMax.append(value)
                    del pixel_max_list_Raman0Order_r[n][index]
                    index, value = max(enumerate(pixel_max_list_White0Order_r[n]),
                                       key=operator.itemgetter(1))
                    White0OrderMax.append(value)
                    del pixel_max_list_White0Order_r[n][index]
            Raman0OrderAvgMax = np.mean(Raman0OrderMax)
            White0OrderAvgMax = np.mean(White0OrderMax)

            total_particle_count = 133
            rings = [22, 29, 35, 39, 82, 86, 90, 101, 117]
            dims = [19, 38, 53, 69, 71, 74, 99, 109, 118, 119]
            asymmetrics = [2, 3, 10, 13, 17, 33, 34, 36, 40, 45, 46, 49, 56, 61, 62, 66, 68, 73, 80, 84, 93, 96, 106, 108, 121] # Note: (3) might also be dim
            junks = [1, 5, 6, 20, 23, 25, 32, 47, 48, 52, 54, 55, 59, 60, 83, 89, 94, 102, 110, 112, 114, 123, 124, 127, 128, 129, 131] # Note: (54, 55, 60) interesting (quadruple NP), (83, 89, 131) double NP
            spots = [x for x in range(0, total_particle_count) if x not in rings+dims+asymmetrics+junks] # Note: 95 is interesting (small spot?)

            identities = []
            for i in range(0, total_particle_count):
                if i in spots:
                    identities.append('spots')
                elif i in junks:
                    identities.append('junks')
                elif i in rings:
                    identities.append('rings')
                elif i in asymmetrics:
                    identities.append('asymmetrics')
                elif i in dims:
                    identities.append('dims')
                else:
                    identities.append('UNCATEGORIZED')

            # OVERRIDE IF NECESSARY
            # scan_number_list = [3]

            # --- Iterate through particle scans and process files
            for n, scan_number in enumerate(scan_number_list):
                for particle_number in particle_number_list[n]:
                    # --- Print name of file being analyzed to track progress
                    print('Processing scan%s/particle%s' % (scan_number,
                                                            particle_number))

                    # --- GET SCAN DATA ---
                    datafile = scan_analyzer.getProcessedScanDataSet(
                                            data, scan_number, particle_number)

                    # --- FIRST WHITE LIGHT IMAGE ---
                    plt.ioff()  # No interactive plots. To suppress plot window
                    data_image = np.array(
                            datafile['Infinity3_First_Processed_Image'])
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image,
                       aspect='auto',
                       extent=[0, data_image.shape[1], 0, data_image.shape[0]],
                       cmap='gray')
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_Inf3.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- RAMAN DF 0 ORDER IMAGE ---
                    data_image = np.array(
                        datafile['Raman_White_Light_0Order_Processed_Image'],
                        dtype='float32')
                    # Stretch X-Axis by 137%
                    data_image_old_w = data_image.shape[1]
                    data_image = cv2.resize(data_image, None, fx = 1.37, fy = 1, interpolation = cv2.INTER_CUBIC)
                    data_image_new_w = data_image.shape[1]
                    data_image = data_image[:, int((data_image_new_w-data_image_old_w)/2):int(((data_image_new_w-data_image_old_w)/2)+data_image_old_w)]
                    
                    # Get Hough transform circles for a more accurate COM
                    data_image_8 = data_image[50:150, 750:850]
                    data_image_8 = np.round(255.0*((data_image_8-np.min(data_image_8))/float(np.max(data_image_8))))
                    data_image_8 = data_image_8.astype(np.uint8)
                    data_image_8_preThresholding = data_image_8
                    
                    # Find an ideal threshold, 
                    pixel_threshold = scan_analyzer.thresholdToTargetPixelCount(data_image_8, 80, 1000, 25)
                    # Globally threshold to the target number of bright pixels
                    ret, data_image_8 = cv2.threshold(data_image_8, pixel_threshold, 255, cv2.THRESH_TOZERO)
                    data_image_8 = cv2.medianBlur(data_image_8,9)
                    data_image_8 = cv2.equalizeHist(data_image_8)
                    
                    # Get a rough estimate of COM
                    ret, data_image_binary = cv2.threshold(data_image_8, 10, 255, cv2.THRESH_BINARY)
                    kernel = np.ones((11,11),np.uint8)
                    for i in range(0,3):
                        data_image_binary = cv2.morphologyEx(data_image_binary, cv2.MORPH_OPEN, kernel)
                    COM_y, COM_x = ndimage.measurements.center_of_mass(data_image_binary)
                    if (COM_y == COM_y) and (COM_x == COM_x):
                        COM_y = int(round(COM_y))
                        COM_x = int(round(COM_x))
                    else:
                        COM_y = int(round(data_image_8.shape[0]/2))
                        COM_x = int(round(data_image_8.shape[1]/2))
                    print("Rough COM: " + str([COM_x, COM_y]))
                    
                    # Find the outer circle of the doughnut by Hough circles method
                    big_circles = cv2.HoughCircles(image = data_image_8, # 8-bit, single channel image. If working with a color image, convert to grayscale first.
                           method = cv2.cv.CV_HOUGH_GRADIENT, # Defines the method to detect circles in images. Currently, the only implemented method is cv2.HOUGH_GRADIENT, which corresponds to the Yuen et al. paper.
                           dp = 1, # This parameter is the inverse ratio of the accumulator resolution to the image resolution (see Yuen et al. for more details). Essentially, the larger the dp gets, the smaller the accumulator array gets.
                           minDist = 140, # Minimum distance between the center (x, y) coordinates of detected circles. If the minDist is too small, multiple circles in the same neighborhood as the original may be (falsely) detected. If the minDist is too large, then some circles may not be detected at all.
                           param1 = 1, # Gradient value used to handle edge detection in the Yuen et al. method.
                           param2 = 2, # Accumulator threshold value for the cv2.HOUGH_GRADIENT method. The smaller the threshold is, the more circles will be detected (including false circles). The larger the threshold is, the more circles will potentially be returned.
                           minRadius = 14, # Minimum size of the radius (in pixels).
                           maxRadius = 17) # Maximum size of the radius (in pixels).
                    if big_circles is not None: 
                        big_circles = np.round(big_circles[0, :]).astype(np.uint8)
                        print("Hough method found big_circles: "+str(big_circles))
#                        for (x, y, r) in big_circles:
#                            # Draw the circle in the output image.
#                            cv2.circle(data_image, (x+750, y+50), r, (White0OrderAvgMax, White0OrderAvgMax, White0OrderAvgMax), thickness=1, lineType=8, shift=0)
#                            cv2.circle(data_image_8, (x, y), r, (128, 128, 128), thickness=1, lineType=8, shift=0)
                        if big_circles.shape[0] == 1:
                            COM2_x = big_circles[0][0]
                            COM2_y = big_circles[0][1]
                        else:
                            print("Could not find unique centre from Hough circles method.")
                            COM2_x = int(round(data_image_8.shape[1]/2))
                            COM2_y = int(round(data_image_8.shape[0]/2))
                    else:
                        print("Hough method found no big_circles.")
                        COM2_x = int(round(data_image_8.shape[1]/2))
                        COM2_y = int(round(data_image_8.shape[0]/2))

                    print("Distance between COM and Big Circle: " + str(np.linalg.norm(np.asarray([COM_x, COM_y]) - np.asarray([COM2_x, COM2_y]))))
                    
                    # Convolution method:
                    # First, generate the model 'donut' to look for (i.e. the kernel)
                    kernel = np.zeros((39, 39), dtype='uint8')
                    start_radius = 5
                    stop_radius = 15
                    for r in range(start_radius, stop_radius+1):
                        cv2.circle(kernel, (int(round(kernel.shape[1]/2)), int(round(kernel.shape[0]/2))), (stop_radius+start_radius)/2, 255, thickness=(stop_radius-start_radius), lineType=8, shift=0)
                    kernel = cv2.GaussianBlur(kernel, (7, 7), 2)
                    # Convolute kernel with image
                    image_convolution = cv2.filter2D(src=data_image_8_preThresholding,
                                                     dst = None,
                                                     ddepth=cv2.CV_64F,
                                                     kernel=kernel,
                                                     borderType=cv2.BORDER_REPLICATE)
                    (COM3_y, COM3_x) = np.unravel_index(image_convolution.argmax(), image_convolution.shape)

                    # Make plots and save figures
#                    fig = plt.figure()
#                    plt.imshow(image_convolution)
#                    plt.plot(COM3_x,COM3_y,'r.', mfc='none', mew=1, markersize=2)
#                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
#                                '/scan' + str(scan_number) +
#                                'particle' + str(particle_number) +
#                                '_Convolution.png',
#                                bbox_inches='tight')
#                    plt.close()
#                    
#                    fig = plt.figure()
#                    plt.imshow(np.pad(kernel, ((5, 6), (5, 6)), 'minimum'))
#                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
#                                '/scan' + str(scan_number) +
#                                'particle' + str(particle_number) +
#                                '_DonutKernel.png',
#                                bbox_inches='tight')
#                    plt.close()
#                    
#                    fig = plt.figure()
#                    plt.imshow(data_image_8)
#                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
#                                '/scan' + str(scan_number) +
#                                'particle' + str(particle_number) +
#                                '_DF0OrderForHough.png',
#                                bbox_inches='tight')
#                    plt.close()
#                    
#                    fig = plt.figure()
#                    plt.imshow(data_image_binary)
#                    plt.plot(COM_x,COM_y,'r.', mfc='none', mew=1, markersize=2)
#                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
#                                '/scan' + str(scan_number) +
#                                'particle' + str(particle_number) +
#                                '_DF0OrderForRoughCOM.png',
#                                bbox_inches='tight')
#                    plt.close()
                    
                    xstart = 775
                    xstop = data_image.shape[1]-775
                    ystart = 75
                    ystop = data_image.shape[0]-75
                    xmean = np.mean(np.asarray([xstart, xstop]))
                    ymean = np.mean(np.asarray([ystart, ystop]))
                    yoffset = int(data_image_8.shape[0]/2) - COM3_y
                    xoffset = -int(data_image_8.shape[1]/2) + COM3_x
                    print("Offset: "+str([xoffset,yoffset]))
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image[ystart-yoffset:ystop-yoffset, xstart+xoffset:xstop+xoffset],
                       aspect='auto',
                       extent=[xstart+xoffset, xstop+xoffset,
                               ystart+yoffset, ystop+yoffset],
                       cmap='gray',
                       vmin=0,
                       vmax=White0OrderAvgMax)
#                    plt.plot(COM_x+750,200-(COM_y+50),'b.', mfc='none', mew=2, markersize=10)
#                    plt.plot(COM2_x+750,200-(COM2_y+50),'g.', mfc='none', mew=2, markersize=10)
                    plt.plot(COM3_x+750,200-(COM3_y+50),'r.', mfc='none', mew=2, markersize=10)
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_DF0Order.png',
                                bbox_inches='tight')
                    plt.close()
                    
                    DF_data_image = data_image
                    
                    # debuggin plot
#                    fig = plt.figure()
#                    ax = plt.subplot(111)
#                    ax.imshow(
#                       threshImg,
#                       aspect='auto',
#                       extent=[0, data_image.shape[1], 0, data_image.shape[0]],
#                       cmap='gray',
#                       vmin=0,
#                       vmax=threshImg.max())
#                    plt.plot(COM_x,COM_y,'ro', mfc='none', mew=1, markersize=1)
#                    plt.plot(bbox[0][0],bbox[0][1],'r.', mfc='none', mew=1, markersize=1)
#                    plt.plot(bbox[1][0],bbox[1][1],'r.', mfc='none', mew=1, markersize=1)
#                    plt.plot(bbox[2][0],bbox[2][1],'r.', mfc='none', mew=1, markersize=1)
#                    plt.plot(bbox[3][0],bbox[3][1],'r.', mfc='none', mew=1, markersize=1)
#                    plt.axis('equal')
#                    ax.set_adjustable('box-forced')
#                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
#                                '/scan' + str(scan_number) +
#                                'particle' + str(particle_number) +
#                                '_DF0OrderThresh.png',
#                                bbox_inches='tight')
#                    plt.close()

                    # --- RAMAN LASER 0 ORDER IMAGE ---
                    data_image = np.array(
                            datafile['Raman_Laser_0Order_Processed_Image'],
                            dtype='float32')
                    # Stretch X-Axis by 137%
                    data_image_old_w = data_image.shape[1]
                    data_image = cv2.resize(data_image, None, fx = 1.37, fy = 1, interpolation = cv2.INTER_CUBIC)
                    data_image_new_w = data_image.shape[1]
                    data_image = data_image[:, int((data_image_new_w-data_image_old_w)/2):int(((data_image_new_w-data_image_old_w)/2)+data_image_old_w)]
                    
                    
                    xstart = 775
                    xstop = data_image.shape[1]-775
                    ystart = 75
                    ystop = data_image.shape[0]-75
                    xmean = np.mean(np.asarray([xstart, xstop]))
                    ymean = np.mean(np.asarray([ystart, ystop]))
                    yoffset = int(data_image_8.shape[0]/2) - COM3_y
                    xoffset = -int(data_image_8.shape[1]/2) + COM3_x
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image[ystart-yoffset:ystop-yoffset, xstart+xoffset:xstop+xoffset],
                       aspect='auto',
                       extent=[xstart+xoffset, xstop+xoffset,
                               ystart+yoffset, ystop+yoffset],
                       cmap='gray',
                       vmin=0,
                       vmax=Raman0OrderAvgMax)
                    weight = 3
                    plt.plot([xoffset+xmean,xoffset+xmean],[yoffset+ystart+15,yoffset+ystart+5],'r-', linestyle = "-", lw=weight)
                    plt.plot([xoffset+xmean,xoffset+xmean],[yoffset+ystart+45,yoffset+ystart+35],'r-', linestyle = "-", lw=weight)
                    plt.plot([xoffset+xstart+15,xoffset+xstart+5],[yoffset+ymean,yoffset+ymean],'r-', linestyle = "-", lw=weight)
                    plt.plot([xoffset+xstart+45,xoffset+xstart+35],[yoffset+ymean,yoffset+ymean],'r-', linestyle = "-", lw=weight)
                    plt.plot(xoffset+xmean,yoffset+ymean,'r.', mfc='none', mew=weight, markersize=1)
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_Raman0Order.png',
                                bbox_inches='tight')
                    plt.close()
                    
                    Raman_data_image = data_image
                    
                    # --- PLOT DF AND RAMAN SIDE-BY-SIDE
                    contrast_enhancement_factor = 0.75 # >0, bigger is darker
                    pixel_min = 0 # lower limit of the colorscale
                    f, (ax1, ax2) = plt.subplots(1, 2, sharex=False, sharey=False)
                    
                    ax1.imshow(
                       DF_data_image[ystart-yoffset:ystop-yoffset, xstart+xoffset:xstop+xoffset],
                       aspect='auto',
                       extent=[xstart+xoffset, xstop+xoffset,
                               ystart+yoffset, ystop+yoffset],
                       cmap='gray',
                       vmin=pixel_min,
                       vmax=White0OrderAvgMax*contrast_enhancement_factor)
                    ax1.plot(COM3_x+750,200-(COM3_y+50),'r.', mfc='none', mew=2, markersize=10)
                    ax1.axis('equal')
                    ax1.set_adjustable('box-forced')
                    ax1.set_title('DF')
                    
                    ax2.imshow(
                       Raman_data_image[ystart-yoffset:ystop-yoffset, xstart+xoffset:xstop+xoffset],
                       aspect='auto',
                       extent=[xstart+xoffset, xstop+xoffset,
                               ystart+yoffset, ystop+yoffset],
                       cmap='gray',
                       vmin=pixel_min,
                       vmax=Raman0OrderAvgMax*contrast_enhancement_factor)
                    weight = 3
                    ax2.plot([xoffset+xmean,xoffset+xmean],[yoffset+ystart+15,yoffset+ystart+5],'r-', linestyle = "-", lw=weight)
                    ax2.plot([xoffset+xmean,xoffset+xmean],[yoffset+ystart+45,yoffset+ystart+35],'r-', linestyle = "-", lw=weight)
                    ax2.plot([xoffset+xstart+15,xoffset+xstart+5],[yoffset+ymean,yoffset+ymean],'r-', linestyle = "-", lw=weight)
                    ax2.plot([xoffset+xstart+45,xoffset+xstart+35],[yoffset+ymean,yoffset+ymean],'r-', linestyle = "-", lw=weight)
                    ax2.plot(xoffset+xmean,yoffset+ymean,'r.', mfc='none', mew=weight, markersize=1)
                    ax2.axis('equal')
                    ax2.set_adjustable('box-forced')
                    ax2.set_title('Raman')
                    
                    f.suptitle('Particle ' + str(particle_number), y = 0.8)
                    
                    f.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_RamanAndDF.png',
                                bbox_inches='tight')
                    plt.close()

                    # Amount by which to crop spectra (from the left)
                    spectrum_left_crop = 300
                    
                    # --- RAMAN LASER SPECTRUM IMAGE ---
                    dset = datafile['Raman_Laser_Spectrum_Processed_Image']
                    data_image = np.array(dset)
                    data_image = data_image[:,spectrum_left_crop:]
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    aspect_factor = (dset.attrs['wavelengths'][-1] - dset.attrs['wavelengths'][spectrum_left_crop])/(data_image.shape[1] - 1)
                    ax.imshow(
                       data_image,
                       aspect=aspect_factor,
                       extent=[dset.attrs['wavelengths'][spectrum_left_crop], dset.attrs['wavelengths'][-1], 0, data_image.shape[0]],
                       cmap='gray')
                    #plt.axis('scaled')
                    #ax.set_adjustable('box-forced')
                    plt.minorticks_on()
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_RamanSpectrum.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- DARK FIELD SPECTRUM IMAGE ---
                    dset = datafile['Raman_White_Light_Spectrum_Processed_Image']
                    data_image = np.array(dset)
                    data_image = data_image[:,spectrum_left_crop:]
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    aspect_factor = (dset.attrs['wavelengths'][-1] - dset.attrs['wavelengths'][spectrum_left_crop])/(data_image.shape[1] - 1)
                    ax.imshow(
                       data_image,
                       aspect=aspect_factor,
                       extent=[dset.attrs['wavelengths'][spectrum_left_crop], dset.attrs['wavelengths'][-1], 0, data_image.shape[0]],
                       cmap='gray')
                    #plt.axis('equal')
                    #ax.set_adjustable('box-forced')
                    plt.minorticks_on()
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_DFSpectrum.png',
                                bbox_inches='tight')
                    plt.close()

            data.close()
            print("Finished processing " + filepath.split('/')[-1] + ' !')
        return

    def btnOutputSelectRadial_clicked(self):
        filepath = QFileDialog.getExistingDirectory(
                self,
                "Choose output directory")
        if filepath:
            self.txtOutputFilepathRadial.setText(filepath)
        else:
            return

    def btnRunRadialProfile_clicked(self):
        filepath = self.txtScanDataFilepath.text()
        processed_filepath = self.txtOutputFilepathRadial.text()
        centre = [self.spnRadialCentreX.value(), self.spnRadialCentreY.value()]
        scan_number = 0
        particle_number = 81

        scan_analyzer = ParticleScanAnalysis()
        data = scan_analyzer.getScanDataFile(filepath)

        print("Preparing list of particles...")
        scan_number_list = []
        particle_number_list = []
        for name in data:
            if name.startswith('scan'):
                scan_number_list.append(len(scan_number_list))
        for name in data['/scan0']:
            if name.startswith('scan_'):
                particle_number_list.append(len(particle_number_list))

        # --- Print name of file being analyzed to track progress
        print('Processing scan%s/particle%s' % (scan_number,
                                                particle_number))

        # --- GET SCAN DATA ---
        datafile = scan_analyzer.getScanDataSetGeneral(
                                data, scan_number, particle_number)

        # --- RAMAN LASER 0 ORDER IMAGE ---
        data_image = np.array(
                datafile['Raman_Laser_0Order_Processed_Image'])
        data_image_crop = data_image[50:(data_image.shape[0]-50),
                                     750:(data_image.shape[1]-750)]

        plt.figure()
        plt.imshow(data_image_crop)

        data_image_radial = scan_analyzer.processRadialProfile(
                                                    data_image_crop, centre)
        fig = plt.figure()
        ax = plt.subplot(111)
        ax.plot(data_image_radial)
        plt.axis('equal')
        plt.xlabel('Radial coordinate (in pixels)')
        plt.ylabel('Counts')
        ax.set_adjustable('box-forced')
        fig.savefig(processed_filepath +
                    '\Raman0Order_scan' + str(scan_number) +
                    'particle' + str(particle_number) + '_radial.png',
                    bbox_inches='tight')
        plt.close()

        data_image_azimuthal = scan_analyzer.processAzimuthalProfile(
                                                    data_image_crop, centre)
        fig = plt.figure()
        ax = plt.subplot(111, projection='polar')
        ax.plot(np.linspace(0, 2*np.pi, data_image_azimuthal.size),
                data_image_azimuthal)
        plt.axis('equal')
        ax.set_adjustable('box-forced')
        fig.savefig(processed_filepath +
                    '\Raman0Order_scan' + str(scan_number) +
                    'particle' + str(particle_number) + '_azimuthal.png',
                    bbox_inches='tight')
        plt.close()
        return

    def btnClose_clicked(self):
        sys.exit()


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
            # processed_image = ((data_image - bias_image - background_image) /
            #                     (white_reference_image - background_image))
        processed_image = ((data_image - background_image) /
                           (white_reference_image - background_image))
        return processed_image

    def processAndor(self, data_image, bias_image=None, background_image=None):
        if bias_image is None:
            bias_image = np.full_like(data_image, 0)
        if background_image is None:
            background_image = np.full_like(data_image, 0)
        # processed_image = data_image - bias_image - background_image
        processed_image = data_image - background_image
        return processed_image

    def getScanDataFile(self, filepath):
        print(filepath)
        datafile = h5py.File(filepath, 'a')
        return datafile

    def getScanDataSetGeneral(self, datafile, scan_number, particle_number):
        try:
            dataset = datafile['/particleScans/scan%s/scan_%s'
                               % (scan_number, particle_number)]
        except:
            try:
                dataset = datafile['/ParticleScannerScan_%s/Particle_%s'
                                   % (scan_number, particle_number)]
            except:
                dataset = datafile['/scan%s/particle%s'
                                   % (scan_number, particle_number)]

        return dataset

    def getScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/particleScans/scan%s/scan_%s'
                           % (scan_number, particle_number)]
        return dataset

    def getProcessedScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/scan%s/particle%s'
                           % (scan_number, particle_number)]
        return dataset

    def reZeroImage(self, image):
        # reZeroedImage = image - image.min()
        reZeroedImage = np.asarray(image).clip(min=0)
        return reZeroedImage

    def processRadialProfile(self, data, centre):
        y, x = np.indices((data.shape))
        r = np.sqrt((x - centre[0])**2 + (y - centre[1])**2)
        r = r.astype(np.int)

        plt.figure()
        plt.imshow(r)

        tbin = np.bincount(r.ravel(), data.ravel())
        nr = np.bincount(r.ravel())
        radialprofile = tbin / nr
        return radialprofile

    def processAzimuthalProfile(self, data, centre):
        y, x = np.indices((data.shape))
        r = np.arctan2((y - centre[1]), (x - centre[0]))
        r_positive = r + np.absolute(np.amin(r))
        r_positive_scaled = r_positive * 100
        r_positive_scaled = r_positive_scaled.astype(np.int)

        tbin = np.bincount(r_positive_scaled.ravel(), data.ravel())
        nr = np.bincount(r_positive_scaled.ravel())

        # azimuthalprofile = tbin / nr
        azimuthalprofile = np.divide(tbin, nr,
                                     out=np.zeros_like(tbin), where=(nr != 0))
        azimuthalprofileLEFT = np.roll(azimuthalprofile, 1)
        azimuthalprofileRIGHT = np.roll(azimuthalprofile, -1)
        azimuthalprofileAVG = (azimuthalprofileLEFT + azimuthalprofileRIGHT)/2
        azimuthalprofile[azimuthalprofile <= 0] = azimuthalprofileAVG[
                                                        azimuthalprofile <= 0]
        return azimuthalprofile

    def findMax10Average(self, data):
        imgThumbMax = []
        dataMaxRemoved = [item for sublist in data for item in sublist]
        for i in range(0, 10):
            index, value = max(enumerate(dataMaxRemoved),
                               key=operator.itemgetter(1))
            imgThumbMax.append(value)
            del dataMaxRemoved[index]
        imgAvgMax = np.mean(imgThumbMax)
        return imgAvgMax
    
    def thresholdToTargetPixelCount(self, image, pixel_threshold, target_n_bright_pixels, tolerence_n_bright_pixels):
        for i in range(1, 100):
            n_bright_pixels = len(image[np.where(image > pixel_threshold)])
            if n_bright_pixels < (target_n_bright_pixels - tolerence_n_bright_pixels):
                pixel_threshold = pixel_threshold - 1
            elif n_bright_pixels > (target_n_bright_pixels + tolerence_n_bright_pixels):
                pixel_threshold = pixel_threshold + 1
            else:
                break
            if i >= 99:
                pixel_threshold = 80
                print("Could not find a good pixel threshold.")
        return pixel_threshold

if __name__ == '__main__':
    app = QApplication(sys.argv)
    MainWindow = MainDialog()
    MainWindow.show()
    sys.exit(app.exec_())
    try:
        from IPython.lib.guisupport import start_event_loop_qt5
        start_event_loop_qt5(app)
    except ImportError:
        app.exec_()
