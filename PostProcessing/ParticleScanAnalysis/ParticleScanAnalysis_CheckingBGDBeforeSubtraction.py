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
        badBgndThreshold = 6.75e+06
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
#            for sname in data['particleScans']:
            for sname in data:
                if sname.startswith('ParticleScannerScan_'):
#                if sname.startswith('scan'):
                    scan_number_list.append(len(scan_number_list))
                    particle_num_sublist = []
                    for pname in data[sname]:
#                    for pname in data['particleScans/'+sname]:
                        if pname.startswith('Particle_'):
#                        if pname.startswith('scan'):
                            particle_num_sublist.append(
                                                    len(particle_num_sublist))
                    particle_number_list.append(particle_num_sublist)
            # del particle_number_list[1][-1] # Prevent last particle crash
            
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
                for i in range(0, 10):
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

            total_particle_count = 152
            rings = [39, 42, 52, 54, 58, 64, 72, 87, 90, 94, 102, 109, 111, 123,
                     128, 129]
            dims = [3, 9, 18, 19, 21, 29, 48, 71, 106, 113, 126, 142]
            asymmetrics = [22, 24, 70, 89, 100, 103, 114, 133, 140, 141, 145]
            junks = [1, 2, 4, 5, 7, 10, 16, 17, 20, 25, 27, 32, 33, 35, 38, 41,
                     45, 49, 50, 55, 56, 63, 66, 67, 73, 74, 75, 76, 77, 80, 81, 82,
                     83, 84, 85, 88, 95, 96, 98, 101, 104, 115, 117, 125, 130, 131, 132, 134,
                     135, 136, 137, 138, 139, 143, 144, 146, 147, 148, 149, 150, 151]
            spots = [x for x in range(0, total_particle_count) if x not in rings+dims+asymmetrics+junks]
            
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

                    # --- RAMAN LASER 0 ORDER IMAGE ---
                    data_image = np.array(
                            datafile['Raman_Laser_0Order_Processed_Image'])
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image[:, 700:(data_image.shape[1]-700)],
                       aspect='auto',
                       extent=[700, data_image.shape[1]-700,
                               0, data_image.shape[0]],
                       cmap='gray',
                       vmin=0,
                       vmax=Raman0OrderAvgMax)
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_Raman0Order.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- RAMAN DF 0 ORDER IMAGE ---
                    data_image = np.array(
                        datafile['Raman_White_Light_0Order_Processed_Image'])
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image[:, 700:(data_image.shape[1]-700)],
                       aspect='auto',
                       extent=[700, data_image.shape[1]-700,
                               0, data_image.shape[0]],
                       cmap='gray',
                       vmin=0,
                       vmax=White0OrderAvgMax)
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_DF0Order.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- RAMAN LASER SPECTRUM IMAGE ---
                    dset = datafile['Raman_Laser_Spectrum_Processed_Image']
                    data_image = np.array(dset)
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image,
                       aspect='auto',
                       extent=[0, data_image.shape[1], 0, data_image.shape[0]],
                       cmap='gray')
                    plt.axis('equal')
                    plt.xticks(np.linspace(0, data_image.shape[1], 8),
                               np.around(
                                   np.linspace(
                                       dset.attrs['wavelengths'][0],
                                       dset.attrs['wavelengths'][-1], 8),
                                   decimals=1))
                    ax.set_adjustable('box-forced')
                    plt.minorticks_on()
                    minorLocator = MultipleLocator(10)
                    ax.yaxis.set_minor_locator(minorLocator)
                    fig.savefig(processed_filepath + '/' + identities[particle_number] +
                                '/scan' + str(scan_number) +
                                'particle' + str(particle_number) +
                                '_RamanSpectrum.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- DARK FIELD SPECTRUM IMAGE ---
                    dset = datafile['Raman_White_Light_Spectrum_Processed_Image']
                    data_image = np.array(dset)
                    fig = plt.figure()
                    ax = plt.subplot(111)
                    ax.imshow(
                       data_image,
                       aspect='auto',
                       extent=[0, data_image.shape[1], 0, data_image.shape[0]],
                       cmap='gray')
                    plt.axis('equal')
                    plt.xticks(np.linspace(0, data_image.shape[1], 8),
                               np.around(
                                   np.linspace(
                                       dset.attrs['wavelengths'][0],
                                       dset.attrs['wavelengths'][-1], 8),
                                   decimals=1))
                    ax.set_adjustable('box-forced')
                    plt.minorticks_on()
                    minorLocator = MultipleLocator(10)
                    ax.yaxis.set_minor_locator(minorLocator)
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
