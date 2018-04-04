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
            # del particle_number_list[1][-1] # Prevent last particle crash

            # --- Iterate through particle scans and process files
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
                    print('HERE WE GO!')
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
                    background_image = np.array(
                            datafile['Raman_Laser_0Order_atBkgndLoc_int'])
                    img = scan_analyzer.processAndor(
                            data_image, bias_image, background_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                            "Raman_Laser_0Order_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    imgMin = zimg.min()
                    imgMax = zimg.max()
                    zimgCrop = (zimg[65:(data_image.shape[0]-68),
                                             780:(data_image.shape[1]-770)])
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
                    pimg.attrs.create("minimum", imgMin)
                    pimg.attrs.create("maximum", imgMax)
                    pimg.attrs.create("integral", imgThumbIntegral)
                    pimg.attrs.create("average of 10 maxima", imgAvgMax)
        
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
                    pimg = pdata.create_dataset(
                        "Raman_White_Light_0Order_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())
                    imgMin = zimg.min()
                    imgMax = zimg.max()
                    zimgCrop = (zimg[65:(data_image.shape[0]-68),
                                             780:(data_image.shape[1]-770)])
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
        
        # =============================================================================
        # --- FIRST WHITE LIGHT IMAGE ---
                    for i in range(0, 14):
                        print('Time slot: ' + str(i))
                        data_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_FirstWhiteLight_Image_'+str(i), datafile)
                        bias_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_Bias_Image_'+str(i), datafile)
                        background_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_FirstBkgndWhiteLight_Image_'+str(i), datafile)
                        img = scan_analyzer.processImage(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                                "Infinity3_First_Processed_Image_"+str(i), data=zimg)
                        pimg.attrs.create(
                                "timestamp", datetime.datetime.now().isoformat())
        
                        # --- SECOND WHITE LIGHT IMAGE ---
                        data_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_SecondWhiteLight_Image_'+str(i), datafile)
                        bias_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_Bias_Image_'+str(i), datafile)
                        background_image = scan_analyzer.getGrayscaleImage(
                                'Infinity3_SecondWhiteLight_atBkgndLoc_Image_'+str(i),
                                datafile)
                        img = scan_analyzer.processImage(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                                "Infinity3_Second_Processed_Image_"+str(i), data=zimg)
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
                        data_image = np.array(datafile['Raman_Laser_0Order_int_'+str(i)])
                        bias_image = np.array(datafile['Raman_Bias_0Order_int_'+str(i)])
                        background_image = np.array(
                                datafile['Raman_Laser_0Order_atBkgndLoc_int_'+str(i)])
                        img = scan_analyzer.processAndor(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                                "Raman_Laser_0Order_Processed_Image_"+str(i), data=zimg)
                        pimg.attrs.create(
                                "timestamp", datetime.datetime.now().isoformat())
                        imgMin = zimg.min()
                        imgMax = zimg.max()
                        zimgCrop = (zimg[65:(data_image.shape[0]-68),
                                                 780:(data_image.shape[1]-770)])
                        imgThumbIntegral = zimgCrop.sum()
                        # Find 10 max values
                        imgThumbMax = []
                        zimgCropMaxRemoved = [item for sublist in zimgCrop for item in sublist]
                        for j in range(0, 10):
                            index, value = max(enumerate(zimgCropMaxRemoved),
                                               key=operator.itemgetter(1))
                            imgThumbMax.append(value)
                            del zimgCropMaxRemoved[index]
                        imgAvgMax = np.mean(imgThumbMax)
                        # Save calculated values
                        pimg.attrs.create("minimum", imgMin)
                        pimg.attrs.create("maximum", imgMax)
                        pimg.attrs.create("integral", imgThumbIntegral)
                        pimg.attrs.create("average of 10 maxima", imgAvgMax)
        
                        # --- RAMAN LASER SPECTRUM ---
                        data_image = np.array(datafile['Raman_Laser_Spectrum_int_'+str(i)])
                        bias_image = np.array(datafile['Raman_Bias_Spectrum_int_'+str(i)])
                        background_image = np.array(
                                datafile['Raman_Laser_Spectrum_atBkgndLoc_int_'+str(i)])
                        img = scan_analyzer.processAndor(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                                "Raman_Laser_Spectrum_Processed_Image_"+str(i), data=zimg)
                        pimg.attrs.create(
                                "timestamp", datetime.datetime.now().isoformat())
                        pimg.attrs.create(
                                "wavelengths", datafile['Raman_Laser_Spectrum_wl_'+str(i)])
        
                        # --- RAMAN WHITE LIGHT ZERO ORDER ---
                        data_image = np.array(
                                datafile['Raman_White_Light_0Order_int_'+str(i)])
                        bias_image = np.array(datafile['Raman_Bias_0Order_int_'+str(i)])
                        background_image = np.array(
                                datafile['Raman_White_Light_Bkgnd_0Order_int_'+str(i)])
                        img = scan_analyzer.processAndor(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                            "Raman_White_Light_0Order_Processed_Image_"+str(i), data=zimg)
                        pimg.attrs.create(
                                "timestamp", datetime.datetime.now().isoformat())
                        imgMin = zimg.min()
                        imgMax = zimg.max()
                        zimgCrop = (zimg[65:(data_image.shape[0]-68),
                                                 780:(data_image.shape[1]-770)])
                        imgThumbIntegral = zimgCrop.sum()
                        # Find 10 max values
                        imgThumbMax = []
                        zimgCropMaxRemoved = [item for sublist in zimgCrop for item in sublist]
                        for j in range(0, 10):
                            index, value = max(enumerate(zimgCropMaxRemoved),
                                               key=operator.itemgetter(1))
                            imgThumbMax.append(value)
                            del zimgCropMaxRemoved[index]
                        imgAvgMax = np.mean(imgThumbMax)
                        # Save calculated values
                        pimg.attrs.create("minimum", imgMin)
                        pimg.attrs.create("maximum", imgMax)
                        pimg.attrs.create("integral", imgThumbIntegral)
                        pimg.attrs.create("average of 10 maxima", imgAvgMax)
        
                        # --- RAMAN WHITE LIGHT SPECTRUM ---
                        data_image = np.array(
                                datafile['Raman_White_Light_Spectrum_int_'+str(i)])
                        bias_image = np.array(
                                datafile['Raman_Bias_Spectrum_int_'+str(i)])
                        background_image = np.array(
                                datafile['Raman_White_Light_Bkgnd_Spectrum_int_'+str(i)])
                        img = scan_analyzer.processAndor(
                                data_image, bias_image, background_image)
                        zimg = scan_analyzer.reZeroImage(img)
                        pimg = pdata.create_dataset(
                                "Raman_White_Light_Spectrum_Processed_Image_"+str(i),
                                data=zimg)
                        pimg.attrs.create(
                                "timestamp", datetime.datetime.now().isoformat())
                        pimg.attrs.create(
                          "wavelengths", datafile['Raman_White_Light_Spectrum_wl_'+str(i)])

            processedData.close()
            data.close()
            print("Finished processing " + filepath.split('/')[-1] + ' !')

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
            except OSError:
                shutil.os.mkdir(processed_filepath)

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
                    fig.savefig(processed_filepath +
                                '\scan' + str(scan_number) +
                                'particle' + str(particle_number) + 'time!!' +
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
                       cmap='gray')
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath +
                                '\scan' + str(scan_number) +
                                'particle' + str(particle_number) + 'time!!' +
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
                       cmap='gray')
                    plt.axis('equal')
                    ax.set_adjustable('box-forced')
                    fig.savefig(processed_filepath +
                                '\scan' + str(scan_number) +
                                'particle' + str(particle_number) + 'time!!' +
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
                    fig.savefig(processed_filepath +
                                '\scan' + str(scan_number) +
                                'particle' + str(particle_number) + 'time!!' +
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
                    fig.savefig(processed_filepath +
                                '\scan' + str(scan_number) +
                                'particle' + str(particle_number) + 'time!!' +
                                '_DFSpectrum.png',
                                bbox_inches='tight')
                    plt.close()
                    
                    for i in range(0, 14):
                        # --- FIRST WHITE LIGHT IMAGE ---
                        plt.ioff()  # No interactive plots. To suppress plot window
                        data_image = np.array(
                                datafile['Infinity3_First_Processed_Image_'+str(i)])
                        fig = plt.figure()
                        ax = plt.subplot(111)
                        ax.imshow(
                           data_image,
                           aspect='auto',
                           extent=[0, data_image.shape[1], 0, data_image.shape[0]],
                           cmap='gray')
                        plt.axis('equal')
                        ax.set_adjustable('box-forced')
                        fig.savefig(processed_filepath +
                                    '\scan' + str(scan_number) +
                                    'particle' + str(particle_number) + 'time' + str(i).zfill(2) +
                                    '_Inf3.png',
                                    bbox_inches='tight')
                        plt.close()
    
                        # --- RAMAN LASER 0 ORDER IMAGE ---
                        data_image = np.array(
                                datafile['Raman_Laser_0Order_Processed_Image_'+str(i)])
                        fig = plt.figure()
                        ax = plt.subplot(111)
                        ax.imshow(
                           data_image[:, 700:(data_image.shape[1]-700)],
                           aspect='auto',
                           extent=[700, data_image.shape[1]-700,
                                   0, data_image.shape[0]],
                           cmap='gray')
                        plt.axis('equal')
                        ax.set_adjustable('box-forced')
                        fig.savefig(processed_filepath +
                                    '\scan' + str(scan_number) +
                                    'particle' + str(particle_number) + 'time' + str(i).zfill(2) +
                                    '_Raman0Order.png',
                                    bbox_inches='tight')
                        plt.close()
    
                        # --- RAMAN DF 0 ORDER IMAGE ---
                        data_image = np.array(
                            datafile['Raman_White_Light_0Order_Processed_Image_'+str(i)])
                        fig = plt.figure()
                        ax = plt.subplot(111)
                        ax.imshow(
                           data_image[:, 700:(data_image.shape[1]-700)],
                           aspect='auto',
                           extent=[700, data_image.shape[1]-700,
                                   0, data_image.shape[0]],
                           cmap='gray')
                        plt.axis('equal')
                        ax.set_adjustable('box-forced')
                        fig.savefig(processed_filepath +
                                    '\scan' + str(scan_number) +
                                    'particle' + str(particle_number) + 'time' + str(i).zfill(2) +
                                    '_DF0Order.png',
                                    bbox_inches='tight')
                        plt.close()
    
                        # --- RAMAN LASER SPECTRUM IMAGE ---
                        dset = datafile['Raman_Laser_Spectrum_Processed_Image_'+str(i)]
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
                        fig.savefig(processed_filepath +
                                    '\scan' + str(scan_number) +
                                    'particle' + str(particle_number) + 'time' + str(i).zfill(2) +
                                    '_RamanSpectrum.png',
                                    bbox_inches='tight')
                        plt.close()
                        
                         # --- DARK FIELD SPECTRUM IMAGE ---
                        dset = datafile['Raman_White_Light_Spectrum_Processed_Image_'+str(i)]
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
                        fig.savefig(processed_filepath +
                                    '\scan' + str(scan_number) +
                                    'particle' + str(particle_number) + 'time' + str(i).zfill(2) +
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
        reZeroedImage = image - image.min()
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
