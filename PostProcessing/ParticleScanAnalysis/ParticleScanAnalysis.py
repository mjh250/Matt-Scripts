# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:33:43 2017

@author: mjh250
"""

from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QFileDialog
import matplotlib.pyplot as plt

import h5py
import numpy as np
import datetime
import sys
import window
import shutil
import time


class MainDialog(QMainWindow, window.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self)

        self.btnRun.clicked.connect(self.btnRun_clicked)
        self.btnScanDataSelect.clicked.connect(
                                    self.btnScanDataSelect_clicked)
        self.btnOutputSelect.clicked.connect(self.btnOutputSelect_clicked)

        self.btnRunSWA.clicked.connect(self.btnRunSWA_clicked)
        self.btnScanDataSelectSWA.clicked.connect(
                                    self.btnScanDataSelectSWA_clicked)
        self.btnOutputSelectSWA.clicked.connect(
                                    self.btnOutputSelectSWA_clicked)

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
                    print('Processing scan%s/scan_%s' % (scan_number,
                                                         particle_number))

                    # --- GET SCAN DATA & INITIALIZE PROCESSED FILE ---
                    datafile = scan_analyzer.getScanDataSet(data, scan_number,
                                                            particle_number)
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

                    # --- LASER BEAM IMAGE ---
                    data_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_Laser_Beam_Image', datafile)
                    bias_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_Bias_Image', datafile)
                    img = scan_analyzer.processImage(data_image, bias_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                            "Infinity3_Laser_Beam_Processed_Image", data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())

                    # --- LASER BEAM IMAGE AT BACKGROUND LOCATION ---
                    data_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_Laser_Beam_Image_atBkgndLoc', datafile)
                    bias_image = scan_analyzer.getGrayscaleImage(
                            'Infinity3_Bias_Image', datafile)
                    img = scan_analyzer.processImage(data_image, bias_image)
                    zimg = scan_analyzer.reZeroImage(img)
                    pimg = pdata.create_dataset(
                          "Infinity3_Laser_Beam_Processed_Image_atBkgndLoc",
                          data=zimg)
                    pimg.attrs.create(
                            "timestamp", datetime.datetime.now().isoformat())

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

    def btnScanDataSelectSWA_clicked(self):
        filepaths = QFileDialog.getOpenFileNames()[0]
        if filepaths:
            self.txtScanDataFilepathSWA.setText(', '.join(filepaths))
        else:
            return

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
        filepaths = self.txtScanDataFilepathSWA.text().split(', ')
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
            for name in data:
                if name.startswith('scan'):
                    scan_number_list.append(len(scan_number_list))
            for name in data['/scan0']:
                if name.startswith('particle'):
                    particle_number_list.append(len(particle_number_list))

            # --- Iterate through particle scans and process files
            for scan_number in scan_number_list:
                for particle_number in particle_number_list:
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
                                '\Inf3_scan' + str(scan_number) +
                                'particle' + str(particle_number) + '.png',
                                bbox_inches='tight')
                    plt.close()

                    # --- RAMAN LASER 0 ORDER IMAGE ---
                    data_image = np.array(
                            datafile['Raman_Laser_0Order_Processed_Image'])
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
                                '\Raman0Order_scan' + str(scan_number) +
                                'particle' + str(particle_number) + '.png',
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
                    fig.savefig(processed_filepath +
                                '\RamanSpectrum_scan' + str(scan_number) +
                                'particle' + str(particle_number) + '.png',
                                bbox_inches='tight')
                    plt.close()

            data.close()
            print("Finished processing " + filepath.split('/')[-1] + ' !')
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

    def getScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/particleScans/scan%s/scan_%s'
                           % (scan_number, particle_number)]
        return dataset

    def getProcessedScanDataSet(self, datafile, scan_number, particle_number):
        dataset = datafile['/scan%s/particle%s'
                           % (scan_number, particle_number)]
        return dataset

    def reZeroImage(self, image):
        reZeroedImage = image - min(image.flatten())
        return reZeroedImage

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
