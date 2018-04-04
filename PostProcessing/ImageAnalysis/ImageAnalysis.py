# -*- coding: utf-8 -*-
"""
Created on Sat Apr 22 2017

@author: Matthew Horton (mjh250)
"""

import sys
import window
import numpy as np
import matplotlib.pyplot as plt

from PyQt5.QtWidgets import QMainWindow
from PyQt5.QtWidgets import QApplication
from PyQt5.QtWidgets import QMessageBox
from PyQt5.QtWidgets import QFileDialog


class MainDialog(QMainWindow, window.Ui_MainWindow):
    def __init__(self, parent=None):
        super(MainDialog, self).__init__(parent)
        self.setupUi(self)

        self.btnRun.clicked.connect(self.btnRun_clicked)
        self.btnClose.clicked.connect(self.btnClose_clicked)
        self.btnBackgroundSelect.clicked.connect(
                                    self.btnBackgroundSelect_clicked)
        self.btnImageSelect.clicked.connect(self.btnImageSelect_clicked)
        self.btnOutputSelect.clicked.connect(self.btnOutputSelect_clicked)
        self.btnCrossSectionSelect.clicked.connect(
                                    self.btnCrossSectionSelect_clicked)
        self.btnRunCrossSection.clicked.connect(
                                    self.btnRunCrossSection_clicked)
        self.btnCrossSectionOutputSelect.clicked.connect(
                                    self.btnCrossSectionOutputSelect_clicked)
        self.btnSaveCrossSection.clicked.connect(
                                    self.btnSaveCrossSection_clicked)
        self.btnCropImageFilepath.clicked.connect(
                                    self.btnCropImageFilepath_clicked)
        self.btnCropOutputFilepath.clicked.connect(
                                    self.btnCropOutputFilepath_clicked)
        self.btnCropRun.clicked.connect(self.btnCropRun_clicked)
        self.btnCropSave.clicked.connect(self.btnCropSave_clicked)

    def btnBackgroundSelect_clicked(self):
        filepath = QFileDialog.getOpenFileName()[0]
        if filepath:
            self.txtBackgroundFilepath.setText(filepath)
        else:
            return

    def btnImageSelect_clicked(self):
        filepath = QFileDialog.getOpenFileName()[0]
        if filepath:
            self.txtImageFilepath.setText(filepath)
        else:
            return

    def btnOutputSelect_clicked(self):
        filepath = QFileDialog.getSaveFileName()[0]
        if filepath:
            self.txtOutputFilepath.setText(filepath)
        else:
            return

    def btnCrossSectionSelect_clicked(self):
        filepath = QFileDialog.getOpenFileName()[0]
        if filepath:
            self.txtCrossSectionFilepath.setText(filepath)
        else:
            return

    def btnCrossSectionOutputSelect_clicked(self):
        filepath = QFileDialog.getSaveFileName()[0]
        if filepath:
            self.txtCrossSectionOutputFilepath.setText(filepath)
        else:
            return

    def btnRun_clicked(self):
        outputPath = self.txtOutputFilepath.text()
        imagePath = self.txtImageFilepath.text()
        backgroundPath = self.txtBackgroundFilepath.text()
        # Import image and background image
        try:
            bg = np.genfromtxt(backgroundPath, skip_footer=30)
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error loading image from background filepath.',
                        QMessageBox.Ok, self).exec_()
            return
        try:
            my_data = np.genfromtxt(imagePath, skip_footer=30)
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error loading image from image filepath.',
                        QMessageBox.Ok, self).exec_()
            return
        # Process image, plot, and save.
        try:
            # Subtract background from image
            ndata = my_data-bg
            # Plot background corrected image
            fig1 = plt.figure(figsize=(20, 20))
            plt.imshow(np.transpose(ndata), cmap='gray')
            plt.show()
            # Save the background corrected image
            if outputPath:
                fig1.savefig(outputPath)
            else:
                QMessageBox(QMessageBox.Warning,
                            'Warning',
                            '''Output saved without user specified output
                            filepath.''',
                            QMessageBox.Ok, self).exec_()
                fig1.savefig(outputPath)
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error processing images, plotting, or saving.',
                        QMessageBox.Ok, self).exec_()
            return

    def btnRunCrossSection_clicked(self):
        crossSectionPath = self.txtCrossSectionFilepath.text()
        boolAxisIsX = self.radioX.isChecked()
        boolAxisIsY = self.radioY.isChecked()
        pos = self.spinboxCutPosition.value()
        # Import cross section image
        try:
            imCS = np.genfromtxt(crossSectionPath, delimiter=',')
            # imCS = np.genfromtxt(crossSectionPath, skip_footer=30)
            # imCS = np.delete(imCS, 0, 1)  # Throw away first column
        except Exception, e:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error loading file from image filepath. ' + repr(e),
                        QMessageBox.Ok, self).exec_()
            return
        # Process image and plot.
        try:
            if boolAxisIsX:
                cs = imCS[pos, :]
            elif boolAxisIsY:
                cs = imCS[:, pos]
            self.widgetDisplay.canvas.axes.plot(cs)
            self.widgetDisplay.canvas.draw()
            self.widgetDisplay.canvas.show()
        except Exception, e:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error taking cross section or plotting. ' + str(e),
                        QMessageBox.Ok, self).exec_()
            return
        return cs

    def btnSaveCrossSection_clicked(self):
        crossSectionOutputPath = self.txtCrossSectionOutputFilepath.text()
        print(crossSectionOutputPath)
        cs = self.btnRunCrossSection_clicked()
        print(cs)
        try:
            if crossSectionOutputPath:
                np.savetxt(crossSectionOutputPath, cs, delimiter=',')
            else:
                QMessageBox(QMessageBox.Warning,
                            'Warning',
                            '''Output saved without user specified output
                            filepath.''',
                            QMessageBox.Ok, self).exec_()
                np.savetxt(crossSectionOutputPath, cs, delimiter=',')
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error saving cross section.',
                        QMessageBox.Ok, self).exec_()
            return

# _______________________ UNDER CONSTRUCTION

    def btnCropImageFilepath_clicked(self):
        filepath = QFileDialog.getOpenFileName()[0]
        if filepath:
            self.txtCropImageFilepath.setText(filepath)
        else:
            return

    def btnCropOutputFilepath_clicked(self):
        filepath = QFileDialog.getSaveFileName()[0]
        if filepath:
            self.txtCropOutputFilepath.setText(filepath)
        else:
            return

    def btnCropRun_clicked(self):
        cropPath = self.txtCropImageFilepath.text()
        leftCrop = self.spinboxLeftMargin.value()
        rightCrop = self.spinboxRightMargin.value()
        topCrop = self.spinboxTopMargin.value()
        bottomCrop = self.spinboxBottomMargin.value()
        # Import cross section image
        try:
            imToCrop = np.genfromtxt(cropPath, skip_footer=30)
            imToCrop = np.delete(imToCrop, 0, 1)  # Throw away first column
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error loading file from image filepath.',
                        QMessageBox.Ok, self).exec_()
            return
        # Process image and plot.
        try:
            y, x = imToCrop.shape
            imToCrop = np.rot90(imToCrop)
            imCropped = imToCrop[topCrop:(x-bottomCrop),
                                 leftCrop:(y-rightCrop)]
            self.widgetDisplay.canvas.axes.imshow(imCropped, 'gray')
            self.widgetDisplay.canvas.draw()
            self.widgetDisplay.canvas.show()
        except Exception, e:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error cropping or plotting. ' + str(e),
                        QMessageBox.Ok, self).exec_()
            return
        return imCropped

    def btnCropSave_clicked(self):
        cropOutputPath = self.txtCropOutputFilepath.text()
        print(cropOutputPath)
        imCropped = self.btnCropRun_clicked()
        print(imCropped)
        try:
            if cropOutputPath:
                np.savetxt(cropOutputPath, imCropped, delimiter=',')
            else:
                QMessageBox(QMessageBox.Warning,
                            'Warning',
                            '''Output saved without user specified output
                            filepath.''',
                            QMessageBox.Ok, self).exec_()
                np.savetxt(cropOutputPath, imCropped, delimiter=',')
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error saving cross section.',
                        QMessageBox.Ok, self).exec_()
            return

# ________________________ /UNDER CONSTRUCTION

    def btnClose_clicked(self):
        sys.exit()

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
