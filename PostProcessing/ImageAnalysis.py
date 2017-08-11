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
            imCS = np.genfromtxt(crossSectionPath, skip_footer=30)
            imCS = np.delete(imCS, 0, 1)  # Throw away first column
        except:
            QMessageBox(QMessageBox.Warning,
                        'An error has occurred.',
                        'Error loading file from image filepath.',
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
