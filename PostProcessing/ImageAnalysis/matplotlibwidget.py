# -*- coding: utf-8 -*-
"""
Created on Wed May 17 13:12:48 2017

@author: mjh250
"""

import matplotlib

matplotlib.use("Qt5Agg")
matplotlib.rcParams['backend.qt5']='PyQt5'

from PyQt5.QtWidgets import QVBoxLayout, QSizePolicy, QWidget
from matplotlib.backends.backend_qt5agg import FigureCanvasQTAgg as FigureCanvas
from matplotlib.backends.backend_qt5agg import NavigationToolbar2QT as NavigationToolbar

from matplotlib.figure import Figure

class MatplotlibCanvas(FigureCanvas):
    def __init__(self):
        self.figure = Figure()
        self.axes = self.figure.add_subplot(111)
        self.axes.hold(False)
        FigureCanvas.__init__(self, self.figure)
        FigureCanvas.setSizePolicy(self,
        QSizePolicy.Expanding,
        QSizePolicy.Expanding)
        FigureCanvas.updateGeometry(self)
        
class MatplotlibWidget(QWidget):
     def __init__(self, parent = None):
        QWidget.__init__(self, parent)
        self.canvas = MatplotlibCanvas()
        self.mpl_toolbar = NavigationToolbar(self.canvas, self)
        self.vbl = QVBoxLayout()
        self.vbl.addWidget(self.canvas)
        self.vbl.addWidget(self.mpl_toolbar)
        self.setLayout(self.vbl)