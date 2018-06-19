# -*- coding: utf-8 -*-
"""
Created on Tue Jan 23 15:38:07 2018

@author: Matthew
"""

from nplab.utils.gui import QtWidgets
from nplab.instrument.camera.Andor import Andor
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.camera.camera_with_location import CameraWithLocation
from nplab.instrument.spectrometer.shamdor import Shamdor
from nplab.instrument.stage.prior import ProScan
from ctypes import *


# Set up camera with click to move stage control
cam = LumeneraCamera(1)
#cam.show_gui(blocking=False)
stage = ProScan("COM9")
CWL = CameraWithLocation(cam, stage)
CWL.show_gui(blocking=False)

shamdor = Shamdor()
shamdor.CoolerON()
shamdor.SetTemperature = -90
shamdor.HSSpeed = 2
shamdor.dll.SetImageFlip(c_int(1), c_int(0))
shamdor.shamrock.show_gui(blocking=False)
#shamdor.show_gui(blocking=False)

# Andor
andor = Andor()  # wvl_to_pxl=32.5 / 1600, magnification=30, pxl_size=16)
app = QtWidgets.QApplication([])
ui1 = andor.get_control_widget()
ui2 = andor.get_preview_widget()
print ui1, ui2

ui1.show()
ui2.show()
