# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 16:44:31 2017

@author: Matthew Horton (mjh250)
"""
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.camera.camera_with_location import CameraWithLocation
from nplab.instrument.stage.prior import ProScan
from nplab.instrument.shutter.thorlabs_sc10 import ThorLabsSC10
from nplab.instrument.shutter.BX51_uniblitz import Uniblitz

# Set up camera with click to move stage control
cam = LumeneraCamera(1)
# cam.show_gui(blocking=False)
stage = ProScan("COM4")
CWL = CameraWithLocation(cam, stage)
CWL.show_gui(blocking=False)

# Display laser shutter control
shutter = ThorLabsSC10('COM1')
shutter.show_gui(blocking=False)

# Display white light shutter control
whiteShutter = Uniblitz("COM10")
whiteShutter.show_gui(blocking=False)
