# -*- coding: utf-8 -*-
"""
Created on Wed Nov 01 16:44:31 2017

@author: Matthew Horton (mjh250)
"""
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.stage.camera_stage_mapper_qt import CameraStageMapper
from nplab.instrument.stage.prior import ProScan
from nplab.instrument.shutter.thorlabs_sc10 import ThorLabsSC10
from nplab.instrument.shutter.BX51_uniblitz import Uniblitz

# Set up camera with click to move stage control
cam = LumeneraCamera(1)
cam.show_gui(blocking=False)
stage = ProScan("COM9")
mapper = CameraStageMapper(cam, stage)
mapper.show_gui(blocking=False)

# Display laser shutter control
shutter = ThorLabsSC10('COM1')
shutter.show_gui(blocking=False)

# Display white light shutter control
whiteShutter = Uniblitz("COM7")
whiteShutter.show_gui(blocking=False)
