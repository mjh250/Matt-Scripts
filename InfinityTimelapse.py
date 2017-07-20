# -*- coding: utf-8 -*-
"""
Created on Thu Mar 02 13:54:02 2017

@author: Matthew Horton
"""
from nplab.instrument.camera.lucam import Lucam, API
from nplab.instrument.camera import CameraParameter
from time import sleep
import numpy

def __init__():
    cam = Lucam()
    cam.__init__()
    for pname in cam.PROPERTY.keys():
        setattr(cam, pname, CameraParameter(pname))
        
    return cam

def SingleShot(cam,filename='test.bmp'):
   # Set camera settings
    frameformat = cam.GetFormat()[0]
    snapshot = cam.Snapshot(exposure=cam.GetProperty('exposure')[0], gain=cam.GetProperty('gain')[0]/2,
                            timeout=cam.GetProperty('exposure')[0] + 500.0, format=frameformat)
    # Take and save data                          
    data = cam.TakeSnapshot(snapshot)
    colorData = cam.ConvertFrameToRgb24(frameformat, data.ctypes.data_as(API.pBYTE))
    numpyColorData = numpy.ctypeslib.as_array(colorData, shape = None)
    numpyColorDataT = numpy.rollaxis(numpyColorData,2)
    cam.SaveImage(numpyColorDataT,filename)
    
if __name__ == '__main__':
    cam = __init__()
    for x in range(0, 3):
        print('Taking image '+str(x)+ '...')
        SingleShot(cam, 'Timelapse'+str(x)+'.bmp')
        sleep(2)
        
    print('Done.')