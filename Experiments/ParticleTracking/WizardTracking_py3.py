# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 13:29:18 2018

@author: np-albali
"""
from __future__ import division
from __future__ import print_function
from builtins import str
from builtins import range
from past.utils import old_div
import sys
sys.path.insert(0, "c:/users/np-albali/documents/github/particle_tracking_app")

from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.camera.camera_with_location import CameraWithLocation
from nplab.instrument.spectrometer.spectrometer_aligner import SpectrometerAligner
from nplab.instrument.stage.prior import ProScan
from nplab.instrument.shutter.thorlabs_sc10 import ThorLabsSC10
from nplab.instrument.shutter.BX51_uniblitz import Uniblitz
from nplab.instrument.spectrometer.shamdor import Shamdor
from pyvcam import pvc
from pyvcam.camera import Camera
# import nplab.datafile as df
from particle_tracking_app.particle_tracking_wizard import TrackingWizard
import numpy as np
import datetime
import time


HSSpeed = 2
andor_exposure_time_0Order = 5
andor_exposure_time_spec = andor_exposure_time_0Order*3
donut_z_offset = 0  # z offset to move focus to the donut mode
centre_wavelength = 695
tile_edge_width_to_ignore = 250

# Set up camera with click to move stage control
cam = LumeneraCamera(1)
stage = ProScan("COM9")
CWL = CameraWithLocation(cam, stage)
CWL.show_gui(blocking=False)
CWL.af_step_size = 0.5
# Display laser shutter control
shutter = ThorLabsSC10('COM1')
shutter.show_gui(blocking=False)

# create spectrometer object
spectrometer = OceanOpticsSpectrometer(0)
spectrometer.show_gui(blocking=False)

# Display white light shutter control
whiteShutter = Uniblitz("COM7")
whiteShutter.show_gui(blocking=False)

# openshamdor
shamdor = Shamdor()
shamdor.HSSpeed = HSSpeed  # includes readout rate 2: 50kHz


shamdor.Shutter = (1, 5, 30, 30)
shamdor.SetTemperature = -90
#shamdor.CoolerON()
shamdor.ImageFlip = (1, 0)

# aligner
aligner = SpectrometerAligner(spectrometer, stage)

# pyvcam
pvc.init_pvcam()
pcam = next(Camera.detect_camera())
pcam.open()

# create equipment dict
equipment_dict = {'spectrometer': spectrometer,
                  'laser_shutter': shutter,
                  'white_shutter': whiteShutter,
                  'andor': shamdor,
                  'shamrock': shamdor.shamrock,
                  'aligner': aligner}
wizard = TrackingWizard(CWL, equipment_dict, )
shamdor.show_gui(blocking=False)
wizard.data_file.show_gui(blocking=False, )
wizard.show()
wizard.tile_edge = tile_edge_width_to_ignore


def OpenWhiteLightShutter():
    print("Matt: Opening white light shutter.")
    whiteShutter.open_shutter()


def CloseWhiteLightShutter():
    print("Matt: Closing white light shutter.")
    whiteShutter.close_shutter()


def OpenLaserShutter():
    print("Matt: Opening the laser shutter.")
    shutter.open_shutter()


def CloseLaserShutter():
    print("Matt: Closing the laser shutter.")
    shutter.close_shutter()


def ShiftFocus(z_offset=donut_z_offset):
    print("Matt: Shifting focus to the donut mode.")
    here = CWL.stage.position
    CWL.stage.move(np.array([0, 0, z_offset]) + here)


def SetShamrockSlit(slitWidth=2000):
    print("Matt: Opening spectrometer slit.")
    shamdor.shamrock.SetSlit(slitWidth)
    # Wait a bit for settings to be applied
    time.sleep(1)


def TakeInfinity3Image(imageName="Infinity3_Bias_Image", group=None):
    # Infinity3
    # We're going to take a picture - best make sure we've waited a
    # moment for the focus to return
    time.sleep(0.3)
    print("Matt: Taking Infinity3 image.")
    # take a frame and ignore (for freshness)
    CWL.camera.update_latest_frame()
    image = CWL.camera.color_image()
    if group is None:
        img = wizard.particle_group.create_dataset(imageName, data=image[
                old_div(image.shape[0],2)-50: old_div(image.shape[0],2)+50,
                old_div(image.shape[1],2)-50:old_div(image.shape[1],2)+50])
    else:
        img = group.create_dataset(imageName, data=image[
                old_div(image.shape[0],2)-50: old_div(image.shape[0],2)+50,
                old_div(image.shape[1],2)-50:old_div(image.shape[1],2)+50])
    img.attrs.create("stage_position", CWL.stage.position)
    img.attrs.create("timestamp", np.string_(datetime.datetime.now().isoformat()))


def TakeAndor0Order(imageName="Raman_Bias_0Order", group=None):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.GotoZeroOrder()
    shamdor.shamrock.SetSlit(2000)
    shamdor.set_andor_parameter('HSSpeed', HSSpeed)
    shamdor.set_andor_parameter('Exposure', andor_exposure_time_0Order)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1],
    #                    (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    # wavelengths = shamdor.metadata['x_axis']
    if group is None:
        rint = wizard.particle_group.create_dataset(imageName+"_int",
                                                    data=image)
    else:
        rint = group.create_dataset(imageName+"_int", data=image)
    # wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)
    
    
def TakePVCamImage(imageName="PVCam_Bias", group=None):
    PVCam_ExposureTime = 20
    print("Matt: Taking "+imageName)
    image = pcam.get_frame(exp_time=PVCam_ExposureTime).reshape(pcam.sensor_size[::-1])
    if group is None:
        rint = wizard.particle_group.create_dataset(imageName+"_int",
                                                    data=image)
    else:
        rint = group.create_dataset(imageName+"_int", data=image)
    rint.attrs.create("Integration time", PVCam_ExposureTime)


def TakeAndorSpec(imageName="Raman_Bias_Spectrum", group=None):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.SetWavelength(centre_wavelength)
    shamdor.shamrock.SetSlit(100)
    shamdor.set_andor_parameter('HSSpeed', HSSpeed)
    shamdor.set_andor_parameter('Exposure', andor_exposure_time_spec)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1],
    #                    (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    wavelengths = shamdor.metadata['x_axis']
    if group is None:
        rint = wizard.particle_group.create_dataset(imageName+"_int",
                                                    data=image)
        wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    else:
        rint = group.create_dataset(imageName+"_int", data=image)
        group.create_dataset(imageName+"_wl", data=wavelengths)

    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)


def TakeOOpticsSpec(imageName="OceanOptics_DF_Spectrum", group=None):
    spectrum = spectrometer.read_processed_spectrum()
    wavelengths = spectrometer.read_wavelengths()
    if group is None:
        rint = wizard.particle_group.create_dataset(imageName+"_int",
                                                    data=spectrum)
        wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    else:
        rint = group.create_dataset(imageName+"_int", data=spectrum)
        group.create_dataset(imageName+"_wl", data=wavelengths)

    rint.attrs.create("Integration time", spectrometer.integration_time)


def MoveToBackgroundLoc():
    # Move stage slightly to take background. (MAKE THIS ACTUALLY LOOK
    # FOR A SPOT WITH NO PARTICLES) TODO !!!!!!!!!!!!!!!!
    # For now, move the stage by 3 microns to a (hopefully) empty area
    CWL.stage.move_rel([3, 0, 0])


def ReturnFromBackgroundLoc():
    # Move stage slightly to return from background.
    # TODO: make this smarter
    # For now, move the stage by -3 microns to return to particle
    CWL.stage.move_rel([-3, 0, 0])


def Autofocus(useThumbnail=True):
    CWL.autofocus(use_thumbnail=useThumbnail)


def CentreOnFeature(NP_image=None, ignore_z_pos=True):
    if NP_image is not None:
        CWL.move_to_feature(NP_image, ignore_z_pos)


def ParticleEmissionProfileTimeScan(scan_number=0, particle_number=0):
    NP_image_focused = CWL.thumb_image()  # Assumes initial alignment good
    for particle_subnumber in range(0, 10):
        print('====================================')
        print('Doing time scan iteration number ' + str(particle_subnumber) +
              '.')
        print('====================================')
        wizard.subparticle_group = wizard.particle_group.create_group(
                            'SubParticle_' + str(particle_subnumber))
        StudyParticleEmissionProfile(NP_image_focused,
                                     wizard.subparticle_group)


def ManualParticleEmissionProfileTimeScan(scan_number=0):
    NP_image_focused = CWL.thumb_image()  # Assumes initial alignment good
    wizard.scan_group = wizard.data_file.create_group(
                                    'ParticleScannerScan_' + str(scan_number))
    for particle_number in range(0, 15):
        print('====================================')
        print('Doing time scan iteration number ' + str(particle_number) + '.')
        print('====================================')
        wizard.particle_group = wizard.scan_group.create_group(
                                'Particle_' + str(particle_number))
        StudyParticleEmissionProfile(NP_image_focused, wizard.particle_group)


def BleachingRecoveryTimeScan(scan_number=0, particle_number=0):
    NP_image_focused = CWL.thumb_image()  # Assumes initial alignment good
    numberOfLoops = 1
    for loopNumber in range(0, numberOfLoops):
        print('====================================')
        print('Doing loop iteration number ' + str(loopNumber) + '.')
        print('====================================')
        wizard.loop_group = wizard.particle_group.create_group(
                            'Loop_' + str(loopNumber))
        StudyParticleEmissionProfile(NP_image_focused, wizard.loop_group)
        for particle_subnumber in range(0, 20):
            print('====================================')
            print('Doing time scan iteration number ' +
                  str(particle_subnumber) + '.')
            print('====================================')
            wizard.subparticle_group = wizard.loop_group.create_group(
                                'SubParticle_' + str(particle_subnumber))
            QuickStudyParticleEmissionProfile(NP_image_focused,
                                              wizard.subparticle_group)
        numberOfMinutesToWait = 5
        if loopNumber != (numberOfLoops-1):
            for sleepLoop in range(0, numberOfMinutesToWait):
                Autofocus()
                CentreOnFeature(NP_image_focused)
                print("Waiting for bleached molecules to recover. " +
                      "Sleeping another " +
                      str(numberOfMinutesToWait-sleepLoop) + " minutes.")
                CloseWhiteLightShutter()
                time.sleep(60)  # Allow bleached molecules to recover
                OpenWhiteLightShutter()


def TestScan(NP_image_focused=None, group=None):
    OpenWhiteLightShutter()
    SetShamrockSlit(2000)
    Autofocus()
    if NP_image_focused is None:
        NP_image_focused = CWL.thumb_image()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    TakeInfinity3Image("Infinity3_FirstWhiteLight_Image", group)
    time.sleep(1)


def StudyParticleEmissionProfile(NP_image_focused=None, group=None):
    OpenWhiteLightShutter()
    SetShamrockSlit(2000)
    # NP_image_unfocused = CWL.thumb_image()
    Autofocus()
    if NP_image_focused is None:
        NP_image_focused = CWL.thumb_image()
    CentreOnFeature(NP_image_focused)
    CloseWhiteLightShutter()
    TakeInfinity3Image("Infinity3_Bias_Image", group)
    TakeAndor0Order("Raman_Bias_0Order", group)
    TakeAndorSpec("Raman_Bias_Spectrum", group)
    OpenWhiteLightShutter()
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    TakeInfinity3Image("Infinity3_FirstWhiteLight_Image", group)
    TakeAndor0Order("Raman_White_Light_0Order", group)
    TakeAndorSpec("Raman_White_Light_Spectrum", group)
    TakeOOpticsSpec(group=group)
    Autofocus()
    # CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakeAndor0Order("Raman_Laser_0Order", group)
    TakeAndorSpec("Raman_Laser_Spectrum", group)
    CloseLaserShutter()
    OpenWhiteLightShutter()
    Autofocus()
    ShiftFocus(donut_z_offset)
    TakeAndor0Order("Raman_White_Light_0Order_b", group)
    TakeOOpticsSpec(group=group)
    TakeInfinity3Image("Infinity3_SecondWhiteLight_Image", group)
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    MoveToBackgroundLoc()
    TakeInfinity3Image("Infinity3_FirstBkgndWhiteLight_Image", group)
    TakeAndor0Order("Raman_White_Light_Bkgnd_0Order", group)
    TakeAndorSpec("Raman_White_Light_Bkgnd_Spectrum", group)
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    MoveToBackgroundLoc()
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakeAndor0Order("Raman_Laser_0Order_atBkgndLoc", group)
    TakeAndorSpec("Raman_Laser_Spectrum_atBkgndLoc", group)
    CloseLaserShutter()
    OpenWhiteLightShutter()
    TakeInfinity3Image("Infinity3_SecondWhiteLight_atBkgndLoc_Image", group)
    CloseLaserShutter()
    OpenWhiteLightShutter()
    ReturnFromBackgroundLoc()
    # CentreOnFeature(NP_image_focused)
    # CentreOnFeature(NP_image_unfocused, ignore_z_pos=False)


def QuickStudyParticleEmissionProfile(NP_image_focused=None, group=None):
    OpenWhiteLightShutter()
    SetShamrockSlit(2000)
    # NP_image_unfocused = CWL.thumb_image()
    Autofocus()
    if NP_image_focused is None:
        NP_image_focused = CWL.thumb_image()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    TakeInfinity3Image("Infinity3_FirstWhiteLight_Image", group)
    TakeOOpticsSpec("OceanOptics_DF_Spectrum", group)
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakeAndor0Order("Raman_Laser_0Order", group)
    CloseLaserShutter()
    OpenWhiteLightShutter()
    TakeInfinity3Image("Infinity3_SecondWhiteLight_Image", group)
    # CentreOnFeature(NP_image_focused)
    # CentreOnFeature(NP_image_unfocused, ignore_z_pos=False)


def QuickPVCam(NP_image_focused=None, group=None):
    OpenWhiteLightShutter()
    Autofocus()
    if NP_image_focused is None:
        NP_image_focused = CWL.thumb_image()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    TakeInfinity3Image("Infinity3_FirstWhiteLight_Image", group)
    TakeOOpticsSpec("OceanOptics_DF_Spectrum", group)
    TakePVCamImage("PVCam_White", group)
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakePVCamImage("PVCam_Laser", group)
    CloseLaserShutter()
    OpenWhiteLightShutter()
    TakeInfinity3Image("Infinity3_SecondWhiteLight_Image", group)
    
    
def QuickSpatialScan(NP_image_focused=None, group=None, scan_number=0, particle_number=0):
    step_size = 0.3
    steps = 3
    initial_position = stage.get_position() # array   
    z = initial_position[2]
    xs = np.linspace(0,(steps)*step_size, num = steps, endpoint = True) - old_div((steps)*step_size,2)
    ys = np.linspace(0,(steps)*step_size, num = steps, endpoint = True) - old_div((steps)*step_size,2)
    xs += initial_position[0]
    print(xs)
    ys += initial_position[1]
    counter = -1
    places = []

    for y in ys:
        counter+=1
        if counter%2 ==0:
            for x in xs:
                places.append([x,y,z])
        else:
            for x in xs[::-1]:
                places.append([x,y,z])
    
    places = np.asarray(places) 
  
    NP_image_focused = CWL.thumb_image()  # Assumes initial alignment good
    for i, place in enumerate(places):
        print('====================================')
        print('Doing spatial scan step number ' + str(i) + ' of ' + str(len(places)) +
              '.')
        print('====================================')
        wizard.subparticle_group = wizard.particle_group.create_group(
                            'SubParticle_' + str(i))
        stage.move(place)
        print(place-initial_position)       
        time.sleep(0.5)
        QuickPVCam(NP_image_focused, wizard.subparticle_group)

    stage.move(initial_position)
