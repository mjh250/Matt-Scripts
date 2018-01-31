# -*- coding: utf-8 -*-
"""
Created on Tue Jan 30 13:29:18 2018

@author: np-albali
"""
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.camera.camera_with_location import CameraWithLocation
from nplab.instrument.spectrometer.spectrometer_aligner import SpectrometerAligner
from nplab.instrument.stage.prior import ProScan
from nplab.instrument.shutter.thorlabs_sc10 import ThorLabsSC10
from nplab.instrument.shutter.BX51_uniblitz import Uniblitz
from nplab.instrument.spectrometer.shamdor import Shamdor
from particle_tracking_app.particle_tracking_wizard import TrackingWizard
import numpy as np
import datetime
import time


HSSpeed = 2
andor_exposure_time_0Order = 5
andor_exposure_time_spec = andor_exposure_time_0Order*3
donut_z_offset = -0.0  # z offset to move focus to the donut mode

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
shamdor.CoolerON()
shamdor.ImageFlip = (1, 0)

# aligner
aligner = SpectrometerAligner(spectrometer, stage)

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


def OpenWhiteLightShutter():
    print("Matt: Opening white light shutter.")
    whiteShutter.open_shutter()
    return


def CloseWhiteLightShutter():
    print("Matt: Closing white light shutter.")
    whiteShutter.close_shutter()
    return


def OpenLaserShutter():
    print("Matt: Opening the laser shutter.")
    shutter.open_shutter()
    return


def CloseLaserShutter():
    print("Matt: Closing the laser shutter.")
    shutter.close_shutter()
    return


def ShiftFocus(z_offset=donut_z_offset):
    print("Matt: Shifting focus to the donut mode.")
    here = CWL.stage.position
    CWL.stage.move(np.array([0, 0, z_offset]) + here)
    return


def SetShamrockSlit(slitWidth=2000):
    print("Matt: Opening spectrometer slit.")
    shamdor.shamrock.SetSlit(slitWidth)
    # Wait a bit for settings to be applied
    time.sleep(1)
    return


def TakeInfinity3Image(imageName="Infinity3_Bias_Image"):
    # Infinity3
    # We're going to take a picture - best make sure we've waited a
    # moment for the focus to return
    time.sleep(0.3)
    print("Matt: Taking Infinity3 bias image.")
    # take a frame and ignore (for freshness)
    CWL.camera.update_latest_frame()
    image = CWL.camera.color_image()
    img = wizard.particle_group.create_dataset(imageName, data=image[
            image.shape[0]/2-50: image.shape[0]/2+50,
            image.shape[1]/2-50:image.shape[1]/2+50])
    img.attrs.create("stage_position", CWL.stage.position)
    img.attrs.create("timestamp", datetime.datetime.now().isoformat())
    return


def TakeAndor0Order(imageName="Raman_Bias_0Order"):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.GotoZeroOrder()
    shamdor.shamrock.SetSlit(2000)
    shamdor.SetParameter('Exposure', andor_exposure_time_0Order)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1], (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    wavelengths = shamdor.shamrock.GetWavelength()
    rint = wizard.particle_group.create_dataset(imageName+"_int", data=image)
    wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)
    return


def TakeAndorSpec(imageName="Raman_Bias_Spectrum"):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.SetWavelength(shamdor.shamrock.center_wavelength)
    shamdor.shamrock.SetSlit(100)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1], (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    wavelengths = shamdor.shamrock.GetWavelength()
    rint = wizard.particle_group.create_dataset(imageName+"_int", data=image)
    wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)
    return


def MoveToBackgroundLoc():
    # Move stage slightly to take background. (MAKE THIS ACTUALLY LOOK
    # FOR A SPOT WITH NO PARTICLES) TODO !!!!!!!!!!!!!!!!
    # For now, move the stage by 3 microns to a (hopefully) empty area
    CWL.stage.move_rel([3, 0, 0])
    return


# def CentreOnFeature():
#     CWL.move_to_feature