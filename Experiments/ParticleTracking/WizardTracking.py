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
import nplab.datafile as df
from particle_tracking_app.particle_tracking_wizard import TrackingWizard
import numpy as np
import datetime
import time


HSSpeed = 2
andor_exposure_time_0Order =5
andor_exposure_time_spec = andor_exposure_time_0Order*3
donut_z_offset = 0  # z offset to move focus to the donut mode
centre_wavelength = 675

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


def TakeAndor0Order(imageName="Raman_Bias_0Order"):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.GotoZeroOrder()
    shamdor.shamrock.SetSlit(2000)
    shamdor.SetParameter('HSSpeed', HSSpeed)
    shamdor.SetParameter('Exposure', andor_exposure_time_0Order)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1],
    #                    (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    # wavelengths = shamdor.metadata['x_axis']
    rint = wizard.particle_group.create_dataset(imageName+"_int", data=image)
    # wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)


def TakeAndorSpec(imageName="Raman_Bias_Spectrum"):
    print("Matt: Taking "+imageName)
    shamdor.shamrock.SetWavelength(centre_wavelength)
    shamdor.shamrock.SetSlit(100)
    shamdor.SetParameter('HSSpeed', HSSpeed)
    shamdor.SetParameter('Exposure', andor_exposure_time_spec)
    time.sleep(5)
    # image = np.reshape(shamdor.raw_snapshot()[1],
    #                    (-1, shamdor.shamrock.pixel_number))
    image = shamdor.raw_snapshot()[1]
    wavelengths = shamdor.metadata['x_axis']
    rint = wizard.particle_group.create_dataset(imageName+"_int", data=image)
    wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)
    # rint.attrs.create("Laser power", raman.laser_power)
    rint.attrs.create("Slit size", shamdor.shamrock.slit_width)
    rint.attrs.create("Integration time", shamdor.Exposure)
    # rint.attrs.create("description", raman.scan_desc)


def TakeOOpticsSpec(imageName="OceanOptics_DF_Spectrum"):
    spectrum = spectrometer.read_processed_spectrum()
    wavelengths = spectrometer.read_wavelengths()
    rint = wizard.particle_group.create_dataset(imageName+"_int",
                                                data=spectrum)
    wizard.particle_group.create_dataset(imageName+"_wl", data=wavelengths)

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
    for particle_subnumber in range(0, 15):
        print('====================================')
        print('Doing time scan iteration number ' + str(particle_subnumber) + '.')
        print('====================================')
        wizard.subparticle_group = wizard.particle_group.create_group(
                            'SubParticle_' + str(particle_subnumber))  # THIS LINE CREATES GROUP, BUT DATA NOT GETTING STORED THERE NEED TO CHANGE EACH FUNTION THAT STORES DATA?
        StudyParticleEmissionProfile(NP_image_focused)


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
        StudyParticleEmissionProfile(NP_image_focused)


def StudyParticleEmissionProfile(NP_image_focused=None):
    OpenWhiteLightShutter()
    SetShamrockSlit(2000)
    #NP_image_unfocused = CWL.thumb_image()
    Autofocus()
    if NP_image_focused is None:
        NP_image_focused = CWL.thumb_image()
    CentreOnFeature(NP_image_focused)
    CloseWhiteLightShutter()
    TakeInfinity3Image("Infinity3_Bias_Image")
    TakeAndor0Order("Raman_Bias_0Order")
    TakeAndorSpec("Raman_Bias_Spectrum")
    OpenWhiteLightShutter()
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    TakeInfinity3Image("Infinity3_FirstWhiteLight_Image")
    TakeAndor0Order("Raman_White_Light_0Order")
    TakeAndorSpec("Raman_White_Light_Spectrum")
    TakeOOpticsSpec()
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakeAndor0Order("Raman_Laser_0Order")
    TakeAndorSpec("Raman_Laser_Spectrum")
    CloseLaserShutter()
    OpenWhiteLightShutter()
    TakeInfinity3Image("Infinity3_SecondWhiteLight_Image")
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    MoveToBackgroundLoc()
    TakeInfinity3Image("Infinity3_FirstBkgndWhiteLight_Image")
    TakeAndor0Order("Raman_White_Light_Bkgnd_0Order")
    TakeAndorSpec("Raman_White_Light_Bkgnd_Spectrum")
    Autofocus()
    CentreOnFeature(NP_image_focused)
    ShiftFocus(donut_z_offset)
    MoveToBackgroundLoc()
    CloseWhiteLightShutter()
    OpenLaserShutter()
    TakeAndor0Order("Raman_Laser_0Order_atBkgndLoc")
    TakeAndorSpec("Raman_Laser_Spectrum_atBkgndLoc")
    CloseLaserShutter()
    OpenWhiteLightShutter()
    TakeInfinity3Image("Infinity3_SecondWhiteLight_atBkgndLoc_Image")
    CloseLaserShutter()
    OpenWhiteLightShutter()
    ReturnFromBackgroundLoc()
    #CentreOnFeature(NP_image_focused)
    #CentreOnFeature(NP_image_unfocused, ignore_z_pos=False)
