# -*- coding: utf-8 -*-
"""
Created on Fri Oct 13 13:44:27 2017

@author: np-albali
"""

from F_SERS_rig.ExperimentalScripts import view_spectra, Raman_class
from F_SERS_rig.Equipment.FW212C import FW212C
from F_SERS_rig.ExperimentalScripts.ParticleScanner import ParticleScanner
from nplab.instrument.camera.lumenera import LumeneraCamera
from nplab.instrument.stage.prior import ProScan
from nplab.instrument.stage.camera_stage_mapper import CameraStageMapper
from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer
from nplab.instrument.shutter.BX51_uniblitz import Uniblitz
from nplab.instrument.spectrometer.spectrometer_aligner\
    import SpectrometerAligner
import nplab.datafile as df
import threading
import numpy as np
import time
import datetime


class FarField_ParticleScanner(ParticleScanner):

    def go_to_particles(self, payload_function=lambda: time.sleep(2),
                        background=True, max_n_particles=None):
            """Find particles, then visit each one and execute a payload.

            This function returns immediately as it spawns a background thread.
            The scan can be monitored through traits scan_status,
            scan_progress, scanning. It can be aborted with the abort_scan()
            method.

            By default it simply waits for 2 seconds at each position.
            """
            if self.scanning:
                return

            def worker_function():
                if not self._scan_lock.acquire(False):
                    raise Exception("Tried to start a scan, \
                    but one was in progress!")
                aborted = False
                self.scanning = True
                self.scan_progress = 0
                self.scan_status = "Setting up scan..."
                here = self.csm.camera_centre_position()
                pixel_positions = self.find_particles()
                positions = [self.csm.camera_pixel_to_sample(p) for p in
                             pixel_positions]
                image = self.csm.camera.color_image()
                feature_images = [image
                                  [p[0]-self.border_pixels:p[0]+self.border_pixels,
                                   p[1]-self.border_pixels:p[1]+self.border_pixels]
                                  for p in pixel_positions]  # get feature images
                for index, p in enumerate(positions):
                    if max_n_particles is not None and index >= max_n_particles:
                        print "Terminating scan as we've now scanned enough \
                        particles"
                        break
                    self.scan_status = "Scanning particle \
                                        %d of %d" % (index, len(positions))
                    self.csm.move_to_sample_position(p)
                    time.sleep(0.3)
                    self.current_feature_image = feature_images[index]
                    self.csm.centre_on_feature(feature_images[index])
                    payload_function()
                    self.scan_progress = float(index)/float(len(positions))*100
                    if self._abort_scan_event.is_set():  # for aborting a scan
                        self.scan_status = "Scan Aborted."
                        self._abort_scan_event.clear()
                        aborted = True
                        break
                self.csm.move_to_sample_position(here)
                self.scan_status = "Scan Finished"
                self.scan_progress = 100.0
                print "Scan Finished :)"
                self.scanning = False
                self._scan_lock.release()
                return not aborted
            # execute the above function in the background
            if background:
                self._scan_thread = threading.Thread(target=worker_function)
                self._scan_thread.start()
            else:  # if we elected not to use a thread, just do it!
                return worker_function()

    def pf_align_and_take_z_scan(self, dz=np.arange(-4, 4, 0.4),
                                 datafile_group=None):
        """Set up for a scan of all particles, then return a payload function.

        The "payload function" is suitable for the eponymous argument of
        go_to_particles, and will autofocus, align particle to fibre, and take
        a Z stack.  NB the payload function "wraps up" the arguments neatly so
        we don't need to store things like the depth of the Z stack.
        """
        if datafile_group is None:
            datafile_group = self.new_data_group("particleScans/scan%d",
                                                 self.datafile)

        def align_and_take_z_scan():
            # Initialize shutter positions by making sure they are all closed
            # print("Matt: Initializing shutter positions.")
            # light_shutter.close_shutter()
            # raman.shutter.close_shutter()
            # Define scan parameters:
            HSSpeed = 2
            andor_exposure_time_0Order = 5
            andor_exposure_time_spec = andor_exposure_time_0Order*3
            donut_z_offset = -0.0  # z offset to move focus to the donut mode
            # Open white light shutter and close laser shutter
            print("Matt: Opening white light shutter.")
            light_shutter.open_shutter()
            # Fully open Shamrock slit
            print("Matt: Opening spectrometer slit.")
            raman.sham.SetSlit(2000)
            # Wait a bit for settings to be applied
            time.sleep(1)

            # Do Autofocus
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5, 2.5, 0.5))
            # short integration time for alignment
            # self.aligner.spectrometer.integration_time = 1000. # aligns on OO
            # self.aligner.optimise_2D(tolerance=0.03, stepsize=0.2)
            self.csm.centre_on_feature(self.current_feature_image)
            # Shift focus to the donut mode
            print("Matt: Shifting focus to the donut mode.")
            here = self.csm.stage.position
            self.csm.stage.move(np.array([0, 0, donut_z_offset]) + here)
            # long integration time for measurement
            # self.aligner.spectrometer.integration_time = 3000.

            # Initialize datafile
            print("Matt: Initializing datafile.")  # TODO: Check that this is needed for autofocus
            g = self.new_data_group("scan_%d", datafile_group)
            dset = g.create_dataset("scan",
                                    data=self.aligner.z_scan(dz))
            for key, val in self.aligner.spectrometer.\
                    get_metadata().iteritems():
                dset.attrs.create(key, val)
            dset.attrs.create("stage_position", self.csm.stage.position)
            dset.attrs.create("camera_centre_position",
                              self.csm.camera_centre_position())
            dset.attrs.create("timestamp", datetime.datetime.now().isoformat())
            dset.attrs.create("dz", dz)
            # Close all shutters
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # Set Infinity3 camera settings --- TODO !!!!!!!!!!!!!!!!!!!!!!!!
            # Set Andor camera settings --- TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!
            raman.cam.SetParameter('HSSpeed', HSSpeed)
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take bias images on Infinity3, Andor, and OceanOptics
            # Infinity3
            print("Matt: Taking Infinity3 bias image.")
            # take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_Bias_Image", data=image[
                    image.shape[0]/2-50: image.shape[0]/2+50,
                    image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            # Andor 0 order
            print("Matt: Taking Andor 0 order bias image.")
            raman.sham.GotoZeroOrder()
            raman.sham.SetSlit(2000)
            raman.cam.SetParameter('Exposure', andor_exposure_time_0Order)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Bias_0Order_int", data=image)
            g.create_dataset("Raman_Bias_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Andor spectrum
            print("Matt: Taking Andor spectrum bias image.")
            raman.sham.SetWavelength(raman.centre_Wavelength)
            raman.sham.SetSlit(100)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Bias_Spectrum_int", data=image)
            g.create_dataset("Raman_Bias_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # OceanOptics
            print("Matt: Taking OceanOptics bias image.")
            (oowl, oospec) = spectrometer.read()
            g.create_dataset("OOptics_Bias_Spectrum_int", data=oospec)
            g.create_dataset("OOptics_Bias_Spectrum_wl", data=oowl)
            # Turn on white light
            print("Matt: Turning white light back on.")
            light_shutter.open_shutter()

            # Do Autofocus
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5, 2.5, 0.5))
            # short integration time for alignment
#            self.aligner.spectrometer.integration_time = 1000.
#            self.aligner.optimise_2D(tolerance=0.03, stepsize=0.2)
            self.csm.centre_on_feature(self.current_feature_image)
            # Shift focus to the donut mode
            print("Matt: Shifting focus to the donut mode.")
            here = self.csm.stage.position
            self.csm.stage.move(np.array([0, 0, donut_z_offset]) + here)
            # long integration time for measurement
            self.aligner.spectrometer.integration_time = 3000.

            # Take Infinity3 image
            print("Matt: Taking first Infinity3 white light image.")
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset(
                    "Infinity3_FirstWhiteLight_Image",
                    data=image[image.shape[0]/2-50: image.shape[0]/2+50,
                               image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            # Take white light spectrum (Ocean Optics) TODO !!!!!!!!!!!!!!
            # Take white light spectrum (Andor)
            print("Matt: Taking white light spectrum on Andor")
            raman.sham.SetWavelength(raman.centre_Wavelength)
            raman.sham.SetSlit(100)
            raman.cam.SetParameter('Exposure', andor_exposure_time_spec)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_Spectrum_int",
                                    data=image)
            g.create_dataset("Raman_White_Light_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Take white light image (Andor)
            print("Matt: Taking Andor 0 order white light image.")
            raman.sham.GotoZeroOrder()
            raman.sham.SetSlit(2000)
            raman.cam.SetParameter('Exposure', andor_exposure_time_0Order)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_0Order_int",
                                    data=image)
            g.create_dataset("Raman_White_Light_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)

            # Do Autofocus
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5, 2.5, 0.5))
            # short integration time for alignment
#            self.aligner.spectrometer.integration_time = 1000.
#            self.aligner.optimise_2D(tolerance=0.03, stepsize=0.2)
            self.csm.centre_on_feature(self.current_feature_image)
            # Shift focus to the donut mode
            print("Matt: Shifting focus to the donut mode.")
            here = self.csm.stage.position
            self.csm.stage.move(np.array([0, 0, donut_z_offset]) + here)
            # long integration time for measurement
            self.aligner.spectrometer.integration_time = 3000.

            # Turn off white light
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # Turn on laser
            print("Matt: Opening the laser shutter.")
            raman.shutter.open_shutter()
            # Set Infinity3 exposure/gain very low and image beam profile,
            # then restore old values
            oldExposure = cam.exposure
            oldGain = cam.gain
            cam.exposure = 0
            cam.gain = 0
            print("Matt: Taking Infinity3 image of laser beam profile.")
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_Laser_Beam_Image", data=image[
                    image.shape[0]/2-50: image.shape[0]/2+50,
                    image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            cam.exposure = oldExposure
            cam.gain = oldGain
            # Take image of laser zero-order (Andor)
            print("Matt: Taking Andor 0 order laser image.")
            raman.sham.GotoZeroOrder()
            raman.sham.SetSlit(2000)
            raman.cam.SetParameter('Exposure', andor_exposure_time_0Order)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_0Order_int", data=image)
            g.create_dataset("Raman_Laser_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Take image of laser spectrum (Andor)
            print("Matt: Taking Andor spectrum laser image.")
            raman.sham.SetWavelength(raman.centre_Wavelength)
            raman.sham.SetSlit(100)
            raman.cam.SetParameter('Exposure', andor_exposure_time_spec)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_Spectrum_int", data=image)
            g.create_dataset("Raman_Laser_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Turn off laser
            print("Matt: Closing the laser shutter.")
            raman.shutter.close_shutter()
            # Open the white light shutter
            print("Matt: Turning white light back on.")
            light_shutter.open_shutter()
            # Take second Infinity3 white light image (to track drift)
            print("Matt: Taking second Infinity3 white light image.")
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset(
                    "Infinity3_SecondWhiteLight_Image",
                    data=image[image.shape[0]/2-50:image.shape[0]/2+50,
                               image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())

            # Do Autofocus
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5, 2.5, 0.5))
            # short integration time for alignment
#            self.aligner.spectrometer.integration_time = 1000.
#            self.aligner.optimise_2D(tolerance=0.03, stepsize=0.2)
            self.csm.centre_on_feature(self.current_feature_image)
            # Shift focus to the donut mode
            print("Matt: Shifting focus to the donut mode.")
            here = self.csm.stage.position
            self.csm.stage.move(np.array([0, 0, donut_z_offset]) + here)
            # long integration time for measurement
            self.aligner.spectrometer.integration_time = 3000.

            # Move stage slightly to take background. (MAKE THIS ACTUALLY LOOK
            # FOR A SPOT WITH NO PARTICLES) TODO !!!!!!!!!!!!!!!!
            # For now, move the stage by 3 microns to a (hopefully) empty area
            self.csm.stage.move_rel([3, 0, 0])
            # Take Infinity3 white light image (as evidence stage moved to a
            # good location for background)
            print(("Matt: Taking first Infinity3 white light image for "
                  "background location."))
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset(
                    "Infinity3_FirstBkgndWhiteLight_Image",
                    data=image[image.shape[0]/2-50:image.shape[0]/2+50,
                               image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            # Take white light spectrum at background location (Ocean Optics)
            # - TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # Take white light spectrum at background location (Andor)
            print(("Matt: Taking white light spectrum at background location "
                  "on Andor"))
            raman.sham.SetWavelength(raman.centre_Wavelength)
            raman.sham.SetSlit(100)
            raman.cam.SetParameter('Exposure', andor_exposure_time_spec)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset(
                    "Raman_White_Light_Bkgnd_Spectrum_int", data=image)
            g.create_dataset(
                    "Raman_White_Light_Bkgnd_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Take white light image at background location (Andor)
            print(("Matt: Taking Andor 0 order white light image at "
                  "background location."))
            raman.sham.GotoZeroOrder()
            raman.sham.SetSlit(2000)
            raman.cam.SetParameter('Exposure', andor_exposure_time_0Order)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset(
                    "Raman_White_Light_Bkgnd_0Order_int", data=image)
            g.create_dataset(
                    "Raman_White_Light_0Order_Bkgnd_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)

            # Do Autofocus
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5, 2.5, 0.5))
            # short integration time for alignment
#            self.aligner.spectrometer.integration_time = 1000.
#            self.aligner.optimise_2D(tolerance=0.03, stepsize=0.2)
            self.csm.centre_on_feature(self.current_feature_image)
            # Shift focus to the donut mode
            print("Matt: Shifting focus to the donut mode.")
            here = self.csm.stage.position
            self.csm.stage.move(np.array([0, 0, donut_z_offset]) + here)
            # long integration time for measurement
            self.aligner.spectrometer.integration_time = 3000.

            # Turn off white light
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # Turn on laser
            print("Matt: Opening the laser shutter.")
            raman.shutter.open_shutter()
            # Set Infinity3 exposure/gain very low and image beam profile,
            # then restore old values
            oldExposure = cam.exposure
            oldGain = cam.gain
            cam.exposure = 0
            cam.gain = 0
            print(("Matt: Taking Infinity3 image of laser beam profile at "
                  "background location."))
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset(
                    "Infinity3_Laser_Beam_Image_atBkgndLoc", data=image[
                            image.shape[0]/2-50: image.shape[0]/2+50,
                            image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            cam.exposure = oldExposure
            cam.gain = oldGain
            # Take image of laser zero-order at background location (Andor)
            print(("Matt: Taking Andor 0 order laser image at background "
                  "location."))
            raman.sham.GotoZeroOrder()
            raman.sham.SetSlit(2000)
            raman.cam.SetParameter('Exposure', andor_exposure_time_0Order)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset(
                    "Raman_Laser_0Order_atBkgndLoc_int", data=image)
            g.create_dataset(
                    "Raman_Laser_0Order_atBkgndLoc_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Take image of laser spectrum at background location (Andor)
            print(("Matt: Taking Andor spectrum laser image at background "
                  "location."))
            raman.sham.SetWavelength(raman.centre_Wavelength)
            raman.sham.SetSlit(100)
            raman.cam.SetParameter('Exposure', andor_exposure_time_spec)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(), (-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset(
                    "Raman_Laser_Spectrum_atBkgndLoc_int", data=image)
            g.create_dataset(
                    "Raman_Laser_Spectrum_atBkgndLoc_wl", data=wavelengths)
            rint.attrs.create("Laser power", raman.laser_power)
            rint.attrs.create("Slit size", raman.slit_size)
            rint.attrs.create("Integration time", raman.Integration_time)
            rint.attrs.create("description", raman.scan_desc)
            # Turn off laser
            print("Matt: Closing the laser shutter.")
            raman.shutter.close_shutter()
            # Open the white light shutter
            print("Matt: Turning white light back on.")
            light_shutter.open_shutter()
            # Take second Infinity3 white light image at background location
            # (to track drift)
            print(("Matt: Taking second Infinity3 white light image at "
                  "background location."))
            # We're going to take a picture - best make sure we've waited a
            # moment for the focus to return
            time.sleep(0.3)
            # Take a frame and ignore (for freshness)
            self.csm.camera.update_latest_frame()
            image = self.csm.camera.color_image()
            img = g.create_dataset(
                    "Infinity3_SecondWhiteLight_atBkgndLoc_Image",
                    data=image[image.shape[0]/2-50:image.shape[0]/2+50,
                               image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position", self.csm.stage.position)
            img.attrs.create("timestamp", datetime.datetime.now().isoformat())
            # Turn off all light sources
            print(("Matt: Opening the white light shutter and closing the "
                  "laser shutter."))
            light_shutter.open_shutter()
            raman.shutter.close_shutter()

            # POTENTIALLY USEFUL: LASER FOCUS
            # "first bring the laser into focus"
            # here = stage.position
            # laser_focus = raman.AlignHeNe(dset)  # remove, keep, or offset
            # print "Moving to HeNe Focus (%g)" % (laser_focus)
            # stage.move_rel([0,0,laser_focus])
            # time.sleep(1)

            ##################################################################
            print("Reached the end of Matt's code.")
            datafile_group.file.flush()
            ##################################################################

        return align_and_take_z_scan


if __name__ == '__main__':
    cam = LumeneraCamera(1)
    stage = ProScan("COM9")
    light_shutter = Uniblitz("COM7")
    stage.query("SERVO 1")  # set up stage to use servocontrol
    stage.query("UPR Z 500")  # 500 um per revolution = lab 6
    mapper = CameraStageMapper(cam, stage)
    mapper.frames_to_discard = 2  # Infinity 3 is faster: must discard frames
    spectrometer = OceanOpticsSpectrometer(0)
    spectrometer.set_tec_temperature(-21)
    aligner = SpectrometerAligner(spectrometer, stage)
    raman = Raman_class.Raman_spec()
    filter_wheel = FW212C()
    scanner = FarField_ParticleScanner(mapper, spectrometer, aligner,
                                       df.current(), raman)
    viewer = view_spectra.ScanViewer(scanner.datafile)
    spectrometer.show_gui(blocking=False)
    cam.show_gui(False)
    aligner.edit_traits()
    viewer.edit_traits()

    raman.cam.show_gui(blocking=False)
    raman.CCD_cooling(-90)

    raman.edit_traits()
    raman.cam.verbosity = False

    stop_tile = threading.Event()
    scanner.edit_traits()
    mapper.edit_traits()
