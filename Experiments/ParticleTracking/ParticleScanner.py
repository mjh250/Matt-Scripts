
# -*- coding: utf-8 -*-
"""
"""
from CloudyScripts import camera, prior_stage, camera_stage_mapper, seabreeze, spectrometer_aligner, Uniblitz_Shutter, view_spectra, Raman_class
import traits
from traits.api import HasTraits, Float, Enum, Int, Range, Bool, Button, String, on_trait_change, Array
import traitsui
from traitsui.api import View, Item, HGroup, VGroup, Tabbed, RangeEditor
import cv2
from scipy import ndimage
import time
import h5py
import numpy as np
from numpy import inf
import datetime, re
import matplotlib.pyplot as plt
import threading
import smtplib
import sys

class ParticleScanner(HasTraits):
    median_filter_width = Range(1, 31)
    threshold_block_size = Int(21)
    threshold_level = Range(-255,255,50,mode="slider")
    live_filter = Enum(["None","Denoised","Thresholded"])
    scan_current_view = Button()
    abort_scan_button = Button(label="abort_scan")
    scan_status = String("Not Scanning")
    scan_progress = Range(0., 100., 0.)
    scanning = Bool(False)
    tiled_scan_size=Array(shape=(2,),dtype=np.int)
    start_tiled_scan=Button()
    border_pixels = 15
    
    traits_view = View(
                    Tabbed(
                        VGroup(
                            Item(name="median_filter_width"),
                            Item(name="threshold_block_size"),
                            Item(name="threshold_level"),
                            Item(name="live_filter"),
                            label="Image Processing",
                        ),
                        VGroup(
                            Item(name="scan_status",style="readonly"),
                            Item(name="scan_progress",style="readonly"),
                            Item(name="scanning",style="readonly"),
                            Item(name="scan_current_view"),
                            Item(name="abort_scan_button"),
                            Item(name="tiled_scan_size"),
                            Item(name="start_tiled_scan"),
                            label="Scan Control",
                        ),
                    ),
                    title="Particle Scanner"
                )
    
    """
    Find particles in an image and move to them
    """
    def __init__(self, camera_stage_mapper, spectrometer, spectrometer_aligner, datafile):
        super(ParticleScanner,self).__init__()
        self.csm=camera_stage_mapper
        self.spectrometer = spectrometer
        self.datafile = datafile
        self.aligner = spectrometer_aligner
        self._live_filter_changed() #enable video filter if required
        self._scan_lock = threading.Lock()
        self._abort_scan_event = threading.Event()
        
    def SendCompleteMessage(self,number):    
        gmail_user = "lab6.Raman.fb400@gmail.com"
        gmail_pwd = "NQ3dPv6SXZUEdfTE"
        FROM = 'lab6.Raman.fb400@gmail.com'
        TO = ['cc831@cam.ac.uk'] #must be a list
        SUBJECT = "Scan finished"
        TEXT = "%d particles scanned" %number
        
        # Prepare actual message
        message = """\From: %s\nTo: %s\nSubject: %s\n\n%s
        """ % (FROM, ", ".join(TO), SUBJECT, TEXT)
        try:
            server = smtplib.SMTP("smtp.gmail.com", 587) #or port 465 doesn't seem to work!
            server.ehlo()
            server.starttls()
            server.login(gmail_user, gmail_pwd)
            server.sendmail(FROM, TO, message)
            server.close()
            print 'successfully sent the mail'
        except:
            print "failed to send mail"
        
    def denoise_image(self,img):
        """apply the current denoising filter to the image"""
        if(self.median_filter_width>0):
            if self.median_filter_width % 2 == 0: self.median_filter_width += 1 #guard agains even integers!
            #return cv2.blur(img,self.median_filter_width)
            return cv2.blur(img,(self.median_filter_width,self.median_filter_width))
        else:
            return img
            
    def threshold_image(self,img):
        """apply threshold with the current settings to an image"""
        #return cv2.threshold(self.denoise_image(img),self.threshold_level,255,cv2.THRESH_BINARY)[1]
        img = cv2.adaptiveThreshold(img,255,cv2.ADAPTIVE_THRESH_MEAN_C,cv2.THRESH_BINARY,(int(self.threshold_block_size)/2)*2+1,self.threshold_level)
        kernel = np.ones((self.median_filter_width, self.median_filter_width),np.uint8)
        return cv2.morphologyEx(img, cv2.MORPH_OPEN, kernel, iterations=1) #after thresholding, erode then dilate to kill small blobs/noise

    def camera_filter_function(self,frame):
        img = cv2.cvtColor(frame,cv2.COLOR_RGB2GRAY)
        if self.live_filter == "Denoised":
            img=self.denoise_image(img)
        elif self.live_filter == "Thresholded":
            img = self.threshold_image(self.denoise_image(img))
        return cv2.cvtColor(img,cv2.COLOR_GRAY2RGB)

    def _live_filter_changed(self):
        if self.live_filter == "None":            
            self.csm.camera.filter_function = None
        else:
            self.csm.camera.filter_function = self.camera_filter_function
    @on_trait_change("find_particles")

    def find_particles_in_new_image(self): #necessary to stop extra arguments from Traits messing things up
        self.find_particles()

    def find_particles(self,img=None):
        """find particles in the supplied image, or in the camera image"""
        if img is None:
            ret, frame = self.csm.camera.raw_snapshot()
            img = self.threshold_image(self.denoise_image(
                cv2.cvtColor(frame,cv2.COLOR_RGB2GRAY)))[self.border_pixels:-self.border_pixels,self.border_pixels:-self.border_pixels] #ignore the edges
        labels, nlabels = ndimage.measurements.label(img)
        return [np.array(p)+15 for p in ndimage.measurements.center_of_mass(img, labels, range(1,nlabels+1))] #add 15 onto all the positions

    def go_to_particles(self,payload_function=lambda: time.sleep(2),background=True,max_n_particles=None):
        """Find particles, then visit each one in turn and execute a payload.
        
        This function returns immediately as it spawns a background thread. The
        scan can be monitored through traits scan_status, scan_progress, 
        scanning.  It can be aborted with the abort_scan() method.
        
        By default it simply waits for 2 seconds at each position.
        """
        if self.scanning:
            return

        def worker_function():
            if not self._scan_lock.acquire(False):
                raise Exception("Tried to start a scan, but one was in progress!")
            aborted=False
            self.scanning = True
            self.scan_progress = 0
            self.scan_status = "Setting up scan..."
            here = self.csm.camera_centre_position()
            pixel_positions = self.find_particles()
            positions = [self.csm.camera_pixel_to_sample(p) for p in pixel_positions]
            image = self.csm.camera.color_image()
            feature_images = [image[p[0]-self.border_pixels:p[0]+self.border_pixels,p[1]-self.border_pixels:p[1]+self.border_pixels] for p in pixel_positions] #extract feature images
            for index, p in enumerate(positions):
                if max_n_particles is not None and index >= max_n_particles:
                    print "Terminating scan as we've now scanned enough particles"
                    break
                self.scan_status = "Scanning particle %d of %d" % (index, len(positions))
                self.csm.move_to_sample_position(p)
                time.sleep(0.3)
                self.csm.centre_on_feature(feature_images[index])
                payload_function()
                self.scan_progress = float(index)/float(len(positions))*100
                if self._abort_scan_event.is_set(): #this event lets us abort a scan
                    self.scan_status="Scan Aborted."
                    self._abort_scan_event.clear()
                    aborted=True
                    break
            self.csm.move_to_sample_position(here)
            self.scan_status = "Scan Finished"
            self.scan_progress=100.0
            print "Scan Finished :)"
            self.scanning=False
            self._scan_lock.release()
            return not aborted
        #execute the above function in the background
        if background:
            self._scan_thread = threading.Thread(target=worker_function)
            self._scan_thread.start()
        else: #if we elected not to use a thread, just do it!
            return worker_function()
    @on_trait_change("abort_scan_button")

    def abort_scan(self):
        """Abort a currently-running scan in a background thread."""
        if self._scan_thread is not None and self._scan_thread.is_alive():
            self._abort_scan_event.set()
#        if self._scan_thread is not None:
#            self._scan_thread.join()
#        self._abort_scan_event.clear()

    def tile_scans(self, size, background = True, tile_start_function=None, ts_args=[], ts_kwargs={}, *args, **kwargs):
        def worker_function():
            grid_size = np.array(size)
            here = self.csm.camera_centre_position()
            scan_centres = [self.csm.camera_point_to_sample(np.array([i,j])-grid_size/2) 
			    	for i in range(grid_size[0]) 
			    	for j in ( range(grid_size[1]) if i % 2 == 0 else reversed(range(grid_size[1]))) #snake-style raster scanning
			    ]
            for centre in scan_centres:
                print "Taking a scan with centre %.1f, %.1f um" % tuple(centre)
                self.csm.move_to_sample_position(centre)
                if tile_start_function is not None:
                    tile_start_function(*ts_args, **ts_kwargs)
                ret = self.go_to_particles(background=False *args, **kwargs)
                if not ret:
                    print "Scan aborted!"
                    break
            self.csm.move_to_sample_position(here)
            print "Scan Finished!"
            latest_group = sorted([v for k, v in self.datafile['particleScans'].iteritems() if 'scan' in k], key=lambda g: int(re.search(r"(\d+)$",g.name).groups()[0]))[-1]
            number_of_particles = len([k for k in latest_group.keys() if 'scan_' in k])            
            self.SendCompleteMessage(number_of_particles)
        #execute the above function in the background
        if background:
            self._scan_thread = threading.Thread(target=worker_function)
            self._scan_thread.start()
        else: #if we elected not to use a thread, just do it!
            worker_function()
        
    def _scan_current_view_fired(self):
        self.take_zstacks_of_particles()

    def take_zstacks_of_particles(self, dz=np.arange(-2.5,2.5,0.2), datafile_group = None, *args, **kwargs):
        """visit each particle and scan it spectrally"""
        self.spectrometer.live_view=False
        g = self.new_data_group("particleScans/scan%d",self.datafile) if datafile_group is None else datafile_group
        g.create_dataset("Raman_wavelengths",data=raman.GetWavelength())
        self.save_overview_images(g)
        self.go_to_particles(self.pf_align_and_do_scan(dz,g), *args, **kwargs)

    def _start_tiled_scan_fired(self):
        self.take_zstacks_of_particles_tiled(self.tiled_scan_size)

    def take_zstacks_of_particles_tiled(self, shape, **kwargs):
        """Take z-stacked spectra of all the particles in several fields-of-view.
        
        We essentially run take_zstacks_of_particles for several fields of view,
        tiling them together in to the "shape" specified (2-element tuple). The
        centre of the tiled image is the current position.
        """
        self.spectrometer.live_view=False
        g = self.new_data_group("particleScans/scan%d",self.datafile)
        g.create_dataset("Raman_wavelengths",data=raman.GetWavelength())
        self.tile_scans(shape, 
                        tile_start_function=self.save_overview_images, ts_args=[g],
                        payload_function=self.pf_align_and_do_scan(datafile_group=g,**kwargs))
        
    def new_data_group(self,name="particleScans/scan%d",parent=None):
        if parent is None: parent=self.datafile
        n=0
        while name % n in parent: n+=1
        return parent.create_group(name % n)

    def new_dataset_name(self,g,name):
        n=0
        while name % n in g: n+=1
        return name % n

    def save_overview_images(self, datafile_group):
        self.csm.autofocus_iterate(np.arange(-5,5,0.5))
        """save an unmodified and a thresholded image, as a reference for scans"""
        time.sleep(1)
        self.csm.camera.update_latest_frame()
        img1 = datafile_group.create_dataset(self.new_dataset_name(datafile_group,"overview_image_%d"),data=self.csm.camera.color_image())
        img1.attrs.create("stage_position",self.csm.stage.position())
        img1.attrs.create("camera_centre_position",self.csm.camera_centre_position())
        img1.attrs.create("mapping_matrix_camera_to_sample",self.csm.camera_to_sample)
	img1.attrs.create("timestamp",datetime.datetime.now().isoformat())
        img2 = datafile_group.create_dataset(self.new_dataset_name(datafile_group,"overview_image_%d_thresholded"),data=self.threshold_image(self.denoise_image(
                self.csm.camera.gray_image())))
        img2.attrs.create("stage_position",self.csm.stage.position())
        img2.attrs.create("camera_centre_position",self.csm.camera_centre_position())
        for key, val in self.get(['median_filter_width','threshold_level']).iteritems():
                img2.attrs.create(key,val)
        img2.attrs.create("camera_to_sample_matrix",self.csm.camera_to_sample)
	img2.attrs.create("timestamp",datetime.datetime.now().isoformat())

    def pf_align_and_do_scan(self, dz=np.arange(-4,4,0.4), datafile_group=None):
        """Set up for a scan of all particles, then return a payload function.
        
        The "payload function" is suitable for the eponymous argument of 
        go_to_particles, and will autofocus, align particle to fibre, and take
        a Z stack.  NB the payload function "wraps up" the arguments neatly so
        we don't need to store things like the depth of the Z stack.
        """
        if datafile_group is None:
            datafile_group=self.new_data_group("particleScans/scan%d",self.datafile)

        def align_and_do_scan():
            # --- Initialize shutter positions by making sure they are all closed
            print("Matt: Initializing shutter positions.")
            light_shutter.close_shutter()
            raman.shutter.close_shutter()
            # --- Open white light shutter and close laser shutter
            print("Matt: Opening white light shutter.")
            light_shutter.open_shutter()
            # --- Fully open Shamrock slit
            print("Matt: Opening spectrometer slit.")
            raman.sham.SetSlit(2000)
            # --- Wait a bit for settings to be applied
            time.sleep(1)
            # --- Do Autofocus using OceanOptics Spectrum
            print("Matt: Doing autofocus.")
            self.csm.autofocus_iterate(np.arange(-2.5,2.5,0.5))
            self.aligner.spectrometer.integration_time = 300. #short integration time for alignment
            self.aligner.optimise_2D(tolerance=0.03,stepsize=0.2)
            self.aligner.spectrometer.integration_time = 1000. #long integration time for measurement
            # --- Initialize datafile
            print("Matt: Initializing datafile.")
            g = self.new_data_group("scan_%d",datafile_group)
            dset = g.create_dataset("scan",
                                    data=self.aligner.scan(dz))
            for key, val in self.aligner.spectrometer.get_metadata().iteritems():
                dset.attrs.create(key,val)
            dset.attrs.create("stage_position",self.csm.stage.position())
            dset.attrs.create("camera_centre_position",self.csm.camera_centre_position())
            dset.attrs.create("timestamp",datetime.datetime.now().isoformat())
            dset.attrs.create("dz",dz)
            # --- Close all shutters
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # --- Set Infinity3 camera settings --- TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            # --- Set Andor camera settings --- TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            # --- Take bias images on Infinity3, Andor, and OceanOptics
            # --- Infinity3
            print("Matt: Taking Infinity3 bias image.")            
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_Bias_Image",data=image)
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            # --- Andor 0 order
            print("Matt: Taking Andor 0 order bias image.")  
            raman.sham.GotoZeroOrder()
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Bias_0Order_int", data=image)
            g.create_dataset("Raman_Bias_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Andor spectrum
            print("Matt: Taking Andor spectrum bias image.")  
            raman.sham.SetWavelength(raman.centre_Wavelength)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Bias_Spectrum_int", data=image)
            g.create_dataset("Raman_Bias_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- OceanOptics
            print("Matt: Taking OceanOptics bias image.")  
            (oowl, oospec) = spectrometer.read()
            g.create_dataset("OOptics_Bias_Spectrum_int", data=oospec)
            g.create_dataset("OOptics_Bias_Spectrum_wl", data=oowl)
            # --- Turn on white light
            print("Matt: Turning white light back on.")  
            light_shutter.open_shutter()
            # --- Take Infinity3 image
            print("Matt: Taking first Infinity3 white light image.")
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_FirstWhiteLight_Image",data=image[image.shape[0]/2-50:image.shape[0]/2+50, image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            # --- Take white light spectrum (Ocean Optics) TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            # --- Take white light spectrum (Andor)
            print("Matt: Taking white light spectrum on Andor")
            raman.sham.SetWavelength(raman.centre_Wavelength)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_Spectrum_int", data=image)
            g.create_dataset("Raman_White_Light_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Take white light image (Andor)
            print("Matt: Taking Andor 0 order white light image.")  
            raman.sham.GotoZeroOrder()
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_0Order_int", data=image)
            g.create_dataset("Raman_White_Light_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Turn off white light
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # --- Turn on laser        
            print("Matt: Opening the laser shutter.")
            raman.shutter.open_shutter()
            # --- Set Infinity3 exposure/gain very low and image beam profile. Then restore old values
            oldExposure = cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._get_value()
            oldGain = cam.parameters[cam.parameters[0].list_names().index('GAIN')]._get_value()
            cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._set_value(0) #sometimes need to set to float(-inf)
            cam.parameters[cam.parameters[0].list_names().index('GAIN')]._set_value(0)
            print("Matt: Taking Infinity3 image of laser beam profile.")   
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_Laser_Beam_Image",data=image)
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._set_value(oldExposure)
            cam.parameters[cam.parameters[0].list_names().index('GAIN')]._set_value(oldGain)
            # --- Take image of laser zero-order (Andor)
            print("Matt: Taking Andor 0 order laser image.")  
            raman.sham.GotoZeroOrder()
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_0Order_int", data=image)
            g.create_dataset("Raman_Laser_0Order_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Take image of laser spectrum (Andor)
            print("Matt: Taking Andor spectrum laser image.")  
            raman.sham.SetWavelength(raman.centre_Wavelength)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_Spectrum_int", data=image)
            g.create_dataset("Raman_Laser_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc) 
            # --- Turn off laser         
            print("Matt: Closing the laser shutter.")
            raman.shutter.close_shutter()
            # --- Open the white light shutter
            print("Matt: Turning white light back on.")  
            light_shutter.open_shutter()
            # --- Take second Infinity3 white light image (to track drift)
            print("Matt: Taking second Infinity3 white light image.")
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_SecondWhiteLight_Image",data=image[image.shape[0]/2-50:image.shape[0]/2+50, image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            # --- Move stage slightly to take background. (MAKE THIS ACTUALLY LOOK FOR A SPOT WITH NO PARTICLES) TODO !!!!!!!!!!!!!!!!
            self.csm.stage.move_rel([1,0,0]) # Move the stage by one micron to a (hopefully) empty area
            # --- Take Infinity3 white light image (as evidence stage moved to a good location for background)
            print("Matt: Taking first Infinity3 white light image for background location.")
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_FirstBkgndWhiteLight_Image",data=image[image.shape[0]/2-50:image.shape[0]/2+50, image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            # --- Take white light spectrum at background location (Ocean Optics) TODO !!!!!!!!!!!!!!!!!!!!!!!!!!!!!!!            
            # --- Take white light spectrum at background location (Andor)
            print("Matt: Taking white light spectrum at background location on Andor")
            raman.sham.SetWavelength(raman.centre_Wavelength)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_Bkgnd_Spectrum_int", data=image)
            g.create_dataset("Raman_White_Light_Bkgnd_Spectrum_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Take white light image at background location (Andor)
            print("Matt: Taking Andor 0 order white light image at background location.")  
            raman.sham.GotoZeroOrder()
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_White_Light_Bkgnd_0Order_int", data=image)
            g.create_dataset("Raman_White_Light_0Order_Bkgnd_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Turn off white light
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            # --- Turn on laser          
            print("Matt: Opening the laser shutter.")
            raman.shutter.open_shutter()
            # --- Set Infinity3 exposure/gain very low and image beam profile. Then restore old values
            oldExposure = cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._get_value()
            oldGain = cam.parameters[cam.parameters[0].list_names().index('GAIN')]._get_value()
            cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._set_value(0) #sometimes need to set to float(-inf)
            cam.parameters[cam.parameters[0].list_names().index('GAIN')]._set_value(0)
            print("Matt: Taking Infinity3 image of laser beam profile at background location.")   
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_Laser_Beam_Image_atBkgndLoc",data=image)
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            cam.parameters[cam.parameters[0].list_names().index('EXPOSURE')]._set_value(oldExposure)
            cam.parameters[cam.parameters[0].list_names().index('GAIN')]._set_value(oldGain)
            # --- Take image of laser zero-order at background location (Andor)
            print("Matt: Taking Andor 0 order laser image at background location.")  
            raman.sham.GotoZeroOrder()
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_0Order_atBkgndLoc_int", data=image)
            g.create_dataset("Raman_Laser_0Order_atBkgndLoc_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc)
            # --- Take image of laser spectrum at background location (Andor)
            print("Matt: Taking Andor spectrum laser image at background location.")  
            raman.sham.SetWavelength(raman.centre_Wavelength)
            time.sleep(5)
            image = np.reshape(raman.take_bkg(),(-1, raman.sham.pixel_number))
            wavelengths = raman.GetWavelength()
            rint = g.create_dataset("Raman_Laser_Spectrum_atBkgndLoc_int", data=image)
            g.create_dataset("Raman_Laser_Spectrum_atBkgndLoc_wl", data=wavelengths)
            rint.attrs.create("Laser power",raman.laser_power)
            rint.attrs.create("Slit size",raman.slit_size)
            rint.attrs.create("Integration time",raman.Integration_time)
            rint.attrs.create("description",raman.scan_desc) 
            # --- Turn off laser         
            print("Matt: Closing the laser shutter.")
            raman.shutter.close_shutter()
            # --- Open the white light shutter
            print("Matt: Turning white light back on.")  
            light_shutter.open_shutter()
            # --- Take second Infinity3 white light image at background location (to track drift)
            print("Matt: Taking second Infinity3 white light image at background location.")
            #we're going to take a picture - best make sure we've waited a moment for the focus to return
            time.sleep(0.3)
            self.csm.camera.update_latest_frame() #take a frame and ignore (for freshness)
            image = self.csm.camera.color_image()
            img = g.create_dataset("Infinity3_SecondWhiteLight_atBkgndLoc_Image",data=image[image.shape[0]/2-50:image.shape[0]/2+50, image.shape[1]/2-50:image.shape[1]/2+50])
            img.attrs.create("stage_position",self.csm.stage.position())
            img.attrs.create("timestamp",datetime.datetime.now().isoformat())
            # --- Turn off all light sources
            print("Matt: Closing white light shutter.")
            light_shutter.close_shutter()
            raman.shutter.close_shutter()
            
            ### POTENTIALLY USEFUL: LASER FOCUS
            #"first bring the laser into focus"
            #here = stage.position()            
            #laser_focus = raman.AlignHeNe(dset) ### remove, keep, or create offset
            #print "Moving to HeNe Focus (%g)" % (laser_focus)
            #stage.move_rel([0,0,laser_focus])
            #time.sleep(1)            
            
            ###################################################################################################################
            print("Reached the end of Matt's code.")
            datafile_group.file.flush()   
            ###################################################################################################################
            
        return align_and_do_scan

    def plot_latest_scan(self):
        """plot the spectra from the most recent scan"""
        g = self.latest_scan_group
        for name, scangroup in g.iteritems():
            if re.match(r"scan_\d+",name):
                scan=scangroup['scan']
                spectrum = np.sum(scan,0)
                plt.plot(scan.attrs['wavelengths'],spectrum/scan.shape[0] - scan.attrs['background'])
        plt.show(block=False)
        
if __name__ == "__main__":
    #from nplab.instrument.spectrometer.seabreeze import OceanOpticsSpectrometer as OOSpec
    cam = camera.Camera(0)
    stage = prior_stage.ProScan("COM9")
    light_shutter = Uniblitz_Shutter.UniblitzShutter()
    light_shutter.open_connection()
    stage.query("SERVO 1") #set up stage to use servocontrol
    stage.query("UPR Z 500") #500 um per revolution = lab 6
    mapper = camera_stage_mapper.CameraStageMapper(cam,stage)
    mapper.frames_to_discard = 2 #the Infinity 3 is faster, so we need to throw away more frames
    mapper.autofocus()
    mapper.calibrate(7)
    spectrometer = seabreeze.OceanOpticsSpectrometer(0)
 #   spectrometer = OOSpec(0)
    spectrometer._set_tec_temperature(-21)
  #  spectrometer.set_tec_temperature(-21)
    aligner = spectrometer_aligner.SpectrometerAligner(spectrometer, stage)
    scanner = ParticleScanner(mapper, spectrometer, aligner, h5py.File("test_scans.hdf5",'a'))
    #aligner.spectrum_mask=[500<s and s<600 for s in spectrometer.wavelength]
    viewer = view_spectra.ScanViewer(scanner.datafile)
    ##################################################################################################################cam.edit_traits()
    spectrometer.edit_traits()
#    spectrometer.show_gui(blocking = False)
#    mapper.edit_traits()  
#    scanner.edit_traits()
    aligner.edit_traits()
    viewer.edit_traits()
    
    raman = Raman_class.Raman_spec()
    raman.CCD_cooling(-90)
    #set the Raman integration time
    #raman.exptime = 3 #in seconds
    
    raman.edit_traits()
    raman.cam.verbosity=False
    
    stop_tile = threading.Event()    

    scanner.edit_traits()    
    mapper.edit_traits()

    def tile_scans():
        while not stop_tile.is_set():
            scanner.take_zstacks_of_particles()
            scanner._scan_thread.join()
    #threading.Thread(target=tile_scans).start()
    #len([k for k in scanner.datafile['particleScans/scan1/'].keys() if "scan" in k])
    #you can count the number of particles with this code:
            #len(scanner.datafile['particleScans/scan0'].keys())
            #NB replace scan1 with the number of your latest scan
    
    def closeall():
        cam.close()
        stage.ser.close()
        seabreeze.shutdown_seabreeze()
        light_shutter.close_connection()
        raman.SystemShutdown()        
        
    
    
