# -*- coding: utf-8 -*-
"""
Created on Wed Apr 08 13:20:35 2015

@author: felixbenz
"""

from nplab.instrument.shutter.thorlabs_sc10 import ThorLabsSC10
from nplab.instrument.camera.Andor import Andor
from nplab.instrument.spectrometer.shamrock import Shamrock
from Powermeter import ThorlabsPM100
import time
import visa
import sys
import signal
import numpy as np
import h5py
import scipy.ndimage.filters

import traits
from traits.api import HasTraits, Property, Instance, Float, Int, String, Button, Bool, on_trait_change
import traitsui
from traitsui.api import View, Item, Group, HGroup, VGroup
import threading
import chaco
from chaco.api import ArrayPlotData, Plot
from chaco.chaco_plot_editor import ChacoPlotItem
from enable.component_editor import ComponentEditor


class Raman_spec(HasTraits):
    slit_size = Float(10)
    centre_wavelength = Float(600)
    laser_power = Float(0)
    CCD_temperature = Int(-65)
    Integration_time = Float(0.1) #s
    
    Raman_spectrum = traits.trait_numeric.Array(shape=(None))
    Raman_wavelengths = traits.trait_numeric.Array(shape=(None))
    CoolerStatus = Bool(False)
    
    scan_number = Int(0)
    scan_desc = String('')    
    
    go_to_zero = Button()
    set_wavelength = Button()
    set_slit = Button()
    reset_slit = Button()
    
    read_laser_power = Button()
    laser_shutter = Button()
    
    system_shutdown = Button()
    next_scan = Button()
    
    get_temperature = Button()
    set_temperature = Button()
    wait_for_temperature_to_stabilise = Button()    
    
    acquire_data = Button()
    set_integration_time = Button()
    save_data = Button()    
    
    Raman_plot = Property(trait=Instance(Plot))

    traits_view = View(VGroup(
                        HGroup(
                            HGroup(
                                Item(name="slit_size",label = "Slit Size (um)"),
                                Item(name="set_slit",show_label=False),
                                Item(name="reset_slit",show_label=False)
                            ),
                            HGroup(
                                Item(name="centre_wavelength",label = "Centre Wavelength (nm)"),
                                Item(name="set_wavelength",show_label=False),
                                Item(name="go_to_zero",show_label=False)
                            ),label="Shamrock 303i",show_border=True),
                        HGroup(
                            Item(name="laser_power",label = "Laser power (mW)",style="readonly"),
                            Item(name="read_laser_power",show_label=False),
                            Item(name="laser_shutter",show_label=False),
                            Item(name="system_shutdown",show_label=False),
                        show_border=True),
                        HGroup(
                            Item(name="scan_desc",label="Scan description"),
                            Item(name="scan_number",label="Scan number",width=10,style="readonly"),
                            Item(name="next_scan",width=10,show_label=False)                        
                        ,show_border=True),
                        VGroup(
                            HGroup(
                                VGroup(
                                    Item("CCD_temperature",label="CCD temperature"),
                                    HGroup(Item("get_temperature",show_label=False),Item("set_temperature",show_label=False)),
                                    Item("CoolerStatus",label="CCD cooling"),
                                    Item("wait_for_temperature_to_stabilise",show_label=False)
                                ),
                                VGroup(
                                    HGroup(Item("Integration_time",label="Integration time (s)"),Item("set_integration_time",show_label=False)),
                                    Item("acquire_data",show_label=False),
                                    Item("save_data",show_label=False)
                                ),
                            ),
                            Item(name="Raman_plot",editor=ComponentEditor(),show_label=False)
                            ,label="Andor Newton 970",show_border=True)                        
                    ),resizable=True,width=800,height=800,title="Raman Experiment")
                        
    def __init__(self):
        #set the shutter, can be opened by calling shutter.toggle()"
        self.shutter = ThorLabsSC10('COM1')

        #Set the powermeter"
        rm = visa.ResourceManager()
        inst = rm.open_resource('USB0::0x1313::0x8072::P2004571::INSTR')
        self.powermeter = ThorlabsPM100(inst=inst)
        self.laserpower = None
    
        #Some basic command for the Andor"
        self.cam = Andor()
        self.EMCCDGain = 0
        self.PreAmpGain = 0
        self.exptime = 3
        #0: EM Amplifier
        #1: Conventional Amplifier
        self.preamp = 1
        #0: 3 MH
        #1: 1 MHz
        #2: 50 kHz
        self.readout_rate = 0
        
        #Acquisition settings"
        #Full Image"
        #self.centre_line = 100
        #self.readout_height = 10
        self.cam.SetParameter('ReadMode', 4)
        #self.cam.SetParameter('SingleTrack', self.centre_line,self.readout_height)
        
        self.cam.SetParameter('TriggerMode', 0)
        self.cam.SetParameter('Shutter', 1,0,30,30)
        self.cam.SetParameter('Exposure', self.exptime)
        self.cam.SetParameter('CoolerMode', 0)
        
        self.cam.SetParameter('OutAmp', self.preamp)
        self.cam.SetParameter('PreAmpGain', self.PreAmpGain)
        #self.cam.SetParameter('EMCCDGain', EMCCDGain)

        self.cam.SetParameter('NumHSSpeed', self.preamp, self.readout_rate)
        
        #some basic commands for the shamrock
        self.sham = Shamrock()
        self.centre_Wavelength = 640 #640 for viewing antistokes, 696 otherwise
        self.slit_width = 100
        
        self.sham.SetPixelWidth(16)
        self.sham.SetNumberPixels(1600)
        self.sham.SetWavelength(self.centre_Wavelength)
        time.sleep(5)
        self.sham.SetSlit(self.slit_width)
        time.sleep(5)
        
        self.slit_size = self.sham.GetSlit()
        self.centre_wavelength = float(str(self.sham.GetWavelength())[:5])
        self.laser_power = float(str(self.read_power()*1E3)[:5])
        self.cam.GetParameter('CurrentTemperature')
        self.CCD_temperature=self.cam.CurrentTemperature
        self.Integration_time = self.exptime
        
        #self.Raman_spectrum = np.zeros(1600)
        #self.Raman_wavelengths = np.arange(1600)
        
        #Some operations to prepare the HDF5 file
        self.file = h5py.File("test_file.hdf5",mode="w") 
        self.output_group = self.file.require_group('RamanDFExperiment')
        self.save_group = self.output_group.require_group('scan%d' %self.scan_number)
        
    def read_power(self):
            objective_correction = 1.1809E+02
            laserpower = objective_correction*self.powermeter.read
            return laserpower

    "Cool the CCD"
    def CCD_cooling(self,Tset):
        self.cam.SetParameter('CoolerMode', 0)
        self.cam.CoolerON()
        self.cam.SetParameter('SetTemperature', Tset)
        self.cam.GetParameter('CurrentTemperature')
        self.CCD_temperature=self.cam.CurrentTemperature
        self.CCD_temp_setpoint = Tset
         
        self.cooling_thread = threading.Thread(target=self.CCD_cooling_function)
        self.cooling_thread.start()
    
    def CCD_cooling_function(self):
        """this function should only EVER be executed by _wait_for_temperature_to_stabilise_fired."""
        while self.cam.GetParameter('CurrentTemperature') is not 'DRV_TEMP_STABILIZED':
            self.CCD_temperature=self.cam.CurrentTemperature
            print "Temperature is: %g (Tset: %g)" % (self.cam.CurrentTemperature, self.CCD_temp_setpoint)
            sys.stdout.flush()
            time.sleep(30)
        return "Cooling finished!"    
    
    def GetWavelength(self):
        self.sham.GetCalibration()
        wavelengths = self.sham.wl_calibration
        return wavelengths
        
    def take_bkg(self):
        return self.cam.capture()[0]
        
    def take_signal(self):
        self.shutter.open_shutter()
        time.sleep(0.1)
        self.laserpower = self.read_power()
        self.shutter.close_shutter()
        return self.cam.capture()[0]
        
    def take_fast_kinetic(self):
        
        #fast raman
        self.cam.SetParameter('AcquisitionMode', 3)
        self.cam.SetParameter('ReadMode', 0)
        self.readout_rate = 2 #50kHz
        self.cam.SetParameter('NumHSSpeed', self.preamp, self.readout_rate)
        self.exptime = 0.01
        self.cam.SetParameter('Exposure', self.exptime) #set exposure time    
        self.cycles = 1000 #number of cycles
        self.cam.SetParameter('NKin', self.cycles) 
        
            
        #get wavelength calibration
        self.sham.GetCalibration()
        self.Raman_wavelengths = self.sham.wl_calibration
            
        #open camera shutter
        self.cam.SetParameter('Shutter', 1,5,30,30) # 5 is open shutter for all series
        
        time.sleep(0.1)
            
        self.cam.GetParameter('AcquisitionTimings')
        self.cycle_time = self.cam.kinetic
        
        self.shutter.open_shutter()
        time.sleep(0.1)
        
        objective_correction = 0.045322
        laserpower = objective_correction*self.read_power()
        self.laserpower = float(str(laserpower*1E3)[:5])
        
        #create datasets for time, power, and spectra
        spectra = np.array([])
        self.times = np.arange(0.0,float(self.cycle_time*self.cycles),float(self.cycle_time))
            
        kinetic_data = self.cam.capture()[0]
        self.shutter.close_shutter()
        self.kinetic_scan_data = np.array((kinetic_data),dtype = 'float64')
        
    
    def _go_to_zero_fired(self):
        self.sham.GotoZeroOrder()
        time.sleep(5)
        self.centre_wavelength = self.sham.GetWavelength()    
        
    def _set_wavelength_fired(self):
        self.sham.SetWavelength(self.centre_wavelength)
        
    def _set_slit_fired(self):
        self.sham.SetSlit(self.slit_size)
        
    def _reset_slit_fired(self):
        self.sham.SlitReset()
        time.sleep(5)
        self.slit_size = self.sham.GetSlit()
    
    def _read_laser_power_fired(self):
        self.laser_power = float(str(self.read_power()*1E3)[:5])
        
    def _laser_shutter_fired(self):
        self.shutter.toggle()
        
    def _system_shutdown_fired(self):
        self.SystemShutdown()
    
    def SystemShutdown(self):
        self.sham.Close()
        self.cam.__del__()
        self.file.close()
    
    def _next_scan_fired(self):
        self.scan_number=self.scan_number+1
        self.save_group = self.output_group.require_group('scan%d' %self.scan_number)        
        
    def _set_temperature_fired(self):
        self.cam.SetParameter('SetTemperature', self.CCD_temperature)
    
    def _get_temperature_fired(self):
        self.cam.GetParameter('CurrentTemperature')
        self.CCD_temperature=self.cam.CurrentTemperature
        
    def _set_integration_time_fired(self):
        self.cam.SetParameter('Exposure', self.Integration_time)
    
    def _acquire_data_fired(self):
        self.Raman_spectrum = self.cam.capture()[0]
        self.Raman_wavelengths = self.GetWavelength()
        self._get_Raman_plot()
        
    def _CoolerStatus_changed(self):
        if(self.CoolerStatus is True):
            self.cam.CoolerON()
        elif(self.CoolerStatus is False):
            self.cam.CoolerOFF()
    
    def _wait_for_temperature_to_stabilise_fired(self):
        if(self.CoolerStatus is False):
            self.cam.CoolerON()
            self.CoolerStatus=True
        self.cam.GetParameter('CurrentTemperature')
        self.CCD_temperature=self.cam.temperature
        
        #self._cooling_stop_event = threading.Event()
        self.cooling_thread = threading.Thread(target=self.CCD_cooling_function)
        self.cooling_thread.start()    
                
    @traits.api.cached_property # only draw the graph the first time it's needed
    def _get_Raman_plot(self):
        """make a nice plot of spectrum vs wavelength"""
        self._plot_data = ArrayPlotData(wavelengths=self.Raman_wavelengths,spectrum=self.Raman_spectrum)
        p = Plot(self._plot_data)
        p.plot(("wavelengths","spectrum"),type='line',color='blue')
        return p

    def _get_Raman_wavelengths(self):
        self.sham.GetCalibration()
        wavelengths = self.sham.wl_calibration
        return wavelengths
        
    def _Raman_spectrum_changed(self):
        if hasattr(self, "_plot_data"):
            self._plot_data.set_data("spectrum",self.Raman_spectrum)

    def _Raman_wavelengths_changed(self):
        if hasattr(self, "_plot_data"):
            self._plot_data.set_data("wavelengths",self.Raman_wavelengths)

    def _save_data_fired(self):
        self.save_group.create_dataset("Raman_int", data=self.Raman_spectrum)
        self.save_group.create_dataset("Raman_wl", data=self.Raman_wavelengths)
        self.save_group['Raman_int'].attrs.create("Laser power",self.laser_power)
        self.save_group['Raman_int'].attrs.create("Slit size",self.slit_size)
        self.save_group['Raman_int'].attrs.create("Integration time",self.Integration_time)
        self.save_group['Raman_int'].attrs.create("description",self.scan_desc)
    
    def AlignHeNe(self,zscan):
        peak_z = self.gaussians_fit_z_scan(zscan)
        
        if abs(peak_z[315]) < 6: #in case the fit doesn't work the function isn't moving too far off
            self.HeNeFocus = peak_z[315]
        else:
            self.HeNeFocus = 0
            
        return self.HeNeFocus
    
    
    def gaussians_fit_z_scan(self,z_scan, z=None, background=None, threshold=0.2, smoothing_width=1.5):
        """Fit a Gaussian to the spectrometer counts vs Z curve for each spectrometer pixel.
        
        Parameters:
            z_scan: 2D array (will be converted to numpy.ndarray) of spectra vs
                Z, where Z is the first index and wavelength is the second.
            z: 1D array of Z values for each row of z_scan.  If z_scan is an HDF5
                dataset with attribute 'dz', this will be used by default.
            background: 1D array containing a dark spectrum which is subtracted
                from the z_scan.  Defaults to HDF5 attribute 'background' as above.
                
        Return value: spectrum, peak_z, standard_deviation, background, r
            spectrum: peak height for the fitted Gaussian at each wavelength
            peak_z: Z value at which the peak height occurs
            standard_deviation: width of the Gaussian
            background: revised version of the supplied (or extracted) background,
                if the Gaussian fit suggests there is a constant background.
            r: correlation coefficient between fitted and actual curve, which can
                be a useful diagnostic as to how well the fit worked.
        
        N.B. You can re-run the fit, specifying the returned background - this may
        improve accuracy, though the effect seems minimal now that it uses a
        threshold to find the peak.  If the SNR of your data is particularly low,
        or you might have peaks very close to the edge, it may make sense to
        increase the threshold.
        """
        if z is None: #try to find the Z values of the different spectra from metadata if not specified
            try:
                z=z_scan.attrs.get('dz') #this should work for datasets in an HDF5 file
            except:
                print "Warning: no valid z positions were found, using indices instead"
                z=range(z_scan.shape[0]) #fall back to indices
        if background is None: #similarly, get a head-start on the background from metadata
            try:
                background=z_scan.attrs.get('background') #this should work for datasets in an HDF5 file
            except:
                z=np.zeros(z_scan.shape[1]) #fall back to zero.  NB it will take a *really* long time to converge like this...
        z_scan_array = np.array(z_scan, dtype=np.float) - background #first prepare background-subtracted, thresholded arrays
        scan_s = scipy.ndimage.filters.gaussian_filter1d(z_scan_array, smoothing_width, axis=0, mode='nearest')
        scan_s[scan_s<0.] = 0. #smooth and remove negative points so that we don't get NaNs
        scan_t = (scan_s - scan_s.min(axis=0))/(scan_s.max(axis=0)-scan_s.min(axis=0)) - threshold
        scan_t[scan_t<0.] = 0. #background-subtracted, thresholded version of z_scan_array - used for finding the peak
        #Now the fitting starts: find the centroid and width of the distribution in Z (the threshold helps avoid bias from background):
        peak_z = np.sum(scan_t * z[:,np.newaxis], axis=0)/np.sum(scan_t,axis=0)#the z[np.newaxis,:] is "numpy broadcasting", google it!
        standard_deviation = np.sqrt(np.sum(scan_s * (z[:,np.newaxis] - peak_z)**2, axis=0)/np.sum(scan_s,axis=0)) 
        #Next, construct a Gaussian in Z with the right parameters and use linear least squares to fit offset/peak to the data
        gaussians = np.exp(-(z[:,np.newaxis]-peak_z[np.newaxis,:])**2/(2*standard_deviation**2))
        var_x = np.var(gaussians,axis=0)
        mean_x = np.mean(gaussians,axis=0)
        var_y = np.var(z_scan_array,axis=0)
        mean_y = np.mean(z_scan_array,axis=0)
        covariance = np.mean((gaussians-mean_x) * (z_scan_array - mean_y),axis=0)
        slopes = covariance/var_x
        intercepts = mean_y - slopes * mean_x
        peak_height = slopes
        r = covariance/np.sqrt(var_x*var_y)
        return peak_z    
    
if __name__ == "__main__":
    Raman = Raman_spec()
    Raman.edit_traits()