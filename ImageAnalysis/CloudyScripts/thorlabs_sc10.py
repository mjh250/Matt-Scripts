# -*- coding: utf-8 -*-
"""
Created on Thu Jun 12 18:00:10 2014

@author: alansanders
"""

from traits.api import HasTraits, Str, Button, Bool
from traitsui.api import View, Item, Group, Spring
from serial import Serial

class SC10Error(Exception):
    def __init__(self, msg):
        self.value = msg
        print "SC10 error: {:s}".format(msg)

class SC10(HasTraits):

    comport = Str('COM1')
    toggle_button = Button('Laser Shutter')
    state = Bool(False)
    
    traits_view = View(
        Item('toggle_button', show_label=False),  
        resizable=True, title='Thorlabs SC10'  
    )    
    
    def __init__(self, dummy=False):
        if dummy: self.ser = None
        else: self.ser = Serial(self.comport, baudrate=9600, bytesize=8, parity='N', stopbits=1, timeout=1)        
        self.init_cmds()
        
    def init_cmds(self):
        self.set_params = {
                        'toggle' : ('ens', None),
                        }  
        self.get_params = {
                        'state' : ('ens?', None),
                        }
    
    def get_attr(self, attr):
        if attr not in self.get_params:
            raise SC10Error('invalid attribute')
        self.ser.write(self.get_params[attr][0] + '\r')
        s = ''
        while self.ser.inWaiting() > 0:
            s += self.ser.readline().strip()
        return s
            
    def set_attr(self, attr, value):
        if attr not in self.set_params:
            raise SC10Error('invalid attribute')
        self.ser.write(self.set_params[attr][0]+str(value)+'\r')    
    
    def _toggle_button_fired(self):
        self.set_attr('toggle', '')
        
    def trigger(self):
        self.set_attr('toggle', '')
        
if __name__ == '__main__':
    shutter = SC10()
    shutter.configure_traits()