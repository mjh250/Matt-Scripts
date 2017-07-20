# -*- coding: utf-8 -*-
"""
@author: Felix Benz (fb400)
"""

import visa, time

class UniblitzShutter():
    def __init__(self):
        self.visa_address = str('ASRL7::INSTR')
        self.shutter_state = True
    
    def open_connection(self):
        rm = visa.ResourceManager()
        self.device = rm.open_resource(self.visa_address, baud_rate=9600, read_termination='\r', timeout=1000)   

    def open_shutter(self):
        if self.shutter_state:
            print "Shutter already open"
        else:
            self.device.write('@')
            self.shutter_state = True

    def close_shutter(self):
        if self.shutter_state:
            self.device.write('A')
            self.shutter_state = False
        else:
            print "Shutter already closed"
    
    def close_connection(self):
        self.open_shutter()
        self.device.close()
        
if __name__=='__main__':
    Shutter = UniblitzShutter()
