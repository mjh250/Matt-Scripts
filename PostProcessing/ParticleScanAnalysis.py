# -*- coding: utf-8 -*-
"""
Created on Thu Aug 03 16:33:43 2017

@author: mjh250
"""

import h5py


def getScanData(filepath, scan_number):
    file = h5py.File(filepath, 'r')
    dataset = file['/particleScans/scan%s' % scan_number]
    return dataset

if __name__ == '__main__':
    filepath = 'test_scans.hdf5'
    data = getScanData(filepath, 0)
    print(data['scan_0'].keys())
