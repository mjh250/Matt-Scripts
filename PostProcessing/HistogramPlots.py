# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 18:14:56 2018

@author: mjh250
"""

import numpy as np
import matplotlib.pyplot as plt


if __name__ == '__main__':
    # All vectors show data from measurements in order:
    # 2018/02/17, 2018/02/16, 2017/12/14
    sample_ages = [0, 3, 8]  # Age of each saple at measurement time in days
    total_scans = [255, 152, 133]  # Total particles scanned
    total_particles = np.asarray([180, 82, 84])  # Total particles scanned minus duds
    rings = np.asarray([17, 16, 17])  # Total rings
    faint = np.asarray([0, 3, 0])
    spots = np.asarray([160, 52, 47])
    asymmetric = np.asarray([3, 11, 20])

    rings_percent = (100.0*rings/total_particles)  # Total rings
    faint_percent = (100.0*faint/total_particles)
    spots_percent = (100.0*spots/total_particles)
    asymmetric_percent = (100.0*asymmetric/total_particles)

    # create plot
    n_groups = 3
    fig, ax = plt.subplots()
    index = np.arange(n_groups)
    bar_width = 0.2
    opacity = 0.8

    rects1 = plt.bar(index, spots_percent, bar_width,
                     alpha=opacity,
                     color='#e9724d',
                     label='Spots')

    rects2 = plt.bar(index + bar_width, rings_percent, bar_width,
                     alpha=opacity,
                     color='#d6d727',
                     label='Rings')

    rects3 = plt.bar(index + 2*bar_width, asymmetric_percent, bar_width,
                     alpha=opacity,
                     color='#92cad1',
                     label='Asymmetric')

    rects4 = plt.bar(index + 3*bar_width, faint_percent, bar_width,
                     alpha=opacity,
                     color='#868686',
                     label='Faint')

    plt.xlabel('Sample Age', fontsize=30)
    plt.ylabel('Percentage of total occurences', fontsize=30)
    plt.title('Relative abundance of Raman scattering profiles', fontsize=30)
    plt.xticks(index + 1.5*bar_width,
               ('0 Days old', '3 Days old', '8 Days old'),
               fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=30)

    # plt.tight_layout()
    plt.grid()
    plt.show()
    print('All done!')
