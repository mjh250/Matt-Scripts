# -*- coding: utf-8 -*-
"""
Created on Wed Feb 21 18:14:56 2018

@author: mjh250
"""

import numpy as np
import matplotlib.pyplot as plt

# CB[7] + MB not filtered
# rings : 17
# dims  : 0
# spots : 160
# asymmetric : 3
# Total : 180

# CB[7] only:
# rings : 0
# dims  : 127
# spots : 24
# asymmetric : 2
# Total : 153

# BPT :
# rings : 1
# dims  : 19
# spots : 41
# asymmetric : 8
# Total : 69

# CB[7] + MB filtered
# rings : 12
# dims  : 0
# spots : 124
# asymmetric : 9
# Total : 145

# 100nm + CB[7] + MB:
# rings : 5
# dims  : 1
# spots : 86
# asymmetric : 7
# Total : 99

if __name__ == '__main__':
    # All vectors show data from measurements in order:
    # 2018/02/17, 2018/02/18
    sample_ages = [0, 0]  # Age of each saple at measurement time in days
    total_scans = [255, 143]  # Total particles scanned
    total_particles = np.asarray([180,145])  # Total particles scanned minus duds
    rings = np.asarray([17,12])  # Total rings
    faint = np.asarray([0,0])
    spots = np.asarray([160,124])
    asymmetric = np.asarray([3,9])

    rings_percent = (100.0*rings/total_particles)  # Total rings
    faint_percent = (100.0*faint/total_particles)
    spots_percent = (100.0*spots/total_particles)
    asymmetric_percent = (100.0*asymmetric/total_particles)

    # create plot
    n_groups = 2
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

    plt.xlabel('Sample', fontsize=30)
    plt.ylabel('Percentage of total occurences', fontsize=30)
    plt.title('Relative abundance of Raman scattering profiles', fontsize=30)
    plt.xticks(index + 1.5*bar_width,
               ('CB[7]+MB+80nmAu filtered', 'CB[7]+MB+80nmAu unfiltered'),
               fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=30)

    # plt.tight_layout()
    plt.grid()
    plt.show()
    print('All done!')
