# -*- coding: utf-8 -*-
"""
Created on Wed Mar 07 15:43:54 2018

@author: mjh250
"""

import shutil
import time
import numpy as np
import matplotlib.pyplot as plt
import matplotlib.ticker as mtick
from ParticleScanAnalysis.ParticleScanAnalysis import ParticleScanAnalysis

if __name__ == '__main__':
    print("Initializing...")
    filepaths = ["C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\2017_12_14\\BackgroundSubtracted.hdf5"]
    processed_filepath = "C:\\Users\\mjh250\\Documents\\Local mjh250\\mjh250\\2017_12_14\\Counts"
    if len(filepaths) > 1:
        processed_filepaths = [processed_filepath + str(i)
                               for i in range(len(filepaths))]
    else:
        processed_filepaths = [processed_filepath]

    for fileNum in range(len(filepaths)):
        filepath = filepaths[fileNum]
        processed_filepath = processed_filepaths[fileNum]
        try:
            shutil.rmtree(processed_filepath)
            time.sleep(1)  # Wait for permissions to be freed up
            shutil.os.mkdir(processed_filepath)
        except OSError:
            shutil.os.mkdir(processed_filepath)

        scan_analyzer = ParticleScanAnalysis()
        data = scan_analyzer.getScanDataFile(filepath)

        print("Preparing list of particles...")
        scan_number_list = []
        particle_number_list = []
        for sname in data:
            if sname.startswith('scan'):
                scan_number_list.append(int(sname[4:]))
                particle_num_sublist = []
                for pname in data[sname]:
                    if pname.startswith('particle'):
                        particle_num_sublist.append(int(pname[8:]))
                particle_number_list.append(particle_num_sublist)

    for sublist in particle_number_list:
        sublist.sort()
        scan_number_list.sort()

        intCountsRaman = []
        intCountsDF = []
        avgMaxCountsRaman = []
        avgMaxCountsDF = []

        # OVERRIDE IF NECESSARY
        # scan_number_list = [3]

        # --- Iterate through particle scans and process files
        for n, scan_number in enumerate(scan_number_list):
            intCountsDFSublist = []
            intCountsRamanSublist = []
            avgMaxCountsRamanSublist = []
            avgMaxCountsDFSublist = []
            for particle_number in particle_number_list[n]:
                # --- Print name of file being analyzed to track progress
                print('Processing scan%s/particle%s' % (scan_number,
                                                        particle_number))

                # --- GET SCAN DATA ---
                datafile = scan_analyzer.getProcessedScanDataSet(
                                        data, scan_number, particle_number)
                plt.ioff()  # No interactive plots. To suppress plot window

                # --- RAMAN LASER 0 ORDER IMAGE ---
                data_image = np.array(
                        datafile['Raman_Laser_0Order_Processed_Image'])
                intCountsDFSublist.append(datafile['Raman_Laser_0Order_Processed_Image'].attrs['integral'])
                avgMaxCountsRamanSublist.append(datafile['Raman_Laser_0Order_Processed_Image'].attrs['average of 10 maxima'])

                # --- RAMAN DF 0 ORDER IMAGE ---
                data_image = np.array(
                    datafile['Raman_White_Light_0Order_Processed_Image'])
                intCountsRamanSublist.append(datafile['Raman_White_Light_0Order_Processed_Image'].attrs['integral'])
                avgMaxCountsDFSublist.append(datafile['Raman_White_Light_0Order_Processed_Image'].attrs['average of 10 maxima'])
            intCountsDF.append(intCountsDFSublist)
            intCountsRaman.append(intCountsRamanSublist)
            avgMaxCountsRaman.append(avgMaxCountsRamanSublist)
            avgMaxCountsDF.append(avgMaxCountsDFSublist)

        data.close()
        print("Finished processing " + filepath.split('/')[-1] + ' !')

    #plt.plot(intCountsDF[0])
    #plt.plot(intCountsRaman[0])
    #plt.show()

    # indices of particles that are:
    total_particle_count = 131
    rings = [7, 18, 22, 24, 31, 39, 46, 51, 56, 69, 78, 80, 81, 82, 86, 117, 125]
    dims = []
    asymmetrics = [10, 21, 28, 29, 35, 45, 49, 57, 62, 66, 68, 75, 77, 84, 98, 99, 105, 106, 116, 120]
    junks = [3, 5, 13, 15, 16, 17, 25, 36, 48, 50, 70, 71, 74, 85, 87, 88, 93, 94, 104, 109, 118, 122, 123, 124, 131,
    	1, 6, 20, 23, 32, 34, 38, 47, 52, 53, 54, 55, 59, 60, 83, 89, 110, 112, 113, 114, 127, 128, 129]
    spots = [x for x in range(0, total_particle_count) if x not in rings+dims+asymmetrics+junks]

    ringCountsDF = [intCountsDF[0][x] for x in rings]
    dimCountsDF = [intCountsDF[0][x] for x in dims]
    asymmetricCountsDF = [intCountsDF[0][x] for x in asymmetrics]
    junkCountsDF = [intCountsDF[0][x] for x in junks]
    spotCountsDF = [intCountsDF[0][x] for x in spots]

    ringCountsRaman = [intCountsRaman[0][x] for x in rings]
    dimCountsRaman = [intCountsRaman[0][x] for x in dims]
    asymmetricCountsRaman = [intCountsRaman[0][x] for x in asymmetrics]
    junkCountsRaman = [intCountsRaman[0][x] for x in junks]
    spotCountsRaman = [intCountsRaman[0][x] for x in spots]

    ringMaxCountsDF = [avgMaxCountsDF[0][x] for x in rings]
    dimMaxCountsDF = [avgMaxCountsDF[0][x] for x in dims]
    asymmetricMaxCountsDF = [avgMaxCountsDF[0][x] for x in asymmetrics]
    junkMaxCountsDF = [avgMaxCountsDF[0][x] for x in junks]
    spotMaxCountsDF = [avgMaxCountsDF[0][x] for x in spots]

    ringMaxCountsRaman = [avgMaxCountsRaman[0][x] for x in rings]
    dimMaxCountsRaman = [avgMaxCountsRaman[0][x] for x in dims]
    asymmetricMaxCountsRaman = [avgMaxCountsRaman[0][x] for x in asymmetrics]
    junkMaxCountsRaman = [avgMaxCountsRaman[0][x] for x in junks]
    spotMaxCountsRaman = [avgMaxCountsRaman[0][x] for x in spots]

    # create plot
    n_groups = 2
    fig, ax = plt.subplots()
    plt.figure(tight_layout=True)
    index = np.arange(n_groups)
    bar_width = 0.2
    opacity = 0.8

    rects1 = plt.bar(index, [np.mean(ringMaxCountsDF), np.mean(ringMaxCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#e9724d',
                     label='Rings')

    rects2 = plt.bar(index + bar_width, [np.mean(spotMaxCountsDF), np.mean(spotMaxCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#d6d727',
                     label='Spots')

    rects3 = plt.bar(index + 2*bar_width, [np.mean(asymmetricMaxCountsDF), np.mean(asymmetricMaxCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#92cad1',
                     label='Asymmetrics')

    rects4 = plt.bar(index + 3*bar_width, [np.mean(dimMaxCountsDF), np.mean(dimMaxCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#868686',
                     label='Dims')

    plt.xlabel('Measurement', fontsize=30)
    plt.ylabel('Mean counts of 10 maximal pixels', fontsize=30)
    plt.title('Relative abundance of Raman scattering profiles', fontsize=30)
    plt.xticks(index + 1.5*bar_width,
               ('Dark field', 'Raman'),
               fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=30)

    # plt.tight_layout()
    plt.grid()
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.draw()
    plt.figure(tight_layout=True)
    ax = plt.gca()
    
    rects1 = plt.bar(index, [np.mean(ringCountsDF), np.mean(ringCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#e9724d',
                     label='Rings')

    rects2 = plt.bar(index + bar_width, [np.mean(spotCountsDF), np.mean(spotCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#d6d727',
                     label='Spots')

    rects3 = plt.bar(index + 2*bar_width, [np.mean(asymmetricCountsDF), np.mean(asymmetricCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#92cad1',
                     label='Asymmetrics')

    rects4 = plt.bar(index + 3*bar_width, [np.mean(dimCountsDF), np.mean(dimCountsRaman)], bar_width,
                     alpha=opacity,
                     color='#868686',
                     label='Dims')

    plt.xlabel('Measurement', fontsize=30)
    plt.ylabel('Integrated Counts', fontsize=30)
    plt.title('Relative abundance of Raman scattering profiles', fontsize=30)
    plt.xticks(index + 1.5*bar_width,
               ('Dark field', 'Raman'),
               fontsize=30)
    plt.yticks(fontsize=30)
    plt.legend(fontsize=30, loc=3)

    # plt.tight_layout()
    plt.grid()
    ax.yaxis.set_major_formatter(mtick.FormatStrFormatter('%.2e'))
    plt.draw()
    plt.show()
    
    print('All done!')
