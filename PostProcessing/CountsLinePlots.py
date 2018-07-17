# -*- coding: utf-8 -*-
"""
Created on Thu Apr 26 10:52:13 2018

@author: mjh250
"""

import matplotlib.pyplot as plt
import numpy as np


if __name__ == "__main__":
    filepath8 = 'C://Users//mjh250//Documents//Local mjh250//mjh250//2017_12_14//Figures//'
    print("Grabbing data from " + filepath8)
    IntCounts8_file = open(filepath8+"IntegratedCounts.txt", "r")
    IntCounts8 = IntCounts8_file.read().split()
    IntCounts8 = [float(i) for i in IntCounts8]
    MaxCounts8_file = open(filepath8+"MaxCounts.txt", "r")
    MaxCounts8 = MaxCounts8_file.read().split()
    MaxCounts8 = [float(i) for i in MaxCounts8]
    
    filepath3 = 'C://Users//mjh250//Documents//Local mjh250//mjh250//2018_02_16//Figures//'
    print("Grabbing data from " + filepath3)
    IntCounts3_file = open(filepath3+"IntegratedCounts.txt", "r")
    IntCounts3 = IntCounts3_file.read().split()
    IntCounts3 = [float(i) for i in IntCounts3]
    MaxCounts3_file = open(filepath3+"MaxCounts.txt", "r")
    MaxCounts3 = MaxCounts3_file.read().split()
    MaxCounts3 = [float(i) for i in MaxCounts3]
    
    filepath0 = 'C://Users//mjh250//Documents//Local mjh250//mjh250//2018_02_17//Figures//'
    print("Grabbing data from " + filepath0)
    IntCounts0_file = open(filepath0+"IntegratedCounts.txt", "r")
    IntCounts0 = IntCounts0_file.read().split()
    IntCounts0 = [float(i) for i in IntCounts0]
    MaxCounts0_file = open(filepath0+"MaxCounts.txt", "r")
    MaxCounts0 = MaxCounts0_file.read().split()
    MaxCounts0 = [float(i) for i in MaxCounts0]
    
    DFIntCounts8 = IntCounts8[:len(IntCounts8)/2]
    DFMaxCounts8 = MaxCounts8[:len(MaxCounts8)/2]
    RamanIntCounts8 = IntCounts8[len(IntCounts8)/2:]
    RamanMaxCounts8 = MaxCounts8[len(MaxCounts8)/2:]
    DFIntCounts3 = IntCounts3[:len(IntCounts3)/2]
    DFMaxCounts3 = MaxCounts3[:len(MaxCounts3)/2]
    RamanIntCounts3 = IntCounts3[len(IntCounts3)/2:]
    RamanMaxCounts3 = MaxCounts3[len(MaxCounts3)/2:]
    DFIntCounts0 = IntCounts0[:len(IntCounts0)/2]
    DFMaxCounts0 = MaxCounts0[:len(MaxCounts0)/2]
    RamanIntCounts0 = IntCounts0[len(IntCounts0)/2:]
    RamanMaxCounts0 = MaxCounts0[len(MaxCounts0)/2:]
    days = [0,3,8]
    
    DFIntCounts8 = np.nan_to_num(DFIntCounts8)
    DFMaxCounts8 = np.nan_to_num(DFMaxCounts8)
    RamanIntCounts8 = np.nan_to_num(RamanIntCounts8)
    RamanMaxCounts8 = np.nan_to_num(RamanMaxCounts8)
    DFIntCounts3 = np.nan_to_num(DFIntCounts3)
    DFMaxCounts3 = np.nan_to_num(DFMaxCounts3)
    RamanIntCounts3 = np.nan_to_num(RamanIntCounts3)
    RamanMaxCounts3 = np.nan_to_num(RamanMaxCounts3)
    DFIntCounts0 = np.nan_to_num(DFIntCounts0)
    DFMaxCounts0 = np.nan_to_num(DFMaxCounts0)
    RamanIntCounts0 = np.nan_to_num(RamanIntCounts0)
    RamanMaxCounts0 = np.nan_to_num(RamanMaxCounts0)
                 
    print(DFIntCounts8)
    print(DFMaxCounts8)
    print(RamanIntCounts8)
    print(RamanMaxCounts8)
    print(DFIntCounts3)
    print(DFMaxCounts3)
    print(RamanIntCounts3)
    print(RamanMaxCounts3)
    print(DFIntCounts0)
    print(DFMaxCounts0)
    print(RamanIntCounts0)
    print(RamanMaxCounts0)
    
    fig1 = plt.figure()
    ax1 = fig1.add_subplot(111)
    for i in range(0,3):
        ax1.plot(days, [DFIntCounts0[i], DFIntCounts3[i], DFIntCounts8[i]])
    plt.title('DF Integrated Counts')
    plt.xlabel('Age of sample (days)')
    plt.ylabel('Integrated counts')
    plt.legend(['Rings','Spots','Asymmetrics'])
    plt.tight_layout()
    plt.show()
    
    fig2 = plt.figure()
    ax2 = fig2.add_subplot(111)
    for i in range(0,3):
        ax2.plot(days, [DFMaxCounts0[i], DFMaxCounts3[i], DFMaxCounts8[i]])
    plt.title('DF Mean of 10 Max Counts')
    plt.xlabel('Age of sample (days)')
    plt.ylabel('Integrated counts')
    plt.legend(['Rings','Spots','Asymmetrics'])
    plt.tight_layout()
    plt.show()
    
    fig3 = plt.figure()
    ax3 = fig3.add_subplot(111)
    for i in range(0,3):
        ax3.plot(days, [RamanIntCounts0[i], RamanIntCounts3[i], RamanIntCounts8[i]])
    plt.title('Raman Integrated Counts')
    plt.xlabel('Age of sample (days)')
    plt.ylabel('Integrated counts')
    plt.legend(['Rings','Spots','Asymmetrics'])
    plt.tight_layout()
    plt.show()
    
    fig4 = plt.figure()
    ax4 = fig4.add_subplot(111)
    for i in range(0,3):
        ax4.plot(days, [RamanMaxCounts0[i], RamanMaxCounts3[i], RamanMaxCounts8[i]])
    plt.title('Raman Mean of 10 Max Counts')
    plt.xlabel('Age of sample (days)')
    plt.ylabel('Integrated counts')
    plt.legend(['Rings','Spots','Asymmetrics'])
    plt.tight_layout()
    plt.show()
    