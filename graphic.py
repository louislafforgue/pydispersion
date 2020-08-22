#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Fri Jun 19 13:51:14 2020

@author: lafforgue

module graphique
Basicly take a pulse object in argument and plot what we want thanks to different methods
"""
import matplotlib.pyplot as plt
import numpy as np
#import numpy as np

class BaseGraph:

    def __init__(self, title="Pulse representation"):
        self.title = title
        self.x_label = "X-axis label"
        self.y_label = "Y-axis label"
        self.show_grid = True

    def envelop_plot(self, pulse, repre="freq"):
        # x_values = gather only x_values from our zones
        # y_values = gather only y_values from our zones
        if repre=='freq': 
            plt.plot(pulse.X, pulse.Y)
        else: 
            plt.plot(pulse.X_time, np.abs(pulse.Y_time))
        plt.xlabel(self.x_label)
        plt.ylabel(self.y_label)
        plt.title(self.title)
        plt.grid(self.show_grid)
        plt.show()
    
    def comparedouble(self, pulse1, pulse2): 
        plt.title("Comparaison between two pulses")
        self.fig1= plt.subplot(211) #left top 
        self.fig1.plot(pulse1.X,np.abs(pulse1.Y), label="pulse1")
        self.fig1.plot(pulse2.X,np.abs(pulse2.Y), label="pulse2")
        self.fig1b=self.fig1.twinx()
        self.fig1b.plot(pulse2.X,(pulse2.phase-pulse1.phase), 'r--',label="phase difference" )
        
        self.fig1.set_title("spectral representation ")
        self.fig1.set_ylabel("Field amplitude a.u")
        self.fig1.set_xlabel("angular frequency w")
        self.fig1.legend(["pulse1","pulse2"],loc='upper left')
        self.fig1b.legend(["phase difference"])
        
        self.fig2=plt.subplot(212) #right top 
        self.fig2.plot(pulse1.X_time,np.abs(pulse1.Y_time), label="first")
        self.fig2.plot(pulse2.X_time,np.abs(pulse2.Y_time), label="second")
        self.fig2.legend(["pulse1","pulse2"])
                
            
        self.fig2.set_title("temporal representation with {} fs of difference".format(np.around(pulse2.parameters['delay']-pulse1.parameters['delay'])))
        self.fig2.set_xlim(-3*max(pulse1.parameters['HMD'], pulse2.parameters['HMD']), 3*max(pulse1.parameters['HMD'], pulse2.parameters['HMD']) )
        self.fig2.set_ylabel("Time Field amplitude in a.u")
        self.fig2.set_xlabel("time in femtosecond")
        plt.tight_layout(pad=0.1)
    
    def compare(self, pulse1, pulse2):
        """plotting of the two pulses in all the representation
        split in 4 parts
        -If the temporal representation was not calculated before, the algorithm is calculating it
        """
        plt.title("pulse propagation after a given distance")
        self.fig1= plt.subplot(221) #left top 
        self.fig1.plot(pulse1.X,np.abs(pulse1.Y))
        self.fig1b=self.fig1.twinx()
        self.fig1b.plot(pulse1.X,pulse1.phase, 'r--')
        self.fig1.set_title("pulse1 spectral representation ")
        self.fig1b.set_ylabel("phase")
        self.fig1.set_ylabel("Field amplitude a.u")
        self.fig1.set_xlabel("angular frequency w")
        
        self.fig2=plt.subplot(222) #right top 
        try:
            self.fig2.plot(pulse1.X_time,pulse1.Y_time)
        except AttributeError: 
            pulse1.FreqtoTime()
            self.fig2.plot(pulse1.X_time,pulse1.Y_time)
        self.fig2.set_title("temporal representation {} fs of delay".format(np.around(pulse1.parameters['delay'])))
        self.fig2.set_xlim(-3*pulse1.parameters['HMD'], 3*pulse1.parameters['HMD'] )
        self.fig2.set_ylabel("Time Field amplitude in a.u")
        self.fig2.set_xlabel("time in femtosecond")
        
        
        self.fig3=plt.subplot(223) #left bottom
        self.fig3.plot(pulse2.X,np.abs(pulse2.Y))
        self.fig3b=self.fig3.twinx()
        self.fig3b.plot(pulse2.X,pulse2.phase,'r--')
        self.fig3.set_title("spetral representation pulse2")
        self.fig3b.set_ylabel("phase")
        self.fig3.set_ylabel("Field amplitude a.u")
        self.fig3.set_xlabel("angular frequency w")
        
        self.fig4=plt.subplot(224) #right bottom
        try:
            self.fig4.plot(pulse2.X_time,pulse2.Y_time)
        except AttributeError: 
            pulse2.FreqtoTime()
            self.fig4.plot(pulse2.X_time,pulse2.Y_time)
        self.fig4.set_title("temporal representation pulse2 with {} fs of delay".format(np.around(pulse2.parameters['delay'])))
        self.fig4.set_xlim(-3*pulse2.parameters['HMD'], 3*pulse2.parameters['HMD'] )
        self.fig4.set_ylabel("Time Field amplitude in a.u")
        self.fig4.set_xlabel("time in femtosecond")
        
        plt.tight_layout(pad=0.1)
    def CompareEnvelop(self, pulse1, pulse2):
        """plotting of the envelops of two pulses either in temporal (repre="time") or frequencial (repre="freq") domain"""
        plt.plot(pulse1.X,pulse1.Y)
        plt.plot(pulse2.X,pulse2.Y)
        plt.set_ylabel("Field amplitude a.u")
        plt.set_xlabel("angular frequency w")
        
    def CompareTimeEnvelop(self, pulse1, pulse2):
        """plot the time envelop of two pulses: 
            pulse1 : pulse object
            pulse2 : pulse object 
            """
        plt.plot(pulse1.X_time,np.abs(pulse1.Y_time), label="first")
        plt.plot(pulse2.X_time,np.abs(pulse2.Y_time), label="second")
        plt.legend(["first","second"])
        plt.ylabel("Time Field amplitude a.u")
        plt.xlabel("Time in femtosecond")
    
    def plotTime(self, pulse): 
        """plot the time field and the time envelop of a pulse
        pulse: pulse object"""
        plt.plot(pulse.X_time,pulse.Y_time)
        plt.plot(pulse.X_time,np.abs(pulse.Y_time), 'r--')
        plt.ylabel("Field amplitude a.u")
        plt.xlabel("Time in femtosecond")
    
    
    def plotattribut(self, pulse, X="X", Y="Y", legend_X="frequency 10^15 rad/s", legend_Y="Field amplitude a.u",label="", absolute=False, **kwargs): 
        """plot the asked graph : 
            pulse : pulse object 
            X : string name of the attribut in abcisse
            Y : string name of the attribut in ordonee
            lengend_X: string abcisse legend
            legend_Y: string ordonnee legend
            absolute: bool, take the the absolut value of the Y signal
            """
        X_plot=getattr(pulse, X)
        Y_plot=getattr(pulse, Y)
        if absolute: 
            Y_plot=np.abs(Y_plot)
        plt.plot(X_plot, Y_plot,label="line", **kwargs)
        plt.legend([label])
        plt.xlabel(legend_X)
        plt.ylabel(legend_Y)
            
        
        
        
        
        
#class TimeGraph(BaseGraph):
#
#    def __init__(self):
#        super().__init__()
#        self.title = "Time representation"
#        self.x_label = "time"
#        self.y_label = "intensity"
#        
#class FreqGraph(BaseGraph):
#    
#    def __init__(self):
#        super().__init__()
#        self.title = "Frequency representation"
#        self.x_label = "frequency"
#        self.y_label = "intensity"
###########################################################code test class###############################################
    