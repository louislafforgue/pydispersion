#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:55:37 2020

@author: lafforgue

file to test the different class and methods

test: 
    pulse file
    
        SetParameters
        FreqtoTime
        Data_input (revoir à cause du sqrt)
    
    Components file: 
        matter: refractive_index, GDD-TOD calculation, gaussian_propagation
        grating: GDD, (TOD) calculation
        DoubleprismCompressor: brewter configuaration, GDD TOD calculation

"""
import numpy as np
import pulse
import matplotlib.pyplot as plt
import components
import graphic
from scipy.fftpack import fft, fftfreq, fftshift, ifft 

c=0.299

def test_SetParameters():
    """should return [0,0,0] or close, it gives the error between theoretical and calculation"""
    pulse1=pulse.GaussianSpectrum(2,1,0.7)
    pulse2=pulse.Pulse(np.copy(pulse1.X), np.copy(pulse1.Y))
    pulse2.SetParameters()
    
    return [pulse1.parameters["w0"]-pulse2.parameters["w0"], pulse1.parameters["E0_f"]-pulse2.parameters["E0_f"], pulse1.parameters["FWHM"]-pulse2.parameters["FWHM"]]

def test_FreqtoTime(): 
    """test the difference between the analytical and the numerical fourier transform
    return the relative error"""
    pulse1=pulse.Analytics(shape="Gaussian", w0=2, E0_t=1, HDM=10) # 10 femtosecond pulse
    pulse2=pulse.Pulse(np.copy(pulse1.X), np.copy(pulse1.Y))
    pulse2.SetTimeParameters()
    return [ np.abs((pulse1.parameters["E0_t"]-pulse2.parameters["E0_t"])/pulse1.parameters["E0_t"]),np.abs((pulse1.parameters["HMD"]-pulse2.parameters["HMD"])/pulse2.parameters["HMD"])]

def test_TimetoFreq():
    pulse1=pulse.Analytics(shape="Gaussian", w0=2, E0_t=1, HDM=10) 
    pulse1.SetParameters()
    Y=pulse1.TimetoFreq(update=False)
    pulse2=pulse.Pulse(np.copy(pulse1.X),Y)
    pulse2.SetParameters()
    return [ np.abs((pulse1.parameters["E0_f"]-pulse2.parameters["E0_f"])/pulse1.parameters["E0_f"]),np.abs((pulse1.parameters["FWHM"]-pulse2.parameters["FWHM"])/pulse2.parameters["FWHM"])]
    
    

def Data_input_test(): 
    """return error in nanometer for l0 and wavelenght bandwith 
    réécrire le test par rapport à un pulseation dont on connait les propriétés 
    ex: gaussian pulse"""

    pulse1=pulse.Data_input("test.dat", axex="lambda")
    pulse1.SetParameters()
    pulse1.SetTimeParameters()
    
    pulse_test=pulse.Pulse(np.copy(np.abs(np.abs(pulse1.X_data))), np.copy(np.abs(np.sqrt(pulse1.Y_data))))
    pulse_test.SetParameters()
    print(pulse1)
    print(pulse_test)
    return np.abs(pulse1.parameters["l0"]*1000-pulse_test.parameters["w0"]), pulse1.parameters["l0"]*pulse1.parameters["FWHM"]/(2*np.pi*c)-pulse_test.parameters["FWHM"]/1000/pulse1.parameters["w0"]

def test_calculate_indice():
	"""return the error in the indice calculation with sell meyer equation"""
	glass=components.Matter("SiO2")
	return glass.calculate_indice(0.8)-1.453 

def test_derivatives():
	glass=components("SiO2")
	return glass.derivatives() # to be completed with the data found on the web

def test_propagation(x_f): 
    """plot a the diiference between the théorical expression and the calculate value for a propagation throught BK7 glass after a distance x_f"""

    pulse_in = pulse.Analytics(shape="Gaussian")


    glass=components.Matter('BK7') # creation of the material
    pulse_out=glass.propagation(x_f, pulse_in)
    pulse_out.SetParameters()
    pulse_out.SetTimeParameters()

    gamma=4*np.log(2)/pulse_in.HMD**2
    ksi=2*gamma*glass.GDD
    gamma_x=gamma/(1+ksi**2*x_f**2)-1j*gamma*ksi*x_f/(1+ksi**2*x_f**2)
    Y_theo=np.sqrt(np.pi/gamma)*np.sqrt(np.abs(gamma_x/np.pi))*np.exp(-gamma/(1+ksi**2*x_f**2)*pulse_out.X_time**2)

    plt.plot(pulse_out.X_time, np.abs(pulse_out.Y_time))
    plt.plot(pulse_out.X_time, Y_theo)
    plt.show()
    
    return pulse.Pulse.GetParameters(pulse_out.X_time, np.abs(pulse_out.Y_time))-pulse.Pulse.GetParameters(pulse_out.X_time, Y_theo)

def test_GratingCompressor():
    """return True if th result is right with +-5 fs2/um"""
    compr=components.GratingCompressor()
    result=compr.GDDcalculate(1.064, np.pi/6, 1000)
    return result, result > -7581 and result < -7571

def test_PrismCompressor(error=0.03):
    """return the relativ error for the different distances, 
    return TRUE for the calcul of the GDD if the relative error is inferior to error
    """
    compr=components.DoublePrismCompressor(matter_name="BK7",Brewster=True)
    compr.GDDcalculate(0.8, 5000 ,465005, 203358 )
    ab_theo=5520
    bc_theo=499150
    cd_theo= 5519
    GDD_theo= -441.015+(ab_theo+cd_theo)*compr.glass_prism.GDD
    return np.abs((compr.GDD-GDD_theo)/GDD_theo) < error, [np.abs((compr.ab-ab_theo)/compr.ab), np.abs((compr.bc-bc_theo)/compr.bc), np.abs((compr.cd-cd_theo)/compr.cd)]



#####################################################################################test random#######################################################


