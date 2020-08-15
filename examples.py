"""file containning few example to understand how the module works"""
"""changer les extensions quand ce sera publier """
#from . import components
#from . import graphic
#from . import pulse
from pydisp import components, graphic, pulse
##example 1 : propagation of experimental impulsion after  mm of SiO2
def experimental_propagation():
    pulse1=pulse.Data_input("test.dat", axex="lambda")
    pulse1.SetParameters() #calculation of frequency parameters
    pulse1.SetTimeParameters() #calculation of temporal parameters
    
    glass=components.Matter('SiO2') # creation of the material
    #
    pulse_out=glass.propagation(10000, pulse1) # propagation after 1 mm, z is given in micrometers
    pulse_out.SetParameters() # calculation of the paramaters
    pulse_out.SetTimeParameters()
    #
    graphique=graphic.BaseGraph()
    graphique.compare(pulse1, pulse_out)

##example 2 creation of an gaussian pulse, propagation and plotting
def gaussian_propagation():
    pulse_in = pulse.Analytics(shape="Gaussian", w0=2.5, E0_t=1, HMD=10) #creation of a gaussian pulse with amplitude of 1 and half maximum duration of 10 fs
    pulse_in.SetParameters() #calculation of frequency parameters
    
    glass=components.Matter('SiO2') # creation of the material
    
    pulse_out=glass.propagation(1000, pulse_in) # propagation after 500 micrometers
    pulse_out.SetParameters() # calculation of the paramaters
    pulse_out.SetTimeParameters()
    
    graphique=graphic.BaseGraph()
    graphique.compare(pulse_in, pulse_out)

##Example 3: manipulation of objects: 
def manipulation():
    glass= components.Matter("BK7")
    print(glass.derivatives(0.8)) # print the different parameters k0 k'0 and k''0=GDD for a wavelenght of 800 nm
    
    
    pulse1=pulse.Data_input("test.dat", axex="lambda")
    pulse1.SetParameters() #calculation of frequency parameters
    pulse1.SetTimeParameters() #calculation of temporal parameters
    
    pulse1 # print a dic with all the parameters of the pulse
    
    pulse1.X # frequency abcisse  array
    pulse1.Y # Field spectrum
    pulse1.X_time #time abcisse
    pulse1.Y_time #temporal complexe field 
    
##example 4: dispersion compensation throught a prism compressor
def double_prism():
    pulse_in = pulse.Analytics(shape="Gaussian", w0=2.5, E0_t=1, HMD=10) #creation of a gaussian pulse with amplitude of 1 and half maximum duration of 10 fs
    pulse_in.SetParameters() #calculation of frequency parameters
        
    compressor= components.DoublePrismCompressor(Brewster=True)
    pulse2=compressor.propagation(pulse_in, 10000, 300,100000)
    print(compressor.GDD)
    pulse2.SetParameters() #calculation of frequency parameters
    pulse2.SetTimeParameters() #calculation of temporal parameters
        
    graphique=graphic.BaseGraph()
    graphique.compare(pulse_in, pulse2)
    

########## figure for the report dans the result #####
#pulse1=pulse.Data_input("test.dat", axex="lambda")
#pulse1.SetParameters() #calculation of frequency parameters
#pulse1.SetTimeParameters() #calculation of temporal parameters
##
#pulse2=pulse.Data_input("test.dat", axex="omega", smoothing=False)
#
#graphique=graphic.BaseGraph()
#graphique.plotattribut(pulse2, legend_X="wavelenght in nm",  label="Spectrum")

#graphique.plotattribut(pulse1, X="X_time", Y="Y_time", legend_X="Time in fs", legend_Y="Field Amplitude in a.u", label="real field")
