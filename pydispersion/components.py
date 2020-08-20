#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Thu Jun 18 12:48:47 2020

@author: lafforgue
version 1.0

Amelioration à faire: 
    -troisième ordre en propagation (okay) 
    -insérer la structure de la classe dans les commentaires

containing the medium class, to create un object associated to the medium
for the moment the data are stored in a dictionnary and not a csv file, see later for the csv file reader


"""
import numpy as np
from numpy import cos, sin, arcsin
try:
    from . import pulse
except ImportError: 
    import pulse
    
from scipy.misc import derivative
c=0.299792458

class Matter:
    media={
            'BK7':list([1,[1.03961212,0.00600069867],[0.231792344,0.0200179144],[1.01046945,103.560653]]),
            'SiO2':list([1,[0.6961663,0.0684043**2],[0.4079426,0.1162414**2],[0.8974794,9.896161**2]]),
            'BAF10':list([1,[1.5851495, 0.00926681282],[0.143559835,0.0424489805],[1.08521269, 105.613573]]),
            'SF5': list([1, [1.52481889, 0.011254756], [0.187085527, 0.0588995392], [1.42729015, 129.141675]]),
            'LASF9': list([1, [2.00029547, 0.0121426017],[0.298926886,0.0538736236], [1.80691843, 156.530829]])
            
    } #dict with all the medium, to be completed
    """The matter is implemented with sellmeyer equation as following :
    n^2= A + B_1l^2/(l^2-C_1)+...  -->  [A,[B_1,C_1],[B_2,C_2], ... ] """
    
    def __init__(self, name): 
        """to create a matter
        name = "BK7", "SiO2", ... """
        try:
            self.indice=self.media[name]
        except KeyError: 
            print("The medium choose is not registered on the list, to list the medium use Matter.media /n the medium has been replaced by SiO2 ")
            self.indice=self.media['SiO2']
            
        self.name=name
    
    def propagation(self,z,pulse_in, third_order=False):
        """compute the propagation for a material at a given distance distance z
        z : float,  propagation distance in micrometers
        pulse_in : pulse object, pulse at the entrance of the material
        h : precision for the different derivatives of the refractive indexe
        third_order: bool , True if we want to take into consideration the TOD in the dispersion
        """
        pulse_in.l0=pulse_in.parameters["l0"]
        pulse_in.w0=pulse_in.parameters["w0"]
        self.derivatives(pulse_in.l0) #calculation of k_0 k'_0 and k''_0
        self.pulse_out=pulse.Pulse(np.copy(pulse_in.X),np.copy(pulse_in.Y), np.copy(pulse_in.phase)) #creation of a new object pulsation
        
        self.pulse_out=pulse_in.Copy()
        
        self.pulse_out.w0=pulse_in.w0
        self.phase_prop=(self.k0+self.k_0*(self.pulse_out.w0)+self.GDD/2*(self.pulse_out.X-self.pulse_out.w0)**2)*z # computing the propagation by adding k(w)*z to the phase
        if third_order: 
            self.phase_prop+=self.TOD/6*(self.X_omega_self.w0)**3
        self.pulse_out.phase+=self.phase_prop
        
        self.pulse_out.Y=self.pulse_out.Y*np.exp(1j*self.phase_prop)
        
        self.pulse_out.parameters['delay']=self.k_0*z
        self.pulse_out.FreqtoTime(update=True)
        return self.pulse_out #return a new pulsation corresponding to the one after the propagation
    
    def calculate_indice(self, l):
        """l: float, central wavelenght in micrometers
        calculate the numerical value of the refractive index for a given wavelenght"""
        
        n=self.indice[0]
        for i in range(1,len(self.indice)): #computing sellmeyer equation
            n+=self.indice[i][0]*(l**2)/(l**2-self.indice[i][1])
        return np.sqrt(n)
    
    def derivatives(self,l, h=10**(-10)):
        """ l: float, central wavelenght in micrometer
        return an array with k(w0),k'(w0),k''(w0) """
        self.dn=(self.calculate_indice(l+h)-self.calculate_indice(l))/h
        self.d2n=derivative(self.calculate_indice,l,dx=0.001, n=2, order=3)
        self.d3n=derivative(self.calculate_indice,l,dx=0.001, n=3, order=7)

        
        self.k0=2*np.pi/l
        self.k_0=(self.calculate_indice(l)-l*self.dn)/c
        self.GDD=(l**3)/((c**2)*2*np.pi)*self.d2n
        facteur=-(l**2/(2*np.pi*c))**2/c
        self.TOD=facteur*(3*(self.d2n) + l*self.d3n)
        
        return np.array([self.k0, self.k_0, self.GDD, self.TOD])


#--------------------------------------------------------------------------------------------------------------------
# A tester

class GratingCompressor:
    """class to create a grating compressor with an incident light with an angle gamma and a grating space= d"""
    def __init__(self, d=1):
        """gamma = incident angle of the light
            d = grating density in 1/micrometer"""
        self.d=d
    def GDDcalculate(self, l0, gamma, z):
        """calculate the GDD coefficient
        l0 : centrale wavelenght in micrometers
        gamma : incident angle in radians
        z : distance between the 2 gratings in micrometer
        """
        self.gamma=gamma
        self.z= z
        self.GDD=-l0**3/(np.pi*(c**2)*(self.d**2))*(1-(l0/self.d-np.sin(gamma))**2)**(-3/2)*self.z # calculate GDD with the formula
        self.ng=1+(l0/self.d-np.sin(gamma))*np.sin(gamma) # compute the group refractive index
        return self.GDD

    def propagation(self, pulse_in, gamma, z):
        """compute the propagation throug a simple pair of grating compressor: 
        pulse_in = pulse object
        gamma= angle of incidence
        z : space between the two gratings """
        self.l0=pulse_in.l0
        pulse_out=pulse_in.Copy() #create a new pulse
        pulse_out.phase+=(pulse_in.w0-pulse_in.X)**2*self.GDD*z+ z*self.ng*self.w0 #add to the phase
        return pulse_out
#--------------------------------------------------------------------------------------------------------------------
class DoublePrismCompressor: 
    """creation of a double prism compressor"""
    def __init__(self, alpha=np.pi/3, matter_name="SiO2",  Brewster=False):
        """initialize the prism compressor with the matter and the apex angle
        alpha : float,  apex_angle in radian
        matter_name : string,  matter of the prism (see the list in Matter.medium)
        Brewster : bool, if we are in the Brewster configuration, then alpha will be calculated automaticaly  """

        self.glass_prism=Matter(matter_name)
        self.alpha=alpha  
        self.Brewster=Brewster
        

    def PathCalculate(self, l0, l1, w, h):
        """compute the optical path depending on lambda as a fonction of l (lambda)
        according to the scheme, path -> P= n(l)AB + BC +n(l)CD'
        cf schema in the report and in the library
        l0 : float, central wavelenght in micrometer
        l1 : float, distance beam to apex
        w : float,  x distance between prism
        z : float, y distance between prism
        """
        self.l1, self.w, self.h = l1, w, h #add the varibla to the object
        
        self.nl=self.glass_prism.calculate_indice(l0)
        if self.Brewster: #if we are in brewter configuration
            self.gamma1=np.arctan(self.nl)
            self.alpha=2*arcsin(sin(self.gamma1)/self.nl)
            
        self.delta1=arcsin(sin(self.gamma1)/self.nl) #intern refraction
        self.delta2=self.alpha-self.delta1 #intern reflexion
        self.gamma2=arcsin(self.nl*sin(self.delta2)) # refractive angle 2nd surface)
        self.l1_=self.l1*cos(self.delta1)/cos(self.delta2) #distance to apex for the first prism output
        
        self.l2=(self.h*cos(self.gamma2-self.alpha/2)-self.w*sin(self.gamma2-self.alpha/2))/cos(self.gamma2)-self.l1_ #distance to apex of the second prism entrance


        self.ab=sin(self.alpha)*self.l1/cos(self.delta2) #path in the first prism
        self.bc=(self.w*cos(self.alpha/2)-self.h*sin(self.alpha/2))/cos(self.gamma2) #path lenght in the air between two prism
        self.cd=sin(self.alpha)*self.l2/cos(self.delta2) #path in the second prism
        self.lg=self.ab+self.cd #Path lenght in a glass
        

    def GDDcalculate(self, l0, l1, w, h, gamma1=0): 
        """calculate GDD with the path Calculate using GDD = l0^2/(2*pi*c)*d²p/dl² 
        return the GDD coefficient"""
        self.gamma1=gamma1
        self.glass_prism.derivatives(l0) #implementation of the prism matter
        self.pathCalculate(l0, l1, w, h) #calculation of the different paths
        self.GDD=2*l0**3/(2*np.pi*c**2)*(self.lg*self.glass_prism.d2n-(4*self.bc+self.lg/(self.nl**3))*self.glass_prism.dn**2)
        self.GDD+=self.glass_prism.GDD*self.lg #add the GDD from the glass propagation
        return self.GDD

    
    def propagation(self, pulse_in, l1, w, h, gamma1=0): 
        """compute the propagation with the different parameters calculated before
        the propagation through a dispercive medium like the glass is taken into consideration of the pulse is taken into consideration
        pulse_in : pulse type object 
        l1 : float, distance entrance-apex in micro meter
        w : float,  x-abcisse apex-apex distance in micrometer
        h : float, y-abcisse apex-apex distance in micrometer
        gamma1 : float, incident angle in radian
        """
        self.GDDcalculate(pulse_in.parameters['l0'], l1, w, h, gamma1)
        self.pulse_out=pulse_in.Copy()
        
#add to the phase just for the polynomial version
        self.pulse_out.phase+=self.GDD/2*(pulse_in.X-pulse_in.parameters["w0"])**2
        self.pulse_out.Y*=np.exp(1j*self.GDD/2*(pulse_in.X-pulse_in.parameters["w0"])**2)
        return self.pulse_out

##################################################code_test class#####################################################
