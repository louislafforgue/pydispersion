#!/usr/bin/env python3
# -*- coding: utf-8 -*-
"""
Created on Wed Jun 17 10:24:28 2020

@author: lafforgue

version 1.1

"""

"""

definition of the different class to compute a pulse object 


pulse --- data_input
       |
       -- GaussianSpectrum
       -- analytics
       
pulse object store a pulsation as a 3*N array [X, Y, Phase], frequencial
possibility to add a phase as an array
possibility to add parameters in argument : ex shape= "gaussian"

methods in pulse: 
    # fitting -> report the approximate gaussian fonction with the parameters
    GetParameters -> staticmethod  gives [w0,FWMH, I0] for frequencial or [t0,HMD, I0] for a temporal
    SetParameters -> Calculate and set the frequency parameters in the dictionnary
    SetTimeParameters -> idem with time parameters
    FreqToTime -> compute the inverse fourier transform
    TimeToFreq -> compute the fourier transform
    ZeroAdding 
    save
    copy

    
methods in Data_input: 
    reading -> take a filename and read the the first column X and the second column Y
    replace -> static method for format processing
    convert -> static methods for csv files
    interpolation -> reproduce the spectrum with linear X spacing
    processing -> smooting and remoove the negative values
    

Analytics:  
    Pulses which can be created are in time domain:
        -Gaussian "Gaussian": I0exp(-(t/tau)²/2) 
        -hyperbolique secante "Hsech" : I0/cosh(t/tau)  
        -Lorentizian "Lorentzian" : I0(1+(t/tau)^2)^-1 
    as the parameters as follow [w0, I0 , HMD] with HMD the duration of the pulse (time fWMH)
    In default: the pulse will be a gaussian
    
Example to declare an object: 
    pulse1=pulse.Data_input("test.dat", axex="lambda")
    pulse_in = pulse.Analytics(shape="Gaussian", param=[2.5,1,5])
    
"""
import numpy as np
import csv
from scipy.fftpack import  fftfreq, fftshift, ifft, fft
from scipy import interpolate
from scipy.signal import savgol_filter

c=0.299

#---------------------------------------------------------------------------------------------------------------------
class Pulse: 
    
    def __init__(self, pulse_X,pulse_Y,pulse_phase=[], X_time=[], N_grid=2**15, **kwargs):
        """object construction with 3 arrays plus a dictionnary parameters
        pulse_X = ordered real array, angular frequency in 10^15
        pulse_Y= real positive array
        **kwargs allow to enter manualy some parameters which would be store the dictionnay parameters, ex: shape="gaussian" """
        
        self.N_grid=N_grid #Size for fourier transform
        self.N_signal=np.size(pulse_X)
        
        self.f_max=pulse_X[-1]
        
        self.X=pulse_X
        self.Y=pulse_Y
        
        self.normf=2**np.log2(self.N_grid/self.N_signal)*(self.f_max/(2*np.pi)) #normalisation factor for the fourier transform due to 0-padding
        
        self.df=(self.X[-1]-self.X[0])/(2*np.pi*(np.size(self.X)-1))
        
        
        if len(X_time)==0:
            self.X_time=fftshift(fftfreq(self.N_grid,self.df)) #time from the fourier grid
            self.X_time=self.X_time[N_grid//2-self.N_signal//2:self.N_grid//2+self.N_signal//2] #size N_signal
        
        if len(pulse_phase)!= self.N_signal:
            self.phase=np.zeros(self.N_signal,dtype='complex128') #creation an empty array for the phase
        else: 
            self.phase=pulse_phase
         
        #we can add some parameters to the object pulse
        self.parameters=dict()
        for key, value in kwargs.items():
            self.parameters['{}'.format(key)]=value 
   
    def __str__(self):
        return "This pulsation has the parameters : {}".format(self.parameters)
    
    def __repr__(self):
        """the representation return the dictionnary of the parameters which cam be directly used"""
        return str(self.parameters)
    
    @staticmethod
    def GetParameters(X,Y,precision=1.e-2):
        """calculate the different the parameters
        for frequencial pulse: [w0, I0, FWMH]
        take an intensity profile and return a tab with the 3 parameters
        
        """ 
        I0=np.amax(Y) #store the maximum value
        
        in_fwhm=False #used to 
        i_fwhm=[0,0] #variable containing the indices of the fwhm
        
        i_w0=[] #store the list of the different w0, the list will be order
        for i in range(np.size(Y)):
            
            #calculation of w0
            if Y[i] > I0-precision:
                i_w0.append(i)
            
            #calculation of FWMH   
            if not in_fwhm and Y[i] >= I0/2:
                in_fwhm=True
                i_fwhm[0]=i
            if in_fwhm and Y[i] <= I0/2:
                in_fwhm=False
                i_fwhm[1]=i     
                
        #numerical calculation of w0
        w0=X[int((i_w0[0]+i_w0[-1])/2)]
        
        #numerical calculation of FWMH 
        fwmh=X[i_fwhm[1]]-X[i_fwhm[0]]
        
        return np.array([w0,I0,fwmh])
    
    def SetParameters(self, precision=1.e-2):
        """calculate the parameters and store them in the dictionnary parameters
        store them as {'w0', 'l0', 'E0_f', 'FWMH' }"""
        self.param=self.GetParameters(self.X,self.Y, precision)
        self.parameters['w0']=self.param[0]
        self.parameters['l0']=2*np.pi*c/self.param[0]
        self.parameters['E0_f']=self.param[1]
        self.parameters['FWHM']=self.param[2]
        
        
    def SetTimeParameters(self, precision=1.e-2):
        """calculate the temporal parameters thanks to a fourier transform
        then store the them in the dicionary parameters
        {'delay','E0_t','HMD'}
        """
        self.FreqtoTime()
        self.param=self.GetParameters(self.X_time, np.abs(self.Y_time))
        try:
            self.parameters['delay']+=self.param[0]
        except KeyError: 
            self.parameters['delay']=self.param[0]
        self.parameters['E0_t']=self.param[1]
        self.parameters['HMD']=self.param[2]
        
    def FreqtoTime(self, update=True):
        #centered padding
#        Y_pad=np.hstack([np.zeros(self.N_grid//2-self.N_signal//2),self.Y,np.zeros(self.N_grid//2-self.N_signal//2)])
        #compute the inverse fourier tranform
        Y_signal=self.Y*np.exp(-1j*self.phase)
        Y_time=fftshift(ifft(Y_signal, self.N_grid))*self.normf #size N_grid
        Y_signal=Y_time[self.N_grid//2-self.N_signal//2:self.N_grid//2+self.N_signal//2]  #size N_signal
        
        if update==True: 
            self.Y_time=np.abs(Y_signal)
            self.phase_time=np.angle(Y_signal)
        return Y_signal  #size N_signal

    def TimetoFreq(self, update=True): 
        Y_signal=self.Y_time*np.exp(1j*self.phase_time)
        #centered padding
        Y_time_pad=np.hstack([np.zeros(self.N_grid//2-self.N_signal//2),Y_signal,np.zeros(self.N_grid//2-self.N_signal//2)])
        #fourier transform calculation
        Yf_omega=fftshift(fft(Y_time_pad))/self.normf
        Yf_omega=Yf_omega[self.N_grid//2:self.N_grid//2+self.N_signal]
        if update==True: 
            self.Y=np.abs(Yf_omega)
            self.phase_time=np.angle(Yf_omega)
        
        return Yf_omega #size N_signal


    def Save(self, file_name):
        """save the pulsation in a .dat file named file_name
        à revoir avec les paramètre à entre en première ligne"""
        np.savetxt(file_name, (self.X,self.Y,self.phase), delimiter=',') 
        pass
    def Copy(self):
        """return a copy of the pulse as a new object
        à tester
        """
        return Pulse( np.copy(self.X), np.copy(self.Y), np.copy(self.phase), N_grid=2**15,  **self.parameters)
        
#-----------------------------------------------------------------------------------------
class GaussianSpectrum(Pulse): 
    """creation of gaussian spectrum, add the comments or add in the class analytics"""
    
    def __init__(self, w0=2, E0_f=1, FWHM=0.5, pulse_phase=[], N_grid=2**15, N_omega=2**11, **kwargs):
        """add the comments"""
        
        self.N_grid=N_grid
        self.N_omega=N_omega
        self.w0=w0
        self.max_frequency=self.w0+3*FWHM
        self.FWHM=FWHM
        self.E0_f=E0_f
        self.X=np.linspace(0,self.max_frequency, self.N_omega)
        self.Y=self.E0_f*np.exp(-np.log(2)*4*((self.X-self.w0)/self.FWHM)**2)
        super().__init__(self.X,self.Y,pulse_phase,self.X_time, self.Y_time, shape="gaussian", w0=self.w0, E0_f=self.E0_f, FWHM=self.FWHM, **kwargs) 


class Analytics(Pulse): 
    """creation of an analytical object
    The optional arguments are N (number of step or the time domaine)
    shape : "Gaussian" , "Hsech", "Lorentzian" 
    """
    
    
    def __init__(self, shape="Gaussian", f_width=6, w0=2, E0_t=1, HMD=10, poly_phase=[0] ,repre="time",N_grid=2**15, N_signal=2**11, **kwargs):
        """w0: angular frequency in 10^15 rad/S
        E0_t : float intensity in a.u
        HMD: float Half maximum duration in fs
        pulse_phase= list(float) , list with the polynomial coffecients, ex: [1,2.3,5] = 1+2.3T+5²T
        N_grid: int, more efficient if N=2^k discretization points for the temporal signal
        N_omega: int, number of points for the spectrum signal 
        repre: "time" or "freq", default "time"
        """
        
        self.N_grid=N_grid
        self.N_signal=N_signal
        
        self.w0=w0
        if repre=="freq": 
            self.I0_f=E0_t
            self.FWHM=HMD
        else:
            self.I0_t=E0_t
            self.HMD=HMD
        
#        self.X_time=np.linspace(-3*self.tau,3*self.tau,self.N_signal)
#        self.dt=self.X_time[1]-self.X_time[0]
#        self.X=2*np.pi*fftshift(fftfreq(self.N_grid,self.dt))[self.N_grid//2:self.N_grid//2+N_signal]
        
        if shape=="Hsech":
            self.K=np.pi*4*0.142 #time band product
            self.X_creation()
            self.Y_time=self.Hsech()
            
        elif shape=="Lorentzian": 
            self.K=np.pi*4*0.315
            self.X_creation()
            self.Y_time=self.Lorentzian()
            
        else: 
            self.K=np.pi*4*0.441
            self.X_creation()
            self.Y_time=self.Gaussian()
            
        #we add the parameters from the analytical expression
        self.temporal_phase=np.zeros(N_signal)
        for i in range(len(poly_phase)):
            self.temporal_phase+=poly_phase[i]*(self.X_time**i)
        
        self.Y_time=self.Y_time*np.exp(1j*self.temporal_phase)
        
        super().__init__(self.X,self.Y, shape=shape, w0=self.w0, E0_t=self.I0_t, HMD=self.HMD, delay=0, **kwargs)

        
    def Gaussian(self):
        """create a Gaussian spectrum and pulse shape in the time domain"""
        self.gamma=4*np.log(2)/(self.HMD**2)
        self.Y=self.I0_t*np.sqrt(np.pi/self.gamma)*np.exp(-(self.X-self.w0)**2/(4*self.gamma)) #theoretical fourier transform
        return self.I0_t*np.exp(-np.log(2)*4*(self.X_time/self.HMD)**2)*np.exp(1j*self.X_time*self.w0) 
    
    def Hsech(self):
        """create a HSECH spectrum and pulse shape in the time domain"""
        self.t0=self.HMD/2.634
        self.Y=self.I0_t*self.t0*np.pi/np.cosh(np.pi/2*self.t0*(self.X-self.w0))
        return self.I0_t*(1/np.cosh(self.X_time/self.t0))*np.exp(1j*self.X_time*self.w0)
    
    def Lorentzian(self):
        """create a Lorentzian spectrum and pulse shape in the time domain"""
        self.Y=self.I0_t*self.HMD*np.exp(-self.HMD*np.abs(self.X-self.w0))
        return self.I0_t/(1+(self.X_time/self.HMD)**2)*np.exp(1j*self.X_time*self.w0)
    
    def X_creation(self): 
        self.X=np.linspace(0,self.w0+3*self.K/self.HMD, self.N_signal)
        self.df=(self.w0+3*self.K/self.HMD)/(2*np.pi*(self.N_signal-1))
        self.X_time=fftshift(fftfreq(self.N_grid, self.df))[self.N_grid//2-self.N_signal//2:self.N_grid//2+self.N_signal//2]
        
        
#---------------------------------------------------------------------------------------------------------------   
class Data_input(Pulse):
    """class to create an object pulse from a data_input"""
    
    def __init__(self,file_name,phase=[], axex="lambda", smoothing=True, **kwargs): 
        """creation of an pulse object without data processing
        A perfect data should be given
        the file as to be given as following 1 column: frequency, 2 column Intensity amplitude
        
        the parameters for axex are : 'lambda', 'freq', 'omega'
        smoothing ( True or False): realize a smoothing
        
        peut etre modification du nom de variable x_lambda"""
        self.file_name=file_name
        self.reading() #reading the file
        if axex=="lambda":
            self.X_lambda=self.X_data
            self.Y_lambda=self.Y_data
            self.convert() #convert into angular frequency
            self.interpolation() # interpolate with the new list of angular frequency to have an linear step self.X
        elif axex=="freq": 
            self.X=self.X_data*2*np.pi
            self.Y=self.Y_data*2*np.pi
        elif axex=="omega": 
            self.X=self.X_data
            self.Y=self.Y_data
        else: 
            self.X=self.X_data
            self.Y=self.Y_data
            print(" Be carreful, maybe the X axis is not appropriate and the program will give you wrong result")    
        if smoothing:
            self.processing()
        self.Y=np.sqrt(np.abs(self.Y)) #we take the amplitude of the field
        self.phase=phase
        self.ZeroAdding()
        super().__init__(self.X,self.Y,phase, **kwargs)


        
    def reading(self):
        """reading the file (.dat, .csv')associate to file_name and store the data in Pulse object
        for a csv file : the first line is take off
        for a .dat file: all the file is taken
        """
        y_csv=[]
        x_csv=[]
        
        #reading of csv files
        if self.file_name[-4:]==".csv":

            with open(self.file_name) as csv_file:
                csv_reader = csv.reader(csv_file, delimiter=',')
                line_count = 0
                for row in csv_reader:
                    if line_count==0:
                        line_count +=1
                    else: 
                        x_csv.append(self.replace(row[0])) #adding x 
                        if float(self.replace(row[1])) > 0: #removing negative values from the experimental spectrum
                            y_csv.append(self.replace(row[1])) #adding y
                        else: 
                            y_csv.append(0)
                        line_count += 1
        
        #reading of other files
        else: 
            self.data=np.loadtxt(self.file_name)
            for i in range (np.shape(self.data)[0]):
                x_csv.append(self.data[i][0])
                if float(self.data[i][1]) > 0: #removing negative values from the experimental spectrum
                    y_csv.append(self.data[i][1]) #adding y
                else: 
                    y_csv.append(0)
                       
        #store the spectrum in two arrays 
        self.X_data=np.array(x_csv,dtype=float)
        self.Y_data=np.array(y_csv,dtype=float)

        return [np.array(x_csv,dtype=float),np.array(y_csv,dtype=float)]
    
    def convert(self):
        """convert from wavelenght spectrum to frequency spectrum using jacobian scalling"""
        self.X_lambda=self.X_lambda*1e-3 #conversion of wavelength in micrometers
        self.X_nl=2*np.pi*c/self.X_lambda #X conversion
        self.Y_nl=self.X_lambda**2/(2*np.pi*c)*self.Y_lambda #jacobian scaling cf rick trebino
        
    def interpolation(self): 
        """interpolate between the non linear spaced frequency and the new tab linear spaced to have usefull arrays"""
        self.N_signal=np.size(self.Y_nl) #size of the sample
        self.X=np.linspace(self.X_nl[-1],self.X_nl[0], self.N_signal) #creation of the new frequencies array, linear spaced
        self.df=(self.X[-1]-self.X[0])/self.N_signal
        self.Y=interpolate.interp1d(self.X_nl, self.Y_nl, kind="cubic")(self.X)
        
        
    @staticmethod
    def replace(liste):
        """method used to convert the str "3,5" to a str "3.5" which can be converted by float()
        liste =  srt or list """
        a=liste.find(",")
        if a >0:
            return liste[:a]+'.'+liste[a+1:]
        else: 
            return liste     

    def processing(self, N=51):
        """apply a smoothing to the data
        remove the noise and apply the savgol_filter
        N: parameter for savgol filter
        """
        for i in range(np.size(self.Y)): # remooving the noise at the extremities
            if self.Y[i] <4: 
                self.Y[i]=0
        
        self.Y=savgol_filter(self.Y, N, 3)
        for i in range(np.size(self.Y)): # remooving negatives values for the square root
            if self.Y[i] <0: 
                self.Y[i]=0
    
    def ZeroAdding(self):
        """add 0 to the amplitude for value to extend X abscisse, needed to have a good fourier transform"""
        tab_add= np.arange(0,self.X[0], (self.X[-1]-self.X[0])/np.size(self.X))
        self.X=np.hstack([tab_add, self.X])
        y_add=np.zeros(np.size(tab_add))
        self.Y=np.hstack([y_add, self.Y])
        self.phase=np.hstack([y_add, self.phase])
        self.N_signal=self.Y.size
                        
###########################################################code test class###############################################        