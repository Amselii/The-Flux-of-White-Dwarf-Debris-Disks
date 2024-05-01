import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

R_sun = 695700000 #m
L_sun = 3.828e26 #W

h = 6.62607015e-34 #Js
c = 299792458 #m/s
k_B = 1.380649e-23
sigma = 5.67037442e-08

fig, ax = plt.subplots()

class Reflected_light:
    
    def __init__(self, R_s, T_s, a_i, a_o, d, A, f, inc, ln = None, mJansky = False):
        
        if not isinstance(R_s, (float, int)):
            raise TypeError("R_s has to be of type float or int.")
        if not isinstance(T_s, (float, int)):
            raise TypeError("T_s has to be of type float or int.")
        if not(ln == 'Wavelength' or ln == 'Frequency'):
            raise Exception("Please enter 'Wavelength' or 'Frequency' as ln.")
        if not isinstance(mJansky, bool):
            raise TypeError("mJansky has to be of type bool.")
        
        self.R_s = R_s
        self.T_s = T_s
        self.a_i = a_i
        self.a_o = a_o
        
        self.d = d
        self.inc = inc
        
        self.A = A
        self.f = f
        
        self.ln = ln
        self.mJansky = mJansky
        
    
    def Stellar_blackbody(self, xval):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")
        
        mJy = np.array([1 if self.mJansky == False else 1e29])[0]
        
        if self.ln == "Wavelength":
            lamda = xval[0]
            
            B_l = 2*h*(c**2)/(lamda**5)*1/(np.exp(h*c/(k_B*self.T_s*lamda))-1)
        
            return np.pi/self.d**2*B_l*self.R_s**2*mJy
        
        if self.ln == "Frequency":
            nu = xval[1]
            
            B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*self.T))-1)
            
            return np.pi/self.d**2*B_n*self.R_s**2*mJy
        
    
    def Lambertian(self, xval):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")
        
        mJy = np.array([1 if self.mJansky == False else 1e29])[0]
    
        
        if self.ln == "Wavelength":
            lamda = xval[0]
            
            B_l = 2*h*(c**2)/(lamda**5)*1/(np.exp(h*c/(k_B*self.T_s*lamda))-1)
            
            return self.A*self.f*4/3*(self.R_s)**3*B_l*(1/self.a_i-1/self.a_o)*1/self.d**2*mJy*np.cos(np.deg2rad(self.inc))

        if self.ln == "Frequency":
            nu = xval[1]
            
            B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*self.T))-1)
            
            return self.A*self.f*4/3*(self.R_s)**3*B_n*(1/self.a_i-1/self.a_o)*1/self.d**2*mJy*np.cos(np.deg2rad(self.inc))
    
    
    def HG(self, xval, g):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")
        
        mJy = np.array([1 if self.mJansky == False else 1e29])[0]
        
        if self.ln == "Wavelength":
            lamda = xval[0]
            
            I = []
            
            B_l = 2*h*(c**2)/(lamda**5)*1/(np.exp(h*c/(k_B*self.T_s*lamda))-1)
    
            def Phase_function(psi):
        
                C = np.sin(psi)*np.sin(np.deg2rad(self.inc))
        
                return 1/(4*np.pi)*(1-g**2)*1/(1+g**2-2*g*C)**(3/2)
            
            I.append(sc.integrate.quad(Phase_function, 0, 2*np.pi)[0])
            
            return self.A*self.f*2/(3)*B_l*self.R_s**3*(1/self.a_i-1/self.a_o)*I/self.d**2*mJy
        
        if self.ln == "Frequency":
            nu = xval[1]
            
            I = []
            
            B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*self.T))-1)
            
            def Phase_function(psi):
        
                C = np.sin(psi)*np.sin(np.deg2rad(self.inc))
        
                return 1/(4*np.pi)*(1-g**2)*1/(1+g**2-2*g*C)**(3/2)
            
            I.append(sc.integrate.quad(Phase_function, 0, 2*np.pi)[0])
            
            return self.A*self.f*2/(3)*B_n*self.R_s**3*(1/self.a_i-1/self.a_o)*I/self.d**2*mJy
        
        
    def Plotting(self, xval, xlim, func, ylim, label, log = False):
        
        if not isinstance(xval, list):
             raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
             raise ValueError("xval must contain exactly two variables.")
        if not isinstance(xlim, list):
            raise TypeError("xlim has to be of type list.")
        if len(xlim) != 2:
            raise ValueError("xlim must contain exactly two variables.")
        if not isinstance(ylim, list):
            raise TypeError("ylim has to be of type list.")
        if len(ylim) != 2:
            raise ValueError("ylim must contain exactly two variables.")
        if not isinstance(label, str):
            raise TypeError("label has to be of type str.")
        if not isinstance(log, bool):
            raise TypeError("log has to be of type bool.")

      
        if log == True:
            ax.set_xscale("log")
            ax.set_yscale("log")
            
        if self.mJansky == False:
            
            if self.ln == "Wavelength":               
                
                xVal = xval[0]
                
                formatter = lambda x, pos: f'{(x / 10**(-6)):0.0000001f}'
                
                ax.xaxis.set_major_formatter(formatter)
                ax.set_xlabel(r"$\lambda$ ($\mu$m)", fontsize = 22)

                ax.set_ylabel(r"$F (W/m^{3})$ ", fontsize = 22)
            
            if self.ln == "Frequency":
                ax.set_xlabel(r"$\nu$ (Hz)", fontsize = 22)
                
                xVal = xval[1]
 
                ax.set_ylabel(r"$F \, (W \cdot m^{-2} Hz^{-1})$ ", fontsize = 22) 
        
        else:
            
            if self.ln == "Wavelength":
                
                formatter = lambda x, pos: f'{(x / 10**(-6)):0.0000001f}'
                
                ax.xaxis.set_major_formatter(formatter)
                ax.set_xlabel(r"$\lambda$ ($\mu$m)", fontsize = 22)
                
                xVal = xval[0]
            
            else:
                ax.set_xlabel(r"$\nu$ (Hz)", fontsize = 22)
                xVal = xval[1]
            
            ax.set_ylabel(r"$F \, (mJy)$ ", fontsize = 22)
        
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])
        plt.xticks(fontsize = 20)
        plt.yticks(fontsize = 20)
        
            
        plot = plt.plot(xVal, func, label=f'{label}')   
        
        return plot