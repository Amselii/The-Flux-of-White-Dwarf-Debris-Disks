import numpy as np
import scipy as sc
import matplotlib.pyplot as plt

R_sun = 695700000 #m
L_sun = 3.828e26 #W

h = 6.62607015e-34 #Js
c = 299792458 #m/s
k_B = 1.380649e-23

fig, ax = plt.subplots()
class IR_excess:
    
    def __init__(self, a_i, a_o, R_s, T_s, val, LF = None, ln = None, Jansky = False, mJansky = False):
        
        if not isinstance(a_i, (float, int)):
            raise TypeError("a_i has to be of type float or int.")
        if not isinstance(a_o, (float, int)):
            raise TypeError("a_o has to be of type float or int.")
        if not isinstance(R_s, (float, int)):
            raise TypeError("R_s has to be of type float or int.")
        if not isinstance(T_s, (float, int)):
            raise TypeError("T_s has to be of type float or int.")
        if not isinstance(val, list):
            raise TypeError("val has to be of type list containing the distance d and inclination inc.")
        if len(val) != 2:
            raise ValueError("val must contain exactly two variables.")
        if not(LF == 'Luminosity' or LF == 'Flux'):
            raise Exception("Please enter 'Luminosity' or 'Flux' as LF.")
        if not(ln == 'Wavelength' or ln == 'Frequency'):
            raise Exception("Please enter 'Wavelength' or 'Frequency' as ln.")
        if not isinstance(Jansky, bool):
            raise TypeError("Jansky has to be of type bool.")
        if not isinstance(mJansky, bool):
            raise TypeError("mJansky has to be of type bool.")
        
        self.a_i = a_i
        self.a_o = a_o
        self.R_s = R_s
        self.T_s = T_s
        self.LF = LF
        self.ln = ln
        self.Jansky = Jansky
        self.mJansky = mJansky
    
        if self.LF == "Flux":
            self.val = val
            self.d = val[0]
            self.inc = val[1]
    
    def Stellar_Blackbody(self, xval):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")
                
        if self.Jansky == False and self.mJansky == False:
            
            T = self.T_s
            
            if self.ln == "Wavelength":
                lamda = xval[0]
                
                B_l = 2*h*(c**2)/(lamda**5)*1/(np.exp(h*c/(k_B*T*lamda))-1)
                
                if self.LF == "Luminosity":
                    
                    return 4*np.pi**2*B_l*lamda*self.R_s**2
                
                if self.LF == "Flux":
                    
                    return np.pi/self.d**2*B_l*self.R_s**2 
            
            if self.ln == "Frequency":
                nu = xval[1]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                
                if self.LF == "Luminosity":
                    
                    return 4*np.pi**2*B_n*nu*self.R_s**2
                
                if self.LF == "Flux":
                    
                    return np.pi/(self.d**2)*B_n*self.R_s**2 
        else:
            
            Jy = np.array([1e26 if self.mJansky == False else 1e29])[0]

            T = self.T_s
            
            if self.ln == "Wavelength":
                nu = c/xval[0]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                    
                return np.pi/self.d**2*B_n*self.R_s**2*Jy 
            
            if self.ln == "Frequency":
                nu = xval[1]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                
                return np.pi/(self.d**2)*B_n*self.R_s**2*Jy
    
    def Blackbody_Disk(self, xval, T_d):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")

                
        if self.Jansky == False and self.mJansky == False:
            
            T = T_d
            
            if self.ln == "Wavelength":
                lamda = xval[0]
                
                B_l = 2*h*(c**2)/(lamda**5)*1/(np.exp(h*c/(k_B*T*lamda))-1)
                
                if self.LF == "Luminosity":
                    
                    return 4*np.pi**2*B_l*lamda*(self.a_o**2-self.a_i**2)
                
                if self.LF == "Flux":
                    
                    return np.pi*np.cos(np.deg2rad(self.inc))/self.d**2*B_l*(self.a_o**2-self.a_i**2)
            
            if self.ln == "Frequency":
                nu = xval[1]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                
                if self.LF == "Luminosity":
                    
                    return 4*np.pi**2*B_n*nu*(self.a_o**2-self.a_i**2)
                
                if self.LF == "Flux":
                    
                    return np.pi*np.cos(np.deg2rad(self.inc))/(self.d**2)*B_n*(self.a_o**2-self.a_i**2)
            
        else:
            
            Jy = np.array([1e26 if self.mJansky == False else 1e29])[0]

            T = T_d
            
            if self.ln == "Wavelength":
                nu = c/xval[0]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                    
                return np.pi*np.cos(np.deg2rad(self.inc))/(self.d**2)*B_n*(self.a_o**2-self.a_i**2)*Jy 
            
            if self.ln == "Frequency":
                nu = xval[1]
                
                B_n = 2*h*nu**3/(c**2)*1/(np.exp(h*nu/(k_B*T))-1)
                
                return np.pi*np.cos(np.deg2rad(self.inc))/(self.d**2)*B_n*(self.a_o**2-self.a_i**2)*Jy
        

        
    def Flat_Debris_Disk(self, xval, alb = None, pf = None):
        
        if not isinstance(xval, list):
            raise TypeError("xval has to be of type list.")
        if len(xval) != 2:
            raise ValueError("xval must contain exactly two variables.")

        Albedo = np.array([1 if alb == None else (1-alb)])
        Emissivity = 0.9
        Packing_Fraction = np.array([1 if pf == None else pf])
        Jy = np.array([1e26 if self.mJansky == False else 1e29])[0]

        I = []        
        
        if self.Jansky == False and self.mJansky == False:
        
            if self.ln == "Wavelength":
                lamda = xval[0]
                xVal = lamda
                
                if self.LF == "Luminosity":   
            
                    C = 8*np.pi**2*lamda
                
                if self.LF == "Flux":
                    
                    C = 2*np.pi*np.cos(np.deg2rad(self.inc))/self.d**2
                    
            if self.ln == "Frequency":
                nu = xval[1]
                xVal = nu
                
                if self.LF == "Luminosity":
                
                    C = 8*np.pi**2*nu
                
                if self.LF == "Flux":
                    
                    C = 2*np.pi*np.cos(np.deg2rad(self.inc))/self.d**2
            
            for i in xVal:        
                def f(a):
                    alpha = 4/(3*np.pi*a)*self.R_s
                
                    T = ((alpha/2)**(1/4)*(self.R_s/a)**(1/2)*self.T_s)*(Albedo/Emissivity)**(1/4)
                
                    if self.ln == "Wavelength":    
                    
                        B_l = (2*h*(c**2)/(i**5)*1/(np.exp(h*c/(k_B*T*i))-1))
                        
                        return B_l*a
                    
                    if self.ln == "Frequency":
                        
                        B_n = 2*h*i**3/(c**2)*1/(np.exp(h*i/(k_B*T))-1)
                        
                        return B_n*a
                        
                
                I.append(sc.integrate.quad(f, self.a_i, self.a_o)[0])
        
            return C * np.array(I)*Packing_Fraction
        
        else:
            
            C = 2*np.pi*np.cos(np.deg2rad(self.inc))/self.d**2*Jy
            
            if self.ln == "Wavelength":
                nu = c/xval[0]
                xVal = nu
                    
            if self.ln == "Frequency":
                nu = xval[1]
                xVal = nu
                
            for i in xVal:        
                def f(a):
                    alpha = 4/(3*np.pi*a)*self.R_s
                
                    T = ((alpha/2)**(1/4)*(self.R_s/a)**(1/2)*self.T_s)*(Albedo/Emissivity)**(1/4)
                        
                    B_n = 2*h*i**3/(c**2)*1/(np.exp(h*i/(k_B*T))-1)
                        
                    return B_n*a
                        
                
                I.append(sc.integrate.quad(f, self.a_i, self.a_o)[0])
        
            return C * np.array(I)*Packing_Fraction
    
        
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
            
        if self.Jansky == False and self.mJansky == False:
            
            if self.ln == "Wavelength":               
                
                xVal = xval[0]
                
                formatter = lambda x, pos: f'{(x / 10**(-6)):0.0000001f}'
                
                ax.xaxis.set_major_formatter(formatter)
                ax.set_xlabel(r"$\lambda$ ($\mu$m)")
                
    
                if self.LF == "Luminosity":
                    ax.set_ylabel(r"$L_{\lambda} ( \, W)$ ")
                else:
                    ax.set_ylabel(r"$F_{\lambda} (W/m^{3})$ ")
            
            if self.ln == "Frequency":
                ax.set_xlabel(r"$\nu$ (Hz)")
                
                xVal = xval[1]
                
                if self.LF == "Luminosity":
                    ax.set_ylabel(r"$L_{\nu} (\, W)$ ")                
                else: 
                    ax.set_ylabel(r"$F_{\nu} \, (W \cdot m^{-2} Hz^{-1})$ ") 
        
        else:
            
            if self.ln == "Wavelength":
                
                formatter = lambda x, pos: f'{(x / 10**(-6)):0.0000001f}'
                
                ax.xaxis.set_major_formatter(formatter)
                ax.set_xlabel(r"$\lambda$ ($\mu$m)")
                
                xVal = xval[0]
            
            else:
                ax.set_xlabel(r"$\nu$ (Hz)")
                xVal = xval[1]
            
            ax.set_ylabel(r"$F \, (Jy)$ ")
            if self.mJansky == True:
                ax.set_ylabel(r"$F \, (mJy)$")
        
        ax.set_xlim(xlim[0], xlim[1])
        ax.set_ylim(ylim[0], ylim[1])

        plot = plt.plot(xVal, func, label=f'{label}')   
        
        return plot
