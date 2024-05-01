from Infrared_excess_class import *

l = np.linspace(0.02e-6, 100e-6, 10000) # Wavelength linspace
n = np.linspace(1e10, 1e17, 10000) # Frequency linspace

R_s = 0.013*R_sun # Stellar radius

a_i = 10*R_s  # Inner radius of debris disk
a_o = 28*R_s # Outter radius of debris disk
T_s = 11240 # Stellar temperature
val = [5.3999e17, 0] # Distance to the star and inclination in degrees
LF = "Flux" # Luminosity or Flux calculations
ln = "Wavelength" # Wavelength or Frequency calculations
Jansky = False #True Jy plot otherwise False
mJansky = True #True mJy plot otherwise False


xval = [l, n] # xvalues
xlim = [0.02e-6, 25e-6] # [xmin, xmax]
ylim = [0, 20] # [ymin, ymax]
log = False # True 'logged' plot otherwise False


IR = IR_excess(a_i, a_o, R_s, T_s, val, LF, ln, Jansky, mJansky) # Initializing class

func1 = IR.Stellar_Blackbody(xval = [l, n]) # Flux density of the stellar blackbody function
plot1 = IR.Plotting(xval, xlim, func1, ylim, 'Stellar blackbody', log) # Plotting of func1 as a function of wavelength in mJy

func2 = IR.Blackbody_Disk(xval = [l, n], T_d = 990) # Flux density of the Blackbody disk function
plot2 = IR.Plotting(xval, xlim, func2, ylim, r"Blackbody disk", log)  # Plotting of func2 as a function of wavelength in mJy


func3 = IR.Flat_Debris_Disk(xval = [l, n], alb = 0.2, pf = 0.7) # Flux density of the Flat Debris Disk function
plot3 = IR.Plotting(xval, xlim, func3, ylim, r'Flat disk', log)  # Plotting of func3 as a function of wavelength in mJy

plt.legend(fontsize = 20) # creates a legend

plt.show() # Shows plots
