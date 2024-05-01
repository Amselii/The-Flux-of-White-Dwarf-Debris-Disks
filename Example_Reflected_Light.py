from Reflected_light_class import *

l = np.linspace(0.02e-6, 1000e-6, 100000) # Wavelength linspace
n = np.linspace(1e12, 1e15, 100) # Frequency linspace


R_s = 0.013*R_sun # Stellar radius
a_i = 15*R_s  # Inner radius of debris disk
a_o = 23*R_s # Outter radius of debris disk
T_s = 13060 # Stellar temperature
d   = 9.089789e17 # Distance

A = 0.2 # Albedo
f = 0.47 # Packing fraction
inc = 0 # Inclination in degree
ln = "Wavelength" # Wavelength or Frequency calculations
mJansky = False # True mJy plot otherwise False

xval = [l, n] #x-values
xlim = [0.02e-6, 12e-6] # [xmin, xmax]
ylim = [1e-11, 1e-6] # [ymin, ymax]
log = True # True logged plot otherwise False

F = Reflected_light(R_s, T_s, a_i, a_o, d, A, f, inc, ln, mJansky) # Initializing class

funcS = F.Stellar_blackbody(xval = [l, n]) # Flux density of the stellar blackbody
plotS = F.Plotting(xval, xlim, funcS, ylim, 'Stellar blackbody',log) # Plotting log-log plot of funcS as a function of wavelength

func = F.Lambertian(xval = [l, n]) # Flux density of the Lambertian surface
plot = F.Plotting(xval, xlim, func, ylim, "Lambertian") # Plotting log-log plot of func as a function of wavelength

#HG - Henyey-Greenstein
func2 = F.HG(xval = [l, n], g = 0.5) # Flux density of the reflected light with HG phase function
plot2 = F.Plotting(xval, xlim, func2, ylim, "HG") # Plotting log-log plot of func2 as a function of wavelength

plt.legend(loc = "upper right", fontsize = 20) # Creates a legend
plt.show()