# The Flux of White Dwarf Debris Disks

The two classes IR_excess and Reflected_light found in Infrared_excess_class.py and Reflected_light_class.py respectively, are for plotting the flux density
of white dwarf debris disks.

# About code
## Parameters

a_i: is a float or an integer and describes the inner edge of the debris disk.

a_o: is a float or an integer and describes the outer edge of the debris disk.

R_S: is a float or an integer and describes the stellar radius.

T_S: is a float or an integer and describes the stellar temperature in Kelvin.

val: is a list containing exactly 2 elements, namely the distance to the white dwarf and inclination in degrees.

LF: can either be 'Luminosity' or 'Flux' depending on if one wants to have the functions calculated for luminosity or flux density.

ln: can either be 'Wavelength' or 'Frequency' depending on if one wants to have the functions calculated in wavelength or frequency.

Jansky: is a Boolean and choosing Jansky all functions get calculated in Jansky instead of SI units. It is set to False by default.

mJansky: is a Boolean and choosing mJansky all functions get calculated in milli Jansky instead of SI units. It is set to False by default.

xval: is a list containing exactly 2 elements, namely an array for wavelength and one array for frequency.

alb: is a float or an integer and describes the albedo of the disk. It is set to 0 by default.

pf: is a float or an integer and describes the packing fraction of the disk. It is set to 1 by default.

xlim: is a list containing exactly 2 elements, namely the xmin and xmax for the plot.

func: is the function that is plotted.

ylim: is a list containing exactly 2 elements, namely the ymin and ymax for the plot.

label: is a string for the label of the plot.

log: is a Boolean and the plot is logged on both axes for log being True. It is set to False by default.


d: is a float or an integer and describes the distance to the white dwarf.

A: is a float or an integer and describes the albedo.

f: is a float or an integer and describes the packing fraction.

inc: is a float or an integer and describes the inclination of the disk in degrees.

g: is a float or an integer and describes the asymmetry of the phase function.

## Infrared excess class

This class has three models for studying the infrared excess of a white dwarf debris disk. All functions can calculate the infrared excess in either 
flux density or luminosity for either wavelength or frequency in either SI units, Jansky or milli Jansky. The first infrared excess function
```Py
Stellar_Blackbody(self, xval)
```
calculates the flux density or luminosity of the white dwarf assuming a blackbody. The second infrared excess function
```py
Blackbody_Disk(self, xval, T_d)
``` 
calculates the flux density or luminosity of the debris disk assuming it to be a simple blackbody at a given temperature. The last infrared excess function
```py
Flat_Debris_Disk(self, xval, alb = None, pf = None)
```
calculates the flux density or luminosity of the debris disk using the flat disk model. The last function in the class 
```py
Plotting(self, xval, xlim, func, ylim, label, log = False)
```
is the plotting function, which takes the x- and y-limits, one of the three infrared excess functions, and a label. It can plot the infrared excess in either 
a logged or non-logged plot.

## Reflected light

This class has three models for studying the reflected light coming from a white dwarf debris disk. All functions can calculate the reflected flux density for 
either wavelength or frequency in either SI units or milli Jansky. The first reflected light function
```py
Stellar_Blackbody(self, xval)
```
is the same as for the infrared excess class and calculates the flux density for the white dwarf assuming a blackbody. The second reflected light function
```py
Lambertian(self, xval)
```
calculates the flux density of the debris disk assuming a Lambertian surface. The last reflected light function
```py
HG(self, xval, g)
```
calculates the flux density of the debris disk using the Henyey-Greenstein phase function, which has a parameter g describing the asymmetry of the phase function.
The last function in the class 
```py
Plotting(self, xval, xlim, func, ylim, label, log = False)
```
is the plotting function, which takes the x- and y-limits, one of the three reflected light functions and a label. It can plot the reflected light in either 
a logged or non-logged plot.

# How to use the program

To use the program one has to import the class as such
```py
from Infrared_excess_class import *
```
and then initialize the class.
```py
IR_excess(a_i, a_o, R_s, T_s, val, LF, ln, Jansky, mJansky)
```
The same steps have to be followed for the reflected light class. For variable definitions, one can look at the example files.
