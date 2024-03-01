import numpy as np
from matplotlib import pyplot as plt
import math
from scipy import optimize
gamma = 1.229
thrust_req = 250000 # Thrust requirment in Newtons
t0 = 3300 # (Kelvin)
p0 = 1500000 #15 Bar in (Pa)
pe = 101325 # Sea level pressure (Pa)
R = 397.8  # Calculated from total molar weight and universal gas constant

# Temperature at exit
te = t0 *(pe/p0)**((gamma -1)/gamma)
print("The temerature at the exit is:",te,"K")
#Mach number at exit calculation
Me = ((2/(gamma - 1)) *(((p0/pe)**((gamma -1)/(gamma))) -1))  ** 0.5
print("Mach number at the exit is:",Me)

#Mass flow rate calc

a = ( gamma * R * te) ** 0.5    #Speed of sound at exit
u = Me * a                      #Speed at exit

massflow = thrust_req/u # Massflow rate clauclated from thrust and speed at the exit 
print("Mass flow rate throught the nozzle is:", massflow,"kg/s")

# Area Calculation
dens0 = p0/(R * t0)  # Chamber desnity
print("Density in the chamber is:", dens0,"kg/m^3")
dens_throat = dens0 * ((2/(gamma+1)) ** (1/(gamma -1))) # Desntiy at the throat
print("Density at the throat is:", dens_throat,"kg/m^3")

a_throat = (gamma * R * t0 * (2/ (gamma + 1))) ** 0.5 #speed of sound at throat

Area_throat = massflow / (a_throat * dens_throat)
print("Area of the throat is:", Area_throat,"m^2")
rad_throat = (Area_throat/ math.pi ) ** 0.5
print ("The Radius of the thoat is ", rad_throat,"m")


# Nozzle exit area
Area_exhaust = (Area_throat * (1/Me)) * ((((2/(gamma+1)) * (1 + (( (gamma -1) /2)  * (Me **2)))) ** ( (gamma + 1) / (gamma -1) )) ** 0.5)
print ("The area of the exhasut is:", Area_exhaust,"m^2")
radius_exhaust = (Area_exhaust/math.pi) ** 0.5
print("The radius of the nozzle exit is:", radius_exhaust,"m")

massflow2 = dens_throat* Area_throat * ((gamma * R * ((2/ (gamma + 1)) * t0)) **0.5)
print("The mass flow rate is",massflow2,"kg/s") # Test for conversion of mass rule in compressible isentropic flow

# Design paremters for conical nozzle
R1 = rad_throat * 1.5

expansion_ratio = (radius_exhaust/rad_throat) ** 2
print ("The expansion ratio 'Epsilion'", expansion_ratio)

angle = 15
anglerad = math.radians(angle)
L1 = math.sin(anglerad) * R1

Rn =((1 - math.cos(anglerad)) * R1) + rad_throat

Ln = ((radius_exhaust - rad_throat) + (R1 * (math.cos(anglerad) -1))) / (math.tan(anglerad))

L = L1 + Ln
print ("Length of the nozzle is :",L,"m")

x = np.linspace(-L1 -0.150,L,150)
# Radius as a function of x
def y2(x):
    if (x<=-L1):
        return rad_throat + (R1 *(1 - (np.cos(np.arcsin(x/R1)))))
    elif (-L1<= x <= L1): #Curved scetion of nozzle between 0 and L1
        return rad_throat + (R1 *(1 - (np.cos(np.arcsin(x/R1)))))
    elif (x >= L1):       #Linear section
        return  ((radius_exhaust - Rn) / Ln) * (x - L1) + Rn 
    else:
        return 0
y3 = np.vectorize(y2)

plt.plot(x,y3(x))
plt.axhline(0,0, color='black', linewidth=0.5, linestyle='--')
plt.title('1D Conical Nozzle Geometry ')
plt.xlabel('Length (m)')
plt.ylabel('Radius (m)')
plt.ylim(rad_throat,radius_exhaust + 0.2)
plt.grid(True)
plt.show()

def mfromastar(x,A, Astar, gamma, M0):
    if x < 0:
        M0 = 0.5
    else:
        M0 = 3   
    def f(M):
        return (Astar/A)*(((M**2 * gamma - M**2 + 2)/(gamma + 1))**((gamma + 1)/(gamma - 1)))**0.5 - M
        
    return optimize.fsolve(f, M0)[0]

mfromastarvec = np.vectorize(mfromastar)

A = ((y3(x) **2) *  math.pi) # The radius' being converted to areas
Astar =  Area_throat
M = mfromastarvec(x,A, Astar, gamma, 3 )
plt.plot(x,M)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('Mach Number')
plt.show()

t = t0 * (1+((gamma -1)/2) * M **2)
def temp(M):
    return t0 /( (1+((gamma -1)/2) * M **2))   # Temperature Definition
temps = temp(M)

plt.plot(x, temps)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('Temperature (K)')
plt.show()

def dens(M):
    return dens0 /( (1+((gamma -1)/2) * M **2) **(1/(gamma -1))) #Density at any value
densitys = dens(M)
plt.plot(x, densitys)
plt.grid(True)
plt.xlabel('x')
plt.ylabel('Density (kg m^-3)')
plt.show()

def ps(M):
        return p0 /( (1+((gamma -1)/2) * M **2) **(gamma/(gamma -1)))
pressures = (ps(M))/100000 # Converted to Bar from Pascals for user readability

plt.plot(x, pressures )
plt.grid(True)
plt.xlabel('x')
plt.ylabel('Pressure (Bar)')
plt.show()