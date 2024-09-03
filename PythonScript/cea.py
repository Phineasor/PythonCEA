# CEA calculations for rocket engines.
# fmt: off
"""
@author: phineas
"""

from EngineGeometry import RatL, LT, Lt
import InputValues as IV
from ChamberPressure import ChamberPressure
import Injector as Inj
import mathfunctions as mf
import cantera as ct
import CoolProp.CoolProp as CP
import math
from Bisect import Bisect

in2m = 0.0254 #inch 2 meter conversion, should be moved to a different file
psi2pa = 6894.71 #Psi to pascals


def runCEA():

    #Starting conditions for iteration
    Pc = IV.AmbP 
    Tc = IV.AmbT
    Tol = 10 ** (-5)
    RelError = 1
    Damp = 0.9

    #Creates CombustionGas
    CombustionGas = ct.Solution(IV.yaml)

    i = 0
    while (RelError > Tol) & (i < 250):
        i += 1  # This is here to ensure no infinite loops
        PcOld = Pc #Keeps track of old Pc to find reletive error
        TcOld = Tc #Keeps track of old TC to find reletive error


        # Gets Mdot for Both fuel and Ox sides all orifices, also total Mdot,
        FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
        OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
        Mdot = FuelMdot+OxMdot

        #Calculates OF ratio, technically not efficient to have it here or caculated this way, but eh.
        OF = OxMdot / FuelMdot
        #OF = 3.5
        #Calculates reference Tempature and uses it to offset the enthalpy value for the CombustionGas to account for phase change
        TRef = max( CP.PropsSI("T", "P", Pc, "Q", 1, "Ethanol"), CP.PropsSI("T", "P", Pc, "Q", 1, "O2"))*1.01
        CombustionGas.TPY = TRef, Pc, "O2:"+str(OF)+", C2H5OH:1" #Makes the CombustionGas have the correct OF ratio. and 
        CombustionGas.equilibrate("HP")  # I want to know why this is needed to prevent 0K -Phineas

        #Calculates and applys the enthalpy change used to account for phase change
        hCorrectionF = CP.PropsSI("H", "T", IV.FuelTankT, "P", Pc, "Ethanol") - CP.PropsSI("H", "T", TRef, "P", Pc, "Ethanol")
        hCorrectionO = CP.PropsSI("H", "T", IV.OxTankT, "P", Pc, "O2") - CP.PropsSI("H", "T", TRef, "P", Pc, "O2")
        h = CombustionGas.h + (hCorrectionF * (1 / (1 + OF))) + (hCorrectionO * (OF / (1 + OF)))
        #h = CombustionGas.h
        CombustionGas.HP = h, Pc
        CombustionGas.equilibrate("HP")

        #Calculates new ChamberPressure from the new tempature info, and gas propertys
        gamma = CombustionGas.cp/CombustionGas.cv
        R = ct.gas_constant/CombustionGas.mean_molecular_weight
        Tc = CombustionGas.T
        Pc = ChamberPressure(Tc, Mdot, gamma, R)
 
        #Calculates Max reletive error between ChamberPressure and ChamberTempature
        RelError = max((abs(Pc-PcOld))/(Pc),(abs(Tc-TcOld))/(Tc))

        #Calculates the chamge in ChamberPressure and ChamberTempature, makes sure its not so large it just overshoots everything and the engine "explodes"
        Pc -= Damp*(Pc-PcOld)
        Tc -= Damp*(Tc-TcOld)

    return CombustionGas, Mdot

#This function fionds the axial values for several things, temp pressure adiabatic wall temp, etc
def AxialValues():


CEAout = runCEA()
CombustionGas = CEAout[0]
Mdot = CEAout[1]
R = ct.gas_constant/CombustionGas.mean_molecular_weight
gamma = CombustionGas.cp/CombustionGas.cv
Tc = CombustionGas.T
Pc = ChamberPressure(Tc, Mdot, gamma, R)

BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)

print(CombustionGas())
#print(str(Pc/Inj.PSI2PA)+" : "+str(Pc))



#Crates array of Lengths and an array of radisus
AxialDistances = mf.linspace(0, LT, IV.CellNum)
RadiusVal = [0.0]*IV.CellNum

for i in range(IV.CellNum):
    RadiusVal[i]=RatL(AxialDistances[i])

#Calculaions for Engine conditions along length
Mach = [0.0]*IV.CellNum
Ts = [0.0]*IV.CellNum
Taw = [0.0]*IV.CellNum
Tr = [0.0]*IV.CellNum
P = [0.0]*IV.CellNum 


for i in range(IV.CellNum):
    A = math.pi*RadiusVal[i]**2
    def MachNumber(x):
        At = 0.25*math.pi*IV.Dt**2
    
        exp = (gamma+1)/(gamma-1)
        mult = 0.5*(gamma-1)

        return ((1/x)*(((1+(mult*x**2))/(1+mult))**exp)**0.5)-(A/At)
    
    if AxialDistances[i] < Lt:
        Mach[i] = Bisect(MachNumber, 0.001, 0.99999, 10**-20)
    else:
        Mach[i] = Bisect(MachNumber, 1.00001, 5, 10**-20)

    Ts[i] = Tc*(1+(0.5*(gamma-1))*Mach[i]**2)**(-1)
    P[i] = Pc*(1+((gamma-1)/2)*Mach[i]**2)**(-(gamma)/(gamma-1))

    Pr = (CombustionGas.cp*CombustionGas.viscosity)/(CombustionGas.thermal_conductivity) 
    r = Pr**0.33
    mult = (0.5*(gamma-1))*Mach[i]**2
    Taw[i] = Tc*((1+r*(mult))/(1+mult))

    #print(Taw[i])
    #print(P[i]*(9.869*10**(-6)))

#print("Pr = "+str(Pr))
#print("mu = "+str(CombustionGas.viscosity))
#print("Cp = "+str(CombustionGas.cp))
#print(CombustionGas.thermal_conductivity)
#print("Ox = "+str(264.172*10*OxMdot/CP.PropsSI("D", "T|liquid", IV.OxTankT, "P", IV.OxTankP*psi2pa, "O2")))
#print("Fuel = "+str(264.172*10*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
#print("Fuel gal/s = "+str(264.172*(FuelMdot+BLCMdot)/CP.PropsSI("D", "T|liquid", IV.FuelTankT, "P", IV.FuelTankP*psi2pa, "Ethanol")))
