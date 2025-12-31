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

Pc = IV.AmbP 
Tc = IV.AmbT
Tol = 10 ** (-5)
RelError = 1
Damp = 0.9

#Creates CombustionGas
GOx = ct.Solution(IV.yaml)
#print(CombustionGas.report())

i = 0
while (RelError > Tol) & (i < 250):
    i += 1  # This is here to ensure no infinite loops
    PcOld = Pc #Keeps track of old Pc to find reletive error
    
    GOxT = CP.PropsSI("T", "P", Pc, "Q", 1, "O2")
    GOx.TPY = GOxT, Pc, "O2:1"

    # Gets Mdot for Both fuel and Ox sides all orifices, also total Mdot,
    FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
    OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
    Mdot = FuelMdot+OxMdot

    OF = OxMdot/FuelMdot
    OFV = OF*(CP.PropsSI("D", "P", Pc, "Q", 0, "Ethanol")/GOx.density)

    Tavg = IV.FuelTankT*(1/(1+OFV)) + GOxT*(OFV/(1+OFV))

    γMavg = (1/(1+OF)) + (GOx.cp/GOx.cv)*(OF/(1+OF))

    MMW = 46*(1/(1+OF)) + (GOx.mean_molecular_weight)*(OF/(1+OF))
    R = ct.gas_constant/MMW



    Pc = ChamberPressure(Tavg, Mdot, γMavg, R)
    #print((γMavg*R*Tc)**0.5)
 
    #Calculates Max reletive error between ChamberPressure and ChamberTempature
    RelError = (abs(Pc-PcOld))/(Pc)

    #Calculates the chamge in ChamberPressure and ChamberTempature, makes sure its not so large it just overshoots everything and the engine "explodes"
    Pc -= Damp*(Pc-PcOld)
    #print(GOx.report())
    #print(str(GOx.density) + ", " + str(CP.PropsSI("D", "P", Pc, "Q", 0, "O2")))
    #print(str(Tavg) + ", " + str(GOxT) + ", " + str(IV.FuelTankT))
    #print(Pc)

    i += 1
#print(i)
Davg = (CP.PropsSI("D", "P", Pc, "Q", 0, "Ethanol"))*(1/(1+OF)) + (GOx.density)*(OF/(1+OF))
print(R)
print(Tavg)
print(γMavg)
print(Mdot)
    