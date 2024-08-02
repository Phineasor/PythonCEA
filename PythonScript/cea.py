# CEA calculations for rocket engines.
# fmt: off
"""
@author: phineas
"""

from EngineGeometry import RatL, LT
import InputValues as IV
from ChamberPressure import ChamberPressure
import Injector as Inj
import mathfunctions as mf
import cantera as ct
import CoolProp.CoolProp as CP

AxialDistances = mf.linspace(0, LT, IV.CellNum)
RadiusVal = []

i = 0
while i < len(AxialDistances):
    RadiusVal.append(RatL(AxialDistances[i]))
    i += 1

# print(AxialDistances)
# print("LineBreak")
# print(RadiusVal)

# mp.plot(AxialDistances, RadiusVal)
# mp.show()

Pc = IV.AmbP
# of = 1.4013 #Not defined here, not defined aywhere
Tol = 10 ** (-5)
RelError = 1

FuelMdot = IV.FuelOrificeNum * Inj.Mdot( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, 1931000, IV.FuelTankP*6894.76)
OxMdot = IV.OxOrificeNum * Inj.Mdot(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, 1931000, IV.OxTankP*6894.76)
Mdot = FuelMdot+OxMdot
print(Mdot)



i = 0
while (RelError > Tol) & (i < 500):
    i += 1  # This is here to ensure no infinite loops
    print("----------------------------------------------------------------------------------")
    PcOld = Pc #Keeps track of old Pc to fild reletive error

    #Creates CombustionGas
    CombustionGas = ct.Solution(IV.yaml)

    # Gets Mdot for Both fuel and Ox sides all orifaces, also total Mdot,
    FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP,)
    OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
    Mdot = FuelMdot+OxMdot

    OF = OxMdot / FuelMdot

    #Calculates reference Tempature and uses it to offset the enthalpy value for the CombustionGas to account for phase change
    TRef = max( CP.PropsSI("T", "P", Pc, "Q", 1, "Ethanol") * 1.01, CP.PropsSI("T", "P", Pc, "Q", 1, "O2") * 1.01,)
    CombustionGas.TPY = TRef, Pc, "O2:"+str(OF)+", C2H5OH:1" #Makes the CombustionGas have the correct OF ratio. and 
    CombustionGas.equilibrate("HP")  # I want to know why this is needed to prevent 0K -Phineas
    print(TRef)
    #Calculates and applys the enthalpy change used to account for phase change
    hCorrectionF = CP.PropsSI("H", "T", 285, "P", Pc, "Ethanol") - CP.PropsSI("H", "T", TRef, "P", Pc, "Ethanol")
    hCorrectionO = CP.PropsSI("H", "T", 85, "P", Pc, "O2") - CP.PropsSI("H", "T", TRef, "P", Pc, "O2")
    h = CombustionGas.h + (hCorrectionF * (1 / (1 + OF))) + (hCorrectionO * (OF / (1 + OF)))
    CombustionGas.HP = h, Pc

    #Calculates new ChamberPressure from the new tempature info
    gamma = CombustionGas.cp/CombustionGas.cv
    R = ct.gas_constant/CombustionGas.mean_molecular_weight
    print(CombustionGas.T)
    print(gamma)
    print(R)
    print(Mdot)
    Pc = ChamberPressure(CombustionGas.T, Mdot, gamma, R)

    #Calculates reletive error
    RelError = (abs(Pc-PcOld))/(Pc)
    print(Pc)
    print("----------------------------------------------------------------------------------")



print(CombustionGas())

print(Cp)
