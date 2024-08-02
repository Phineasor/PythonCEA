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

Pc = IV.AmbP #Starting ChamberPressure can be a guess
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


    # Gets Mdot for Both fuel and Ox sides all orifaces, also total Mdot,
    FuelMdot = IV.FuelOrificeNum * Inj.MdotSPIONLY( IV.FuelOrificeCd, IV.FuelOrificeDiameter, IV.Fuel, IV.FuelTankT, Pc, IV.FuelTankP)
    OxMdot = IV.OxOrificeNum * Inj.MdotSPIONLY(IV.OxOrificeCd, IV.OxOrificeDiameter, IV.Ox, IV.OxTankT, Pc, IV.OxTankP)
    Mdot = FuelMdot+OxMdot

    #Calculates OF ratio, technically not efficient to have it here or caculated this way, but eh.
    OF = OxMdot / FuelMdot

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
    gamma2 = CombustionGas.cp_mass/CombustionGas.cv_mass
    R = ct.gas_constant/CombustionGas.mean_molecular_weight
    Tc = CombustionGas.T
    Pc = ChamberPressure(Tc, Mdot, gamma, R)
 
    #Calculates Max reletive error between ChamberPressure and ChamberTempature
    RelError = max((abs(Pc-PcOld))/(Pc),(abs(Tc-TcOld))/(Tc))

    #Calculates the chamge in ChamberPressure and ChamberTempature, makes sure its not so large it just overshoots everything and the engine "explodes"
    Pc -= Damp*(Pc-PcOld)
    Tc -= Damp*(Tc-TcOld)

print(CombustionGas())
print(Tc)
print(Pc/Inj.PSI2PA)
