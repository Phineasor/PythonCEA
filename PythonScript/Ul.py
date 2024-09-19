# fmt: off
import InputValues as IV
from cea import runCEA, AxialValues
import math as m
import Injector as Inj
from EngineGeometry import LT, Dc, RatL
import mathfunctions as mf
import cantera as ct
import CoolProp.CoolProp as CP
from Bisect import Bisect

pa2atm = 9.86923 * 10 ** (-6)
in2m = 0.0254

#Combustion Paramaters
ceaOut = runCEA()
CombustionGas = ceaOut[0]
γ = CombustionGas.cp/CombustionGas.cv
cps = CombustionGas.cp

#Starting Paramaters
x = 0
dL = LT/IV.CellNum
Lold = 0
Rold = (Dc*in2m)/2
Tc = IV.FuelTankT
Tca = [0.0]*IV.CellNum


#Axial Values along engine
Values = AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)
Ts = Values[0]
ps = Values[1]
ρs = Values[2]
Tr = Values[3]
Ms = Values[4]

#Calculates Velocity Based on engine data in values and Mach number also does radii and lengthto save a loop
L = mf.linspace(0, LT, IV.CellNum)
Rad = [0.0]*IV.CellNum
Us = [0.0]*IV.CellNum
for i in range(IV.CellNum):
    Rad[i] = RatL(L[i])*in2m
    L[i] *= in2m
    Us[i] = Ms[i]*(γ*Ts[i]*(ct.gas_constant/CombustionGas.mean_molecular_weight))**0.5



BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
Ularray1 = [0.0]*IV.CellNum
Ularray2 = [0.0]*IV.CellNum
Ularray3 = [0.0]*IV.CellNum
for i in range(IV.CellNum):
    #Calculates Contour length and axial length
    dx = (((L[i]-Lold)**2)+(Rad[i]-Rold)**2)**0.5
    Lold = L[i]
    Rold = Rad[i]
    x += dx

    #This is an implicit formulation derrived from Grissom, Hopefully I can work this out for someone else:
    #This was pain and this still might eb the wrong function
    Tv = CP.PropsSI("T", "P", ps[i], "Q", 1, "Ethanol") #Saturation temp
    xe = 3.53*(Rad[i]*2)*(1+(x/(3.53*Rad[i])+(10**(-5)))**(-1.2))**(-1/1.2)
    Gch = ρs[i]*Us[i]
    Tm = 0.5*(Ts[i]+Tv)

    #Actuall function def for bisecting
    def Ulf(Ul):
        const = ((0.0592)/(2*0.023))**5
        a = (Gch*(Ts[i]/Tm)*((Us[i]-Ul)/Us[i]))**2
        b = const*(1/(xe*Rad[i]*2))
        return a - b
    Ul = Bisect(Ulf, 0, Us[i], 10*(-20))
    Ularray1[i]= Ul


    muc = CP.PropsSI("V", "T|liquid", IV.FuelTankT-15, "P", ps[i], "Ethanol")
    ρc = CP.PropsSI("D", "P", ps[i], "T|liquid", IV.FuelTankT, "Ethanol")

    #Austins Method
    def Ulf2(Ul):
        return m.sqrt((0.0592*BLCMdot*(Gch*Tr[i]/(Tm*Us[i]))**0.8*(Us[i]-Ul)**1.8)/(muc*ρc))-Ul
    Ul2 = Bisect(Ulf2, 0, Us[i], 10*(-20))
    #Ul2 = 1
    Ularray2[i] = Ul2

    Ga = BLCMdot/(2*Rad[i]*m.pi)
    mus = CombustionGas.viscosity


    def Gf(Ul):
        a = ρs[i]*Us[i]
        b = (Ts[i])/Tm
        c = (Us[i] - Ul)/Us[i]
        return a*b*c

    def Ulf3(Ul):
        a = Ga/ρs[i]
        b = (0.0592*ρc*Gf(Ul)*(Us[i]-Ul))/(mus*Ga)
        c = (Gf(Ul)*xe)/mus + (10**(-5))
        return (a*(b*c**(-0.2))**0.5) - Ul


    Ul3 = Bisect(Ulf3, 0, Us[i], 10*(-20))
    Ularray3[i] = Ul3

    #Ul = 10
    print(str(Us[i]) + ", " + str(Ul) + ", " + str(Ul2) + ", " + str(Ul3))
