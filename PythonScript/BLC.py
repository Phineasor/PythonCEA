# No like truly this is a lot
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

#This function calculates the emittance of the gas for some point in the rocket engine based on the diameter and the CombusionGas
def emittance(CombusionGas, R):
    # effective lenght through gas, this is an approx given in grissom, will probably be updated eventually.
    Leff = 0.95 * 2 * R * IV.Aw ** (-0.85)

    # Partial Pressures calculated form mole fraction and total Pressures
    χC = CombusionGas.mole_fraction_dict()["CO2"]
    χH = CombusionGas.mole_fraction_dict()["H2O"]
    ppC = χC * CombusionGas.P * pa2atm
    ppH = χH * CombusionGas.P * pa2atm

    # Optical Densities Based on the effective length through the gas and the partial pressure, should probably be a path integral
    ρOptC = ppC * Leff
    ρOptH = ppH * Leff

    T = CombusionGas.T
    # These are the 4 tempature based parabolas from the table in Grisson for the c and n coefficients this should probably be replaced with a general polyfit function
    cC = (2.5 * 10 ** (-8)) * T**2 + (-0.00005) * T + (0.075)
    nC = (0) * T**2 + (0) * T + (0.6)
    cH = (2.075 * 10 ** (-7)) * T**2 + (0.0001125) * T + (-0.155)
    nH = (-1.2 * 10 ** (-7)) * T**2 + (0.00056) * T + (0.01)

    # Calculates uncorrecected emittance values without the Kp pressure correction
    εfC = 0.231
    εfH = 0.825
    εC = εfC * (1 + (ρOptC / cC) ** (-nC)) ** (-1 / nC)
    εH = εfH * (1 + (ρOptH / cH) ** (-nH)) ** (-1 / nH)

    # Pressure correction factors for the given emittance values at 1 atm
    KpC = 10 ** ( 0.036 * ρOptC ** (-0.489) * (1 + (2 * m.log(CombusionGas.P*pa2atm, 10)) ** (-100 * ppC)) ** (-1 / (100 * ppC)))

    c1 = 0.26 + 0.74 * m.exp(-2.5 * ρOptH)
    c2 = 0.75 + 0.31 * m.exp(-10 * ρOptH)
    KpH = 1 + c1 * ( 1 - (m.exp((1 - CombusionGas.P*pa2atm * (1 + χH)) / c2)))

    # Calculates the emittance correction factor to handle ovlerlap in the 2 specturms
    n = 5.5 * (1 + (1.09 * (ρOptC + ρOptH)) ** (-3.88)) ** (-1 / 3.88)
    Kx = 1 - (abs((2 * χH)/( χH + χC) - 1)) ** n
    Δε = 0.0551*Kx*(1-m.exp(-4*(ρOptC+ρOptH)))*(1-m.exp(-12.5*(ρOptH+ρOptH)))

    # Return the final emittance (Probably)
    return (KpC*εC) + (KpH*εH) - (Δε)

def getUl():    
    return 1


def calcBLC():
    #Starting Paramaters
    Tbl = IV.FuelTankT
    x = 0
    dL = LT/IV.CellNum
    Lold = 0
    Rold = (Dc*in2m)/2

    #Combustion Paramaters
    ceaOut = runCEA()
    CombustionGas = ceaOut[0]
    γ = CombustionGas.cp/CombustionGas.cv

    #Axial Values along engine
    Values = AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)
    Ts = Values[0]
    ps = Values[1]
    ρs = Values[2]
    Tr = Values[3]
    Ms = Values[4]

    L = mf.linspace(0, LT, IV.CellNum)
    Rad = [0.0]*IV.CellNum

    #Creats array of radii and lengths in meters
    for i in range(IV.CellNum):
        Rad[i] = RatL(L[i])*in2m
        L[i] *= in2m

    #Calculates Velocity Based on engine data in values and Mach number also does radii to save a loop
    Us = [0.0]*IV.CellNum
    for i in range(IV.CellNum):
        Us[i] = Ms[i]*(γ*Ts[i]*(ct.gas_constant/CombustionGas.mean_molecular_weight))**0.5 

    BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
    for i in range(IV.CellNum):
        #Calculates Contour length and axial length
        dx = (((L[i]-Lold)**2)+(Rad[i]-Rold)**2)**0.5
        Lold = L[i]
        Rold = Rad[i]
        x += dx
        #Coolant Flow length per circumfrence
        #Γ = BLCMdot/



        #working Probably
        Tv = CP.PropsSI("T", "P", ps[i], "Q", 1, "Ethanol") #Saturation temp
        xe = 3.53*(Rad[i]*2)*(1+(x/(3.53*Rad[i])+0.000000001)**(-1.2))**(-1/1.2)
        Gch = ρs[i]*Us[i]
        Tm = 0.5*(Ts[i]+Tv)

        def Ul(Ul):
            const = ((0.0592)/(2*0.023))**5
            a = (Gch*(Ts[i]/Tm)*((Us[i]-Ul)/Us[i]))**2
            b = const*(1/(xe*Rad[i]*2))
            return a - b
        print(Bisect(Ul, 0, Us[i], 10*(-20)))
        #---------------------------------------------------------
        

        ρc = CP.PropsSI("D", "P", ps[i], "T", IV.FuelTankT-15, "Ethanol")
        muc = CP.PropsSI("V", "P", ps[i], "T", IV.FuelTankT-15, "Ethanol")

        def Ul2(Ul):
            return m.sqrt((0.0592*BLCMdot*((Gch/(Tm*Us[i]))**0.8)*((Us[i]-Ul)**1.8))/(muc*ρc))-Ul
        print(Bisect(Ul2, 0, Us[i], 10*(-20)))
        print(Us[i])
        print("---------------------------")

    return x

ceaOut = runCEA()
BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
CombustionGas = ceaOut[0]
#print(emittance(CombustionGas, 3*0.0254))
#print(((calcBLC()[249]*ceaOut[1])+(BLCMdot*(0.8*calcBLC()[249])*0.9)/9.81)/(BLCMdot+ceaOut[1]))
#print(RatL(LT)/12)
#Values = AxialValues(CombustionGas.T, CombustionGas.P, CombustionGas.density, CombustionGas)
#print(Values[1])
#print(CombustionGas.report())
print(calcBLC())
