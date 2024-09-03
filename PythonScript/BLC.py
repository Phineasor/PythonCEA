# What the fuck.
# No like truly this is a lot
# fmt: off
import InputValues as IV
from cea import runCEA
import math as m
import Injector as Inj
from EngineGeometry import LT, Dc

pa2atm = 9.86923 * 10 ** (-6)

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

def calcBLC():
    #Starting Paramaters
    Tbl = IV.FuelTankT
    x = 0
    dL = LT/IV.CellNum
    L = 0
    D = Dc

    ceaOut = runCEA()
    CombustionGas = ceaOut[1]
    γ = CombustionGas.cp/CombustionGas.cv


    BLCMdot = IV.BLCOrificeNum * Inj.MdotSPIONLY( IV.BLCOrificeCd, IV.BLCOrificeDiameter, IV.Fuel, IV.FuelTankT, ceaOut[0].P, IV.FuelTankP)
    for i in range(IV.CellNum):
        x = i
        return x

CombustionGas = runCEA()[0]
print(emittance(CombustionGas, 3*0.0254))
