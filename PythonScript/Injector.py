# Right now this is incredibally simplue, just uses the SPI model, wont work for all props work for now though
"""
@author: phineas
"""

import CoolProp as CP
import math


def MdotSPIONLY(Cd, D, prop, TProp, PChamber, PTank):
    rho0 = CP.PropsSI("D", "T", TProp, "Q", 0, prop)
    return Cd * 0.25 * math.pi * (D**2) * (2 * rho0 * (PTank - PChamber)) ** 0.5
