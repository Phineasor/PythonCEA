# Function that calculates chamber pressure from choked flow at the known throat area and input gas properties
"""
@author: phineas
"""

from InputValues import Dt
import math


# The formula used here was derived from Sutton on page 59 eq 3-24
inch2meter = 0.0254


def ChamberPressure(Tc, mdot, gamma, R):
    exp = -0.5 * ((gamma + 1) / (gamma - 1))
    u = (gamma * R * Tc) ** 0.5
    At = math.pi * (Dt * 0.0254 / 2) ** 2
    MdotTerm = mdot / (At * gamma)

    return MdotTerm * u * (2 / (gamma + 1)) ** exp
