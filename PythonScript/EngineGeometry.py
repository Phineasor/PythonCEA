# Calculates the needed engine geometry paramaters based on the given infromation in the input file
"""
@author: phineas
"""

from InputValues import (
    ExpRatio,
    ConRatio,
    Dt,
    Rcont,
    ThetaCont,
    Rexp,
    ThetaExp,
    Lchamber,
    Method,
)
import mathfunctions as mf
import math

if Method == 1:
    # Calculates needed Engine dimension values from Grissom
    Dc = Dt * math.sqrt(ConRatio)
    De = Dt * math.sqrt(ExpRatio)
    D2 = Dc - 2 * Rcont * (1 - mf.cosd(ThetaCont))
    D3 = Dt + 2 * Rexp * (1 - mf.cosd(ThetaCont))
    D5 = Dt + 2 * Rexp * (1 - mf.cosd(ThetaExp))
    L2 = Lchamber + Rcont * mf.sind(ThetaCont)
    L3 = L2 + ((D2 - D3) / (2 * mf.tand(ThetaCont)))
    Lt = L3 + Rexp * mf.sind(ThetaCont)
    L5 = Lt + Rexp * mf.sind(ThetaExp)
    LT = L5 + ((De - D5) / (2 * mf.tand(ThetaExp)))
    # End of engine dimension value calculation


def RatL(L):
    # Start of contour calculation
    if L < Lchamber:
        return 0.5 * (Dc)
    if L < L2:
        return 0.5 * (Dc - 2 * (Rcont - math.sqrt(Rcont**2 - (L - Lchamber) ** 2)))
    if L < L3:
        return 0.5 * (D2 - 2 * (L - L2) * (mf.tand(ThetaCont)))
    if L < L5:
        return 0.5 * (Dt + 2 * (Rexp - math.sqrt(Rexp**2 - (L - Lt) ** 2)))
    else:
        return 0.5 * (D5 + 2 * (L - L5) * (mf.tand(ThetaExp)))
