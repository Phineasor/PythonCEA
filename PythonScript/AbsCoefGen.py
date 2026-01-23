# File to generate absorbtion coefficient data from hapi
# fmt: off
"""
@author: phineas
"""
#External modules
import matplotlib.pyplot as plt
import numpy as np
from hapi import *
from tables import *

#Internal modules
from cea import *
import InputValues as IV

#Tells Hapi where the data is stored
db_begin('Data')

#conversion factors
pa2atm = 9.86923*(10**(-6))

#Needed engine data to produce abscoeff
cea = runCEA()
AxVal = AxialValues(cea[0].T, cea[0].P, cea[0].density, cea[0])

AbsCoefAray = np.load('AbsCoefData.npy', allow_pickle=True)
print(AbsCoefAray[240])
'''
AbsCoefAray = np.array([None]*IV.CellNum)

moleFrac = cea[0]['H2O', 'CO', 'CO2'].X
for i in range(IV.CellNum):
        nu,coef = absorptionCoefficient_HT(SourceTables=['H2O', 'CO', 'CO2'], WavenumberStep = 0.001, Environment = {'p':((AxVal[1][i])*pa2atm), 'T':AxVal[0][i]}, OmegaWingHW = 100, HITRAN_units = False, Diluent = {'CO2':moleFrac[0], 'CO':moleFrac[1], 'H2O':moleFrac[2]})
        AbsCoefAray[i] = [nu, coef]
np.save('AbsCoefData', AbsCoefAray)
'''