# File to generate absorbtion coefficient data from hapi
# fmt: off
"""
@author: phineas
"""
#External modules
import matplotlib.pyplot as plt
import numpy as np
from hapi import *
#from tables import *

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

AbsCoefAray2 = np.load('AbsCoefData2.npy', allow_pickle=True)
print(AbsCoefAray2.dtype)
#print(AbsCoefAray[240])
'''
AbsCoefAray = np.array([None]*IV.CellNum)

moleFrac = cea[0]['H2O', 'CO', 'CO2'].X
for i in range(IV.CellNum):
        nu,coef = absorptionCoefficient_HT(SourceTables=['H2O', 'CO', 'CO2'], WavenumberStep = 0.001, Environment = {'p':((AxVal[1][i])*pa2atm), 'T':AxVal[0][i]}, OmegaWingHW = 100, HITRAN_units = False, Diluent = {'CO2':moleFrac[0], 'CO':moleFrac[1], 'H2O':moleFrac[2]})
        AbsCoefAray[i] = [nu, coef]
np.save('AbsCoefData', AbsCoefAray)
'''

#plt.rcParams['figure.dpi'] = 500
#fig, ax = plt.subplots()
#ax.set_facecolor('white')
'''
i = 0
while i < IV.CellNum:
        print(len(AbsCoefAray[i][0]))
        i += 1




AbsCoefAray2 = np.array([[[0.0]*len(AbsCoefAray[0][0]), [0.0]*len(AbsCoefAray[0][0])]]*250)
count = 0
i = 0
while i < 250:
        j = 0
        while j < 2:
                k = 0
                while k < len(AbsCoefAray[0][0]):
                        AbsCoefAray2[i][j][k] = AbsCoefAray[i][j][k]
                        k += 1
                        count += 1
                        
                        if count % 100000 == 0:
                                print(count)
                j += 1
        i += 1

'''
#np.save('AbsCoefAray2', AbsCoefAray2)