import matplotlib.pyplot as plt

import numpy as np
from pylab import show,plot,subplot,xlim,ylim,title,legend,xlabel,ylabel
from hapi import *
from tables import *


fetch('CO2',2,1,0,8000)
fetch('CO',5,1,0,8000)
fetch('H2O',1,1,0,8000)
nu,coef = absorptionCoefficient_Voigt(SourceTables=['CO2', 'CO', 'H2O'], WavenumberStep = 0.001, Environment = {'p':20.,'T':3200.}, OmegaWingHW = 100, HITRAN_units = True, Diluent = {'CO2':0.1559, 'CO':0.31797, 'H2O':0.52607})
#nu2,coef2 = absorptionCoefficient_HT(SourceTables='CO2_2', WavenumberStep = 0.001, Environment = {'p':20.,'T':3500.})
#nu3,coef3 = absorptionCoefficient_HT(SourceTables='CO2_3', WavenumberStep = 0.001, Environment = {'p':20.,'T':3500.})

#nu2,coef2 = absorptionCoefficient_HT(SourceTables='CO2', Diluent={'air':1.0}, WavenumberStep = 0.0001)


plt.rcParams['figure.dpi'] = 1000
fig, ax = plt.subplots()
ax.set_facecolor('white')

'''
i = 1
while i < 13:
    fetch('CO2',2,i ,0,250000)
    nu,coef = absorptionCoefficient_HT(SourceTables='CO2', WavenumberStep = 0.001, Environment = {'p':20.,'T':3200.}, OmegaWingHW = 100, HITRAN_units = False, Components = [(2, i)], GammaL = 'gamma_self')
    ax.plot(nu,coef, color=np.random.rand(3))
    i += 1
'''

ax.plot(nu,coef, color="blue")
#ax.plot(nu2,coef2, color="green")
#ax.plot(nu3,coef3, color="blue")

ax.grid(False)
plt.title("3D engine model")
plt.legend(fancybox=False, shadow=False, framealpha=1,fontsize='small',loc='lower left')
plt.savefig("test18.png")
print(len(nu))
#getHelp(PROFILE_VOIGT)
