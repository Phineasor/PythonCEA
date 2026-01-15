import matplotlib.pyplot as plt

import numpy as np
from pylab import show,plot,subplot,xlim,ylim,title,legend,xlabel,ylabel
from hapi import *


#fetch('CO2',2,1,0,5000)
#fetch('CO2_2',2,2,0,5000)
#fetch('CO2_3',2,3,0,5000)
#nu,coef = absorptionCoefficient_HT(SourceTables='CO2', WavenumberStep = 0.001, Environment = {'p':20.,'T':3500.})
#nu2,coef2 = absorptionCoefficient_HT(SourceTables='CO2_2', WavenumberStep = 0.001, Environment = {'p':20.,'T':3500.})
#nu3,coef3 = absorptionCoefficient_HT(SourceTables='CO2_3', WavenumberStep = 0.001, Environment = {'p':20.,'T':3500.})

#nu2,coef2 = absorptionCoefficient_HT(SourceTables='CO2', Diluent={'air':1.0}, WavenumberStep = 0.0001)


plt.rcParams['figure.dpi'] = 500
fig, ax = plt.subplots()
ax.set_facecolor('white')


i = 1
while i < 7:
    fetch('CO',5,i ,0,25000)
    nu,coef = absorptionCoefficient_HT(SourceTables='CO', WavenumberStep = 0.001, Environment = {'p':20.,'T':3200.})
    ax.plot(nu,coef, color=np.random.rand(3))
    i += 1


#ax.plot(nu,coef, color="red")
#ax.plot(nu2,coef2, color="green")
#ax.plot(nu3,coef3, color="blue")

ax.grid(False)
plt.title("3D engine model")
plt.legend(fancybox=False, shadow=False, framealpha=1,fontsize='small',loc='lower left')
plt.savefig("test9.png")

#getHelp(PROFILE_VOIGT)
