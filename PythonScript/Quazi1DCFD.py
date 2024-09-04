# I think this will end up being a quazi 1D CFD for providing more precice chemical reactions along the nozzle
# from EngineGeometry import RatL, LT, Lt
import InputValues as IV

# from ChamberPressure import ChamberPressure
##import Injector as Inj
# import mathfunctions as mf
import cantera as ct

# import CoolProp.CoolProp as CP
# import math
# from Bisect import Bisect
import numpy

OF = 1.6
TRef = 400
Pc = 100000
CombustionGas = ct.Solution(IV.yaml)
CombustionGas.TPY = TRef, Pc, "O2:" + str(OF) + ", C2H5OH:1"
print(CombustionGas.forward_rates_of_progress.size)
print(CombustionGas.n_reactions)
print(CombustionGas.advance_coverages(0.001))
