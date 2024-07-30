# CEA calculations for rocket engines.
"""
@author: phineas
"""

from EngineGeometry import RatL, LT
from InputValues import CellNum
import mathfunctions as mf
import matplotlib.pyplot as mp
import cantera as ct
import CoolProp.CoolProp as CP

AxialDistances = mf.linspace(0, LT, CellNum)
RadiusVal = []

i = 0
while i < len(AxialDistances):
    RadiusVal.append(RatL(AxialDistances[i]))
    i += 1

# print(AxialDistances)
# print("LineBreak")
# print(RadiusVal)

# mp.plot(AxialDistances, RadiusVal)
# mp.show()

Pc = 1918481.3852
of = 1.4013

TRef = max(
    CP.PropsSI("T", "P", Pc, "Q", 1, "Ethanol") * 1.01,
    CP.PropsSI("T", "P", Pc, "Q", 1, "O2") * 1.01,
)
print(TRef)

gas1 = ct.Solution("PCRL-Mech1.yaml")

gas1.TPY = TRef, Pc, "O2:1.4013, C2H5OH:1"
gas1.equilibrate("HP")


hCorrectionF = CP.PropsSI("H", "T", 285, "P", Pc, "Ethanol") - CP.PropsSI(
    "H", "T", TRef, "P", Pc, "Ethanol"
)
hCorrectionO = CP.PropsSI("H", "T", 85, "P", Pc, "O2") - CP.PropsSI(
    "H", "T", TRef, "P", Pc, "O2"
)
hCant = gas1.h

h = hCant + (hCorrectionF * (1 / (1 + of))) + (hCorrectionO * (of / (1 + of)))
print(h)
print(gas1.h)


gas1.HP = h, Pc
gas1.equilibrate("HP")


print(gas1())

gamma = gas1.cp / gas1.cv
