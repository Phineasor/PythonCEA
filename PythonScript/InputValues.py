# General InputValues
yaml = "PCRL-Mech1.yaml"


# Input numbers in the following section are in inches and degrees(expet the ratios they are ratios)
ExpRatio = 3
ConRatio = 9
Dt = 1
Rcont = 0.5
ThetaCont = 45
Rexp = 0.5
ThetaExp = 15
Lchamber = 3
# End of input section for chamber and nozzle geometry

# Method of radius calculation (linear/parabolic, 1/2)
Method = 1
# End of Method input

# Number of cells in the 1D CEA
CellNum = 2500
# End of cells

# Injector charataristics, Ox holes, fuel holes BLC holes
FuelOrificeDiameter = 0.0550
FuelOrificeNum = 8
FuelOrificeCd = 0.7
FuelAngle = 0  # NOT CURRENTLY USED

OxOrificeDiameter = 0.0591
OxOrificeNum = 8
OxOrificeCd = 0.7
OxAngle = 0  # NOT CURRENTLY USED

BLCOrificeDiameter = 0.0156
BLCOrificeNum = 24
BLCOrificeCd = 0.7
BLCAngle = 0  # NOT CURRENTYL USED this will probably be less helpfull than the other 2 im cookin
# End of Injector pro

# Tank charataristics
FuelTankP = 350  # PSI
FuelTankT = 290  # K
Fuel = "Ethanol"

OxTankP = 350  # PSI
Ox = "O2"
OxTankT = 85  # K
# End Of tank charataristics

# Ambiant Conditions
AmbP = 101325  # In pascales, sorry not PSI here
AmbT = 300  # Ambient temp in K
