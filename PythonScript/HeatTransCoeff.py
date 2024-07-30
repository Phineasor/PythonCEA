"""
@author: phineas
"""

import csv
import math

in2m = 0.0254


def sind(theta):
    return math.sin((math.pi / 180) * theta)


def cosd(theta):
    return math.cos((math.pi / 180) * theta)


def tand(theta):
    return math.tan((math.pi / 180) * theta)


theta1 = 45
thetaCon = theta1
theta2 = 15
thetaExp = theta2
expRatio = 3
contRatio = 9
Dt = 1 * in2m

rCon = r1 = 0.5 * in2m
rExp = r2 = 0.5 * in2m
Lc = 3 * in2m
L1 = Lc

Dc = Dt * math.sqrt(contRatio)
De = Dt * math.sqrt(expRatio)
D2 = Dc - (2 * rCon * (1 - cosd(thetaCon)))
D3 = Dt + (2 * rExp * (1 - cosd(thetaCon)))
D5 = Dt + (2 * rExp * (1 - cosd(thetaExp)))

L2 = Lc + (rCon * sind(thetaCon))
L3 = L2 + ((De - D5) / (2 * tand(thetaExp)))
Lt = L3 + (rExp * sind(thetaCon))
L5 = Lt + (rExp * sind(thetaExp))

print(cosd(90))


def getD(L):
    if L < L1:
        return Dc
    if L < L2:
        return Dc - 2 * (r1 - sind((r1**2) - ((L - L1) ** 2)))
    if L < L3:
        return D2 - 2 * (L - L2) * tand(theta1)
    if L < L5:
        return Dt + 2 * (r2 - math.sqrt((r2**2)) - ((L - Lt) ** 2))
    else:
        return D5 + 2 * (L - L5) + tand(theta2)


rho = []
arclength = []
U = []
D = []
Dd = []
h = []
hCOSTANT = 0.023 * ((((6.2182 ** (-10)) ** 0.2) * 2305.3429) / (0.50919**0.6))
hCOSTANT = 11.4569

with open("test4.csv", mode="r") as file:
    csv_reader = csv.reader(file, delimiter=",")
    line_count = 0
    for row in csv_reader:
        if line_count == 0:
            #            print(f'Column names are {", ".join(row)}')
            line_count += 1
        else:
            if (line_count >= 255) & (line_count < 504):
                #                print(f'\t{row[4]} is arc, and {row[24]} is rho')
                arclength.append(row[4])
                U.append(row[17])
                Dd.append(row[23])
                rho.append(row[24])
                line_count += 1
            else:
                if line_count >= 255:
                    Dd.append(row[23])
                line_count += 1
    print(f"Processed {line_count} lines.")
    print(len(arclength))
    print(len(rho))
    print(len(U))
    #    print(hCOSTANT)
    print(len(Dd))
# i=0
# while i < len(Dd):
#    print(Dd[i])
#    i+=1


i = 0
while i < len(arclength):
    D.append(getD(float(arclength[i])))
    #    print(D[i])
    i += 1
print("break")
i = 0
while i < len(arclength):
    #    print(arclength[i])
    i += 1
print("break")
i = 0
while i < len(D):
    h.append(
        (
            hCOSTANT
            * ((float(rho[i]) * float(U[i])) ** 0.8)
            / ((float(Dd[i * 2]) * 2) ** 0.2)
        )
    )
    #    print(Dd[i*2])
    print(h[i])
    i += 1
