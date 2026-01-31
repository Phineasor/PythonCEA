import math as m
import numpy as np
from EngineGeometry import RatL, LT
import InputValues as IV
from Bisect import Bisect

import matplotlib.pyplot as plt


#This fuiction gets the ray data from an arbitrary location across the engine in terms of cell number, and an arbitrary pair of angles measures from pointing at the top of the engine
#the other directionality is irrelivent because the engine is symetric from any plane around the axis
#We will obtain the power values at linspace location across the engine section, however we will quite probably need to polynomialy fit data to get some points inbetween

def getRay(x, theta1, theta2):
    #determines the distance and the radial distance from the z axis at each locatio. Y will be straight twards the wall in a positive direction, X is measured from the injector face in a posive manor following into the end
    AxialDistances = np.linspace(0, LT, IV.CellNum)
    RadiusVal = [0.0]*IV.CellNum
    for i in range(IV.CellNum):
        RadiusVal[i] = RatL(AxialDistances[i])

    #Next we need to obtain the vector normal to the plane that we are working with for the angle of the chamber wall. This is found with the cross product

    dx = 10**(-6) #value for finding the 2d derivative in the xy plane
    
    #will use the slope of the point farther down the engine unless is the last cell, then it will be the one behiend it
    if not(x == (IV.CellNum-1)):
        p1 = np.array([AxialDistances[x], RadiusVal[x], 0])
        p2 = np.array([AxialDistances[x]+dx, RatL(AxialDistances[x]+dx), 0])
        
        Slope = p2-p1
        NormalVector = np.cross(Slope, np.array([0, 0, 1]))
        UNormalV = NormalVector/np.linalg.norm(NormalVector)
    else: 
        p1 = np.array([AxialDistances[x], RadiusVal[x], 0])
        p2 = np.array([AxialDistances[x]-dx, RatL(AxialDistances[x]-dx), 0])
    
        Slope = p1-p2
        NormalVector = np.cross(Slope, np.array([0, 0, 1]))
        UNormalV = NormalVector/np.linalg.norm(NormalVector)
    #p1 and p2 are flipped to produce the negetive NormalVector such that the vectors used to construct it can be used in a change of basis matrix
    #we now have the unit normal vector that is orthogonal to the wall at the locaiton, this is used to cefine the zy plane reziding in this region
    
    #we now need the 3d change of basis matrix so that we can transform these vectors vetween the engine centered frame, and the local wall frame, this is important for doing the rotation matrix on the vector to point in a new direction.
    #conviniently we know all three basis vectors, the normal vector for the x axis, adn the two vectors used to construct it in the array
    USlope = Slope/np.linalg.norm(Slope)
    W = np.array([
        [(UNormalV[0]), (USlope[0]), (0)],
        [(UNormalV[1]), (USlope[1]), (0)],
        [(UNormalV[2]), (USlope[2]), (1)]
    ])
    #now in theory this is the transformation matrix so lets normalize the vector that we need to rotate
    RotateVector = np.linalg.inv(W) @ UNormalV
    
    
    #now we need to apply two rotation matricies to rotate this vector dependent on the theta1 and theta2
    RmatPitch   = np.array([
        [(m.cos(-theta1)), (-m.sin(-theta1)), (0)],
        [(m.sin(-theta1)), (m.cos(-theta1)), (0)],
        [(0), (0), (1)]
    ])
    RmatYaw = np.array([
        [(m.cos(theta2)), (0), (m.sin(theta2))],
        [(0), (1), (0)],
        [(-m.sin(theta2)), (0), (m.cos(theta2))]
    ])
    #the theta1s in RmatYaw are to ensure that for 0-90 theta2 the value remains positive, this is not a requirement but it makes it more convinient
    
    #now we just compute the fully rotated vector
    RotatedVector = RotateVector @ RmatPitch @ RmatYaw #get rotated idiot
    RotatedVectorUnrot = W @ RotatedVector 
    
    ChamberZ = lambda x, z: m.sqrt((RatL(x))**2-(z)**2)
    ChamberZn= lambda x, z: -m.sqrt((RatL(x))**2-(z)**2)
    line = lambda t, point, slope:[(point[0]+t*slope[0]), (point[1]+t*slope[1]), (point[2]+t*slope[2])]
    func1 = lambda t, point, slope: line(t, point, slope)[1]-ChamberZ(line(t, point, slope)[0], line(t, point, slope)[2])
    func2 = lambda t, point, slope: line(t, point, slope)[1]-ChamberZn(line(t, point, slope)[0], line(t, point, slope)[2])

    tval = Bisect(func2, (10**(-6)), 3, (10**(-10)), p1, RotatedVectorUnrot)
    intersect = np.array(line(tval, p1, RotatedVectorUnrot))
   

    
    ray = 0
    return [p1, RotatedVectorUnrot, intersect, (intersect[1]**2+intersect[2]**2)**0.5, tval]



print(getRay(0, (10*(m.pi/180)), (20*(m.pi/180))))
print(getRay(0, (89*(m.pi/180)), (0*(m.pi/180))))
print(getRay(0, (85*(m.pi/180)), (0*(m.pi/180))))
#print(getRay(249, 0, 0))


''''
#RAYMARCH testing 
plt.rcParams['figure.dpi'] = 500
fig = plt.figure(figsize=(10,10))
ax = fig.add_subplot(111,projection='3d')
ax.set_facecolor('black')


num = 100
Edist = np.linspace(0, LT, num)
i = 0
while(i < 360):
    j = 0
    list1 = [0]*num
    list2 = [0]*num
    list3 = [0]*num
    while(j < num):
        list1[j] = Edist[j]
        list2[j] = m.sin((m.pi/180)*i)*RatL(Edist[j])
        list3[j] = m.cos((m.pi/180)*i)*RatL(Edist[j])
        j += 1
    ax.plot(list1,list2,list3, color='white',linestyle='-',linewidth=1) 
    i += 8

    
    






ax.set_xlim3d(-LT/2, LT/2)
ax.set_ylim3d(-LT/2, LT/2)
ax.set_zlim3d(-LT/2, LT/2)
ax.grid(False)
plt.title("3D engine model")
plt.legend(fancybox=False, shadow=True, framealpha=1,fontsize='small',loc='lower left')
plt.show()

'''