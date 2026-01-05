import math as m
import mathfunctions as mf
import numpy as np
from EngineGeometry import RatL, LT
import InputValues as IV



#This fuiction gets the ray data from an arbitrary location across the engine in terms of cell number, and an arbitrary pair of angles measures from pointing at the top of the engine
#the other directionality is irrelivent because the engine is symetric from any plane around the axis
#We will obtain the power values at linspace location across the engine section, however we will quite probably need to polynomialy fit data to get some points inbetween

def getRay(x, theta1, theta2):
    #determines the distance and the radial distance from the z axis at each locatio. Y will be straight twards the wall in a positive direction, X is measured from the injector face in a posive manor following into the end
    AxialDistances = mf.linspace(0, LT, IV.CellNum)
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
    
    
    
    #The rest of this computation can be compleated in the 2d plane with a surface normal vector that is the result fo the cross product of the rotated vector in 3d space, and the vector pointing up and dows [1, 0, 0]
    #first we need to get what will become the y bassis vector for a new coordinate space
    xbase = np.array([1, 0, 0])
    ybase = np.array([0, RotatedVectorUnrot[1], RotatedVectorUnrot[2]])/np.linalg.norm(np.array([0, RotatedVectorUnrot[1], RotatedVectorUnrot[2]]))
    zbase = np.cross(xbase, ybase)/np.linalg.norm(np.cross(xbase, ybase))
    #we nneed full xyz for continiuty although the z will be dropped for the intersection compute
    
    #Transformation matrix to go from the vectors view to this new one, it should be purely in xy now
    W = np.array([
        [(xbase[0]), (ybase[0]), (zbase[0])],
        [(xbase[1]), (ybase[1]), (zbase[1])],
        [(xbase[2]), (ybase[2]), (zbase[2])]
    ])
    xyPvec = np.linalg.inv(W) @ RotatedVectorUnrot
    
    #now the data in is the new coordinate systems such that it aligns nicely with the RatL function
    xyPoint = np.array([p1[0], p1[1]])
    xySlope = np.array([xyPvec[0], xyPvec[1]])
    
    
    
    
    
    ray = 0
    return [xyPoint, xySlope]



print(getRay(0, (10*(m.pi/180)), (20*(m.pi/180))))
#print(getRay(249, 0, 0))