#! /user/bin/python
# -*- coding:UTF-8 -*-
# filename：3DSpheres.py
from abaqus import *
from abaqusConstants import *
from caeModules import *
import os
import numpy as np
import math
import numpy as np
from visualization import *
from odbAccess import *
import math
session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry= COORDINATE)
Mdb()
myModel=mdb.Model(name='#name your model')
Modelname = myModel.name
# rectangular region
ConcreteLength=50.0  
ConcreteWidth=50.0    
ConcreteHeight=50.0    

# basic information
ConcreteVolume=ConcreteLength*ConcreteWidth*ConcreteHeight                 
AggRatio=0.3                                                   
TargetVolume =ConcreteVolume*AggRatio;                          
TotalAggVolume= 0.0                                            
AggGap = 1.01
BaseGapMax = 0.99
BaseGapMin = 0.01


AggLimite = 5000                                               
Dmax_1 = 4.75                                               
Dmin_1 = 2.36                                               
IterLimite = 10000                                         


def first_agg(point):
    x1 = point[0]
    y1 = point[1]
    z1 = point[2]
    r1 = point[3]
    sign = True
    x_z = x1 + r1
    x_f = x1 - r1
    y_z = y1 + r1
    y_f = y1 - r1
    z_z = z1 + r1
    z_f = z1 - r1
    if x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength or y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin *  ConcreteWidth or z_z > BaseGapMax *ConcreteHeight or z_f < BaseGapMin* ConcreteHeight :
        sign = False
    return sign

# Interaction check
def interact_judgement(points, point):
    x1 = point[0]
    y1 = point[1]
    z1 = point[2]
    r1 = point[3]
    sign = True
    x_z = x1 + r1
    x_f = x1 - r1
    y_z = y1 + r1
    y_f = y1 - r1
    z_z = z1 + r1
    z_f = z1 - r1
    if x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength or y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin *  ConcreteWidth or z_z > BaseGapMax *ConcreteHeight or z_f < BaseGapMin* ConcreteHeight :
        sign = False
    for ii in points:
        x2 = ii[0]
        y2 = ii[1]
        z2 = ii[2]
        r2 = ii[3]
        # distance calculation   
        distance = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if distance < AggGap *(r1+r2) :
            sign = False
            break
    return sign

 
k=0
Aggpoint = []             
Aggpoints = []            
AggData = []             

i=0
for Aggnum in range(AggLimite):
    radius= pow( ( np.random.random((1,1))[0][0]*(pow(Dmax_1,0.5)-pow(Dmin_1,0.5)) ) + pow(Dmin_1,0.5) , 2) /2.0 
    for iter in range (IterLimite):
        if iter < IterLimite - 1 :
            x1=np.random.uniform(0+radius,ConcreteLength-radius)      
            y1=np.random.uniform(0+radius,ConcreteWidth-radius)      
            z1=np.random.uniform(0+radius,ConcreteHeight-radius)      
            
            Aggpoint = (x1,y1,z1,radius)  #x y z r
            
            if len(Aggpoints) == 0:
                Aggpoints.append(Aggpoint)
                AggData.append([i+1,x1,y1,z1,radius])
                GenerateVolume=4*3.14*radius*radius*radius/3
                TotalAggVolume = TotalAggVolume + GenerateVolume
                Process = TotalAggVolume*100/TargetVolume
                Friction = TotalAggVolume*100/ConcreteVolume
                i = i+1
                print('Success!')
                print('Current Volume = {:.3f}, Process = {:.3f}%, Friction = {:.3f} %Current Step = {}'.format(TotalAggVolume,Process,Friction,Aggnum))
                print('Current Iter time = {:.3f}'.format(iter))
                break
            elif interact_judgement(Aggpoints,Aggpoint):
                Aggpoints.append(Aggpoint)
                AggData.append([i+1,x1,y1,z1,radius])
                GenerateVolume=4*3.14*radius*radius*radius/3
                TotalAggVolume = TotalAggVolume + GenerateVolume
                Process = TotalAggVolume*100/TargetVolume
                Friction = TotalAggVolume*100/ConcreteVolume
                i = i+1
                print('Success!')
                print('Current Volume = {:.3f}, Process = {:.3f}%, Friction = {:.3f}% Current Step = {}'.format(TotalAggVolume,Process,Friction,Aggnum))
                print('Current Iter time = {:.3f}'.format(iter))
                break
        else:
            print('Fail!')
            print('Generate Failure')
            break
    if (TotalAggVolume - TargetVolume) >= 0.001:
        break

b=AggData

# shape information
NumLayers = 5
NumMiddleLayer = 8
NumNode1st = 1
NumNode2nd = 8
NumNode3rd = 8
NumNode4th = 8
NumNode5th = 1
AngleGap = 180.0 / (NumLayers-1)                          
NumTotNodes = NumNode1st+NumNode2nd+NumNode3rd+NumNode4th+NumNode5th

NumAgg = len(AggData)                                            
# aggregate information
RandomNodesPor = np.zeros((NumAgg,NumTotNodes,4))       
RandomNodesCar = np.zeros((NumAgg,NumTotNodes,4))       
#Shell to Solid
NumFaces = NumMiddleLayer*2 + (NumLayers-2-1) * NumMiddleLayer*2
Faces = np.zeros((NumFaces,13))                        
FacesITZ = np.zeros((NumFaces,13)) 

AggRadius = np.zeros((NumAgg,1))                          
AggCentroid = np.zeros((NumAgg,4))                      



RandomNodesITZCar = np.zeros((NumAgg,NumTotNodes,4))       
AggITZData = []                                            

for Aggnum in range(len(AggData)):
    FaceFlage = 0                    
    FaceFlageITZ = 0
    # Create aggregate part
    SizeLength=200
    mdb.models[Modelname].ConstrainedSketch(name='__profile__',sheetSize=1)
    mdb.models[Modelname].sketches['__profile__'].rectangle(point1=(-SizeLength/2, SizeLength/2), 
    point2=(SizeLength/2, -SizeLength/2))
    mdb.models[Modelname].Part(dimensionality=THREE_D, name='PolyAgg-'+str(Aggnum), type=DEFORMABLE_BODY)
    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].BaseSolidExtrude(depth=SizeLength/2, sketch=
    mdb.models[Modelname].sketches['__profile__'])
    del mdb.models[Modelname].sketches['__profile__']
    del mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].features['Solid extrude-1']
    # Create ITZregion part
    mdb.models[Modelname].ConstrainedSketch(name='__profile__',sheetSize=1)
    mdb.models[Modelname].sketches['__profile__'].rectangle(point1=(-SizeLength/2, SizeLength/2), 
    point2=(SizeLength/2, -SizeLength/2))
    mdb.models[Modelname].Part(dimensionality=THREE_D, name='PolyAggITZregion-'+str(Aggnum), type=DEFORMABLE_BODY)
    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].BaseSolidExtrude(depth=SizeLength/2, sketch=
    mdb.models[Modelname].sketches['__profile__'])
    del mdb.models[Modelname].sketches['__profile__']
    del mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].features['Solid extrude-1']
    RadiusVio = 0.15         # Node fluctuation coefficient along the radius (relative to radius [0,1])
    AngleVio = 9           # Node fluctuation range along the angle (theta & phi)
    RadiusVioAdjust = 1     # Determine whether to fluctuate along the radius direction [0 or 1]
    AngleVioAdjust = 1      # Determine whether to fluctuate along the angle direction [0 or 1]
    radius = AggData[Aggnum][4]  
    for layer in range(NumLayers):                     
        theta = layer * AngleGap/180*math.pi       # Use theta to divide each layer, for 5-layer aggregate, each layer is 45°
        ############################################## First layer (top layer) node [1 node] #############################################
        if abs(layer-0) <= 1e-10:                           # Convergence check
            for Node in range(NumNode1st):
                # Coeff = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                phi = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][layer][0] = layer+1    # Node number of top layer 1: 1
                RandomNodesPor[Aggnum][layer][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           # Control radius of top layer 1
                RandomNodesPor[Aggnum][layer][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta of top layer 1
                RandomNodesPor[Aggnum][layer][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta of top layer 1
        ##############################################'' Last layer node [1 node] #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            for Node in range(NumNode5th):
                # Coeff = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                phi = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][0] = (layer-1)*NumNode2nd+NumNode1st+1                
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])    
        else:
            for Node in range(NumMiddleLayer):
                phi = Node * AngleGap/180 * math.pi         
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][0] = Node+(layer-1)*NumNode2nd+NumNode1st+1                
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])          
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])   
    '''*********************************************************** Convert spherical coordinates por to Cartesian Car coordinates *********************************************************'''
    # x = r*sin(theta)*cos(phi);      y = r*sin(theta)*sin(phi);        z = r*cos(theta)
    for Node in range(NumTotNodes):
        RandomNodesCar[Aggnum][Node][0] = RandomNodesPor[Aggnum][Node][0]
        RandomNodesCar[Aggnum][Node][1] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.cos(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][2] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.sin(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][3] = RandomNodesPor[Aggnum][Node][1]*math.cos(RandomNodesPor[Aggnum][Node][2])
    '''**************************************************************************************************************'''


    a1 = [0,0,0]
    a2 = [0,0,0]
    a3 = [0,0,0]
    a4 = [0,0,0]
    a5 = [0,0,0]
    a6 = [0,0,0]
    b1 = (0,0,0)
    b2 = (0,0,0)
    b3 = (0,0,0)
    c1 = (0,0,0)
    c2 = (0,0,0)
    c3 = (0,0,0)

    for layer in range(NumLayers):

        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][layer][1]
            a1[1] = RandomNodesCar[Aggnum][layer][2]
            a1[2] = RandomNodesCar[Aggnum][layer][3]

            b1 = tuple(a1)

            for Node in range(NumNode2nd):

                a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                b2 = tuple(a2)

                a3[0] = RandomNodesCar[Aggnum][Node+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node+NumNode1st+1][3]
                b3 = tuple(a3)

                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)

                    a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)

                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))

        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]

            b1 = tuple(a1)

            for Node in range(NumNode4th):

                a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)

                a3[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)

                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)

                    a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)

                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))

        elif abs(layer - 3) >= 1e-10:
            for Node in range(NumMiddleLayer):
                a1[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a1[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a1[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                b1 = tuple(a1)          
                a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)          
                a3[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)          
                a4[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a4[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a4[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                c1 = tuple(a4)          
                a5[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
                a5[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
                a5[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
                c2 = tuple(a5)          
                a6[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a6[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a6[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                c3 = tuple(a6)          
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      
                    a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      
                
                    a5[0] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          
                    a6[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))


    for layer in range(NumLayers):
        ############################################# Close the top layer face #############################################
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][layer][1]
            a1[1] = RandomNodesCar[Aggnum][layer][2]
            a1[2] = RandomNodesCar[Aggnum][layer][3]
            # Convert list coordinates to tuple [first vertex, connect with all nodes below]
            b1 = tuple(a1)
            # Second layer 8 points loop, connect with vertex
            for Node in range(NumNode2nd):
                # Connect point 2
                a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]      # x
                a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]      # y
                a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]      # z
                b2 = tuple(a2)
                # Connect point 3
                a3[0] = RandomNodesCar[Aggnum][Node+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node+NumNode1st+1][3]
                b3 = tuple(a3)
                # After one round of connection, the last connected face is 192, need to update the second point to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 to the second point, numbered 2
                    a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS default method to define the surface:
                Faces[FaceFlage][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FacesITZ[FaceFlageITZ][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FaceFlage += 1
                FaceFlageITZ += 1
                
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]
            b1 = tuple(a1)
            for Node in range(NumNode4th):
                a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)
                a3[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                Faces[FaceFlage][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FacesITZ[FaceFlageITZ][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FaceFlage += 1
                FaceFlageITZ += 1
                
        else:
            for Node in range(NumMiddleLayer):
                a1[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a1[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a1[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                b1 = tuple(a1)          
                a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)          
                a3[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)          
                a4[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a4[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a4[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                c1 = tuple(a4)          
                a5[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
                a5[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
                a5[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
                c2 = tuple(a5)          
                a6[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a6[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a6[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                c3 = tuple(a6)          
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      
                    a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      
                
                    a5[0] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          
                    a6[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          
                Faces[FaceFlage][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FacesITZ[FaceFlageITZ][:] = [layer+1, Node+1, Node+2, 1, Node+1, Node+2, layer+2, Node+1, Node+2, 1, Node+1, Node+2, 1]
                FaceFlage += 1
                FaceFlageITZ += 1

    AggRadius[Aggnum][0] = AggData[Aggnum][4]       
    AggCentroid[Aggnum][0] = AggData[Aggnum][1]
    AggCentroid[Aggnum][1] = AggData[Aggnum][2]
    AggCentroid[Aggnum][2] = AggData[Aggnum][3]


for layer in range(NumLayers):
    for Node in range(NumTotNodes):
        AggData[Aggnum][layer][Node] = 1

# Generating wireframe for polyhedral aggregates

a = [0, 0, 0]
b = (0, 0, 0)
for layer in range(NumLayers):
    for Node in range(NumMiddleLayer):
        a[0] = RandomNodesCar[Aggnum][Node + layer * NumMiddleLayer + NumNode1st][1]
        a[1] = RandomNodesCar[Aggnum][Node + layer * NumMiddleLayer + NumNode1st][2]
        a[2] = RandomNodesCar[Aggnum][Node + layer * NumMiddleLayer + NumNode1st][3]
        b = tuple(a)
        mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
            ON, points=((b, b), (b, b), (b, b)))

# Code block for defining surface faces, lines, and nodes
FaceFlage = 0
FaceFlageITZ = 0

for layer in range(NumLayers):
    for Node in range(NumTotNodes):
        if layer == 0:
            Faces[FaceFlage][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FacesITZ[FaceFlageITZ][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FaceFlage += 1
            FaceFlageITZ += 1
        elif layer == (NumLayers - 1):
            Faces[FaceFlage][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FacesITZ[FaceFlageITZ][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FaceFlage += 1
            FaceFlageITZ += 1
        else:
            Faces[FaceFlage][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FacesITZ[FaceFlageITZ][:] = [layer + 1, Node + 1, Node + 2, 1, Node + 1, Node + 2, layer + 2, Node + 1, Node + 2, 1, Node + 1, Node + 2, 1]
            FaceFlage += 1
            FaceFlageITZ += 1

# Generate final polyhedral shape using wireframe definition
mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].Surface(mnemonic='PolyAggSurface')

for layer in range(NumLayers):
    mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].Set(faces=
        mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].faces.findAt(((1, 2, 3),)), name=
        'Set-' + str(layer + 1))

for layer in range(NumLayers):
    mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].Set(nodes=
        mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].nodes.findAt(((1, 2, 3),)), name=
        'Set-' + str(layer + 1))

# Store polyhedral shape data
mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].Set(name='PolyAggData', faces=
    mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].faces)

for layer in range(NumLayers):
    mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].Set(name='PolyAggCentroid', nodes=
        mdb.models[Modelname].parts['PolyAgg-' + str(Aggnum)].nodes.findAt(((1, 2, 3),)))

