! /user/bin/python
# -*- coding:UTF-8 -*-
# filename：3DSpheres.py
from abaqus import *
from abaqusConstants import *
from caeModules import *
import os
import numpy as np
import math
from visualization import *
from odbAccess import *
import math
session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry= COORDINATE)
Mdb()
myModel=mdb.Model(name='# your model name')
Modelname = myModel.name
# rectangular region
ConcreteLength=50.0  # mm, Length of rectangular concrete
ConcreteWidth=50.0    # mm, Width of rectangular concrete
ConcreteHeight=50.0    # mm, Height of rectangular concrete

# Modeling based on the percentage of aggregate, not the number of aggregates
ConcreteVolume=ConcreteLength*ConcreteWidth*ConcreteHeight                 # Volume of concrete
AggRatio=0.35                                                   # Aggregate ratio
TargetVolume =ConcreteVolume*AggRatio;                         # Target volume of aggregate
TotalAggVolume= 0.0                                            # Cumulative aggregate volume
AggGap = 1.01
BaseGapMax = 0.99
BaseGapMin = 0.01

# AGG iteration method
AggLimite = 20000                                               # Maximum iterations for Agg
# Adopt Fuller grading and compare with the CCR article
Dmax_1 = 4.75                                               
Dmin_1 = 2.36                                                                                      
IterLimite = 10000                                       

# Check if the first sphere is inside:
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

# Function to check if the newly generated sphere intersects with previously generated spheres
def interact_judgement(points, point):
    # points: existing spheres, point: newly generated sphere, judge if they intersect
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
        # distance calculate   Judgment standard: The distance between the two centers is greater than the sum of their radii
        distance = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if distance < AggGap *(r1+r2) :
            sign = False
            break
    return sign

# generate spheres
# SphereNum=30        # Number of spherical aggregates
k=0
Aggpoint = []             # Accumulated data of generated spherical aggregates
Aggpoints = []            # Intermediate variable
AggData = []             # Data of spherical aggregates


i=0
for Aggnum in range(AggLimite):
    radius= pow( ( np.random.random((1,1))[0][0]*(pow(Dmax_1,0.5)-pow(Dmin_1,0.5)) ) + pow(Dmin_1,0.5) , 2) /2.0  # Aggregate radius
    '''*********************************Aggregate placement judgment*****************************************************************************'''
    for iter in range (IterLimite):
        if iter < IterLimite - 1 :
            # Generate random centroid coordinates xyz
            x1=np.random.uniform(0+radius,ConcreteLength-radius)     # x-coordinate of the centroid, mm   
            y1=np.random.uniform(0+radius,ConcreteWidth-radius)      # y-coordinate of the centroid, mm  
            z1=np.random.uniform(0+radius,ConcreteHeight-radius)     # z-coordinate of the centroid, mm    
            # Store the random xyz and the longest radius of the generated aggregate:
            Aggpoint = (x1,y1,z1,radius)  # x y z r(radius)
            # Judge if the generated aggregate intersects with others:
            if len(Aggpoints) == 0:
                Aggpoints.append(Aggpoint)
                # Store the generated aggregate into the total data
                AggData.append([i+1,x1,y1,z1,radius])
                # Calculate the volume of this aggregate
                GenerateVolume=4*3.14*radius*radius*radius/3
                # Record the total generated volume:
                TotalAggVolume = TotalAggVolume + GenerateVolume
                Process = TotalAggVolume*100/TargetVolume
                Friction = TotalAggVolume*100/ConcreteVolume
                # Update i
                i = i+1
                print('Success!')
                print('Current Volume = {:.3f}, Process = {:.3f}%, Friction = {:.3f} %Current Step = {}'.format(TotalAggVolume,Process,Friction,Aggnum))
                print('Current Iter time = {:.3f}'.format(iter))
                break
            elif interact_judgement(Aggpoints,Aggpoint):
                # Store the qualified aggregate coordinates and radius into the total
                Aggpoints.append(Aggpoint)
                # Store the generated aggregate into the total data
                AggData.append([i+1,x1,y1,z1,radius])
                # Calculate the volume of this aggregate
                GenerateVolume=4*3.14*radius*radius*radius/3
                # Record the total generated volume:
                TotalAggVolume = TotalAggVolume + GenerateVolume
                Process = TotalAggVolume*100/TargetVolume
                Friction = TotalAggVolume*100/ConcreteVolume
                # Update i
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
'''Generate aggregate********************************************************************************'''
# create in Abaqus
b=AggData


## Define a single aggregate:
NumLayers = 5
NumMiddleLayer = 8
NumNode1st = 1
NumNode2nd = 8
NumNode3rd = 8
NumNode4th = 8
NumNode5th = 1
AngleGap = 180.0 / (NumLayers-1)                          # Angle size of equal division
NumTotNodes = NumNode1st+NumNode2nd+NumNode3rd+NumNode4th+NumNode5th

NumAgg = len(AggData)                                             # Number of polygonal aggregates
# Information of each aggregate in the polygonal aggregate set:
RandomNodesPor = np.zeros((NumAgg,NumTotNodes,4))       # Spherical coordinates (point number, radius, theta, phi)
RandomNodesCar = np.zeros((NumAgg,NumTotNodes,4))       # Cartesian coordinates (point number, x, y, z)
# Surface definition when converting SHELL to solid:
NumFaces = NumMiddleLayer*2 + (NumLayers-2-1) * NumMiddleLayer*2
Faces = np.zeros((NumFaces,13))                         # Saved area and, record order as (face number, first point(3), second point(3), third point(3), center point(3))
FacesITZ = np.zeros((NumFaces,13)) 

AggRadius = np.zeros((NumAgg,1))                          # Aggregate diameter set
AggCentroid = np.zeros((NumAgg,4))                      # Aggregate number and centroid coordinates (xyz) for each generated aggregate

# ITZ region modeling:
RandomNodesITZCar = np.zeros((NumAgg,NumTotNodes,4))       # ITZ region Cartesian coordinates (point number, x, y, z)
AggITZData = []                                            # ITZ data

for Aggnum in range(len(AggData)):
    FaceFlage = 0                    # Record the face number generated, starting from 0, a total of 48 faces
    FaceFlageITZ = 0
    # Create Aggregate part
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
    RadiusVio = 0.15         # Node fluctuation coefficient along the radius (radial) [relative to radius[0,1]]
    ## Angle fluctuations can affect convergence, it is recommended to keep around 12, or adjust the fluctuation range, otherwise generation may fail
    AngleVio = 9           # Node fluctuation range along the angle (theta & phi)
    RadiusVioAdjust = 1     # Determine whether to fluctuate along the radius direction [0 or 1]
    AngleVioAdjust = 1      # Determine whether to fluctuate along the angle direction [0 or 1]
    radius = AggData[Aggnum][4]  # Radius is already generated in advance
    for layer in range(NumLayers):                      # Distribute nodes randomly for each layer
        theta = layer * AngleGap/180*math.pi       # Use theta to divide each layer, for 5 layers of aggregates, the angle for each layer is 45°
        ############################################## First layer (top layer) node [1 node] #############################################
        if abs(layer-0) <= 1e-10:                           # Convergence check
            for Node in range(NumNode1st):
                #Coe = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                fai = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][layer][0] = layer+1    # Node number for the top layer 1st point: 1
                RandomNodesPor[Aggnum][layer][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           # Control radius for the top layer 1st point
                RandomNodesPor[Aggnum][layer][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta for the top layer 1st point
                RandomNodesPor[Aggnum][layer][3] = fai + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta for the top layer 1st point
        ##############################################''Last layer node [1 node] #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            for Node in range(NumNode5th):
                # Coe = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                fai = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][0] = (layer-1)*NumNode2nd+NumNode1st+1                # Node number for the bottom layer 26th point: 1
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           # Control radius for the bottom layer 26th point
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta for the bottom layer 26th point
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3] = fai + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])    # Theta for the bottom layer 26th point
        ##################################################### Middle three layers [8 nodes] ###############################################
        else:
            for Node in range(NumMiddleLayer):
                # Coe = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                fai = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][0] = Node+(layer-1)*NumNode2nd+NumNode1st+1                # Node number for the bottom layer 26th point: 1
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])          # Control radius for the bottom layer 26th point
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Theta for the bottom layer 26th point
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3] = fai + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])   # Theta for the bottom layer 26th point
    '''***********************************************************Convert spherical coordinates to Cartesian coordinates*********************************************************'''
    #x = r*sin(theta)*cos(fai);      y = r*sin(theta)*sin(fai);        z = r*cos(theta)
    for Node in range(NumTotNodes):
        RandomNodesCar[Aggnum][Node][0] = RandomNodesPor[Aggnum][Node][0]
        RandomNodesCar[Aggnum][Node][1] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.cos(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][2] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.sin(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][3] = RandomNodesPor[Aggnum][Node][1]*math.cos(RandomNodesPor[Aggnum][Node][2])
    '''**************************************************************************************************************'''
    '''***********************************************************Connect 26 points into lines*********************************************************'''
    # Define an empty list a to assign points to the list, define an empty tuple b to convert the list into a tuple for line connection, c is used for middle layer connection
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
    '''***********************************************************Points -> Lines: Connect the previous layer points with the next layer: 123, 145, 156, 167, 178, 189, 192*********************************************************'''
    for layer in range(NumLayers):
        ############################################# Top layer (1 point) connected to the next layer #############################################
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][layer][1]
            a1[1] = RandomNodesCar[Aggnum][layer][2]
            a1[2] = RandomNodesCar[Aggnum][layer][3]
            # Convert the list coordinates to a tuple [the first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # The second layer has 8 points, connected to the top vertex
            for Node in range(NumNode2nd):
                # Connect point 2
                a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                b2 = tuple(a2)
                # Connect point 3
                a3[0] = RandomNodesCar[Aggnum][Node+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node+NumNode1st+1][3]
                b3 = tuple(a3)
                # After one round of connection, the last connection surface is 192, so the second point needs to be updated to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connection point b3 to the second node, numbered 2
                    a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect the three points into three lines to form a closed surface:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        ############################################# The bottom point [1 point] is connected to the previous layer (penultimate layer [8 points]): #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]
            # Convert the list coordinates to a tuple [the first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # The penultimate layer has 8 points, connected to the bottom point:
            for Node in range(NumNode4th):
                # Connect point 2
                a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)
                # Connect point 3
                a3[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)
                # After one round of connection, the last connection surface is 192, so the second point needs to be updated to b3:
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connection point b3 to the second node, numbered 2
                    a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect the three points into three lines to form a closed surface:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        ############################################# Middle layers (2,3) -> (1 2) [excluding 3] 8 nodes correspond to 8 node connections:########################################################
                a1[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a1[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a1[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                b1 = tuple(a1)          # Node 2
                a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)          # Node 10
                a3[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)          # Node 11
                # Generate the triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the single upper +1 node
                a4[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a4[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a4[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                c1 = tuple(a4)          # Node 2
                a5[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
                a5[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
                a5[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
                c2 = tuple(a5)          # Node 3
                a6[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a6[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a6[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                c3 = tuple(a6)          # Node 11
                # After one round of connection, update the two triangles ** Node 9 does not need to be changed, so b1 c1 does not need to be updated:
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    ## Connect 9 - 17 -10 (at this time NODE = 7)
                    a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      # Node 17
                    a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      # Node 10
                    ## Connect 9 - 10 -2 (at this time NODE = 7)                  
                    a5[0] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          # Node 2
                    a6[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          # Node 10
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
for layer in range(NumLayers):
    ############################################# Close the top layer #############################################
    if abs(layer - 0) <= 1e-10:
        a1[0] = RandomNodesCar[Aggnum][layer][1]
        a1[1] = RandomNodesCar[Aggnum][layer][2]
        a1[2] = RandomNodesCar[Aggnum][layer][3]
        # Convert the list coordinates to a tuple [the first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # The second layer has 8 points, connected to the top vertex
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
            # After one round of connection, the last connection surface is 192, so the second point needs to be updated to b3:
            if abs(Node - (NumNode2nd-1)) <= 1e-10:
                a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the edges of the triangle's three sides through findat, and then generate the surface:
            mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (8 in total), record the center point coordinates of each surface (extract the three vertices of the triangle and then calculate the center point coordinates):
            Faces[FaceFlage][0] = FaceFlage + 1                 # Record the first surface, starting from number 1
            # First point coordinates
            Faces[FaceFlage][1] = a1[0]                         # x-coordinate of the first point
            Faces[FaceFlage][2] = a1[1]                         # y-coordinate of the first point
            Faces[FaceFlage][3] = a1[2]                         # z-coordinate of the first point
            # Second point coordinates
            Faces[FaceFlage][4] = a2[0]                         # x-coordinate of the second point
            Faces[FaceFlage][5] = a2[1]                         # y-coordinate of the second point
            Faces[FaceFlage][6] = a2[2]                         # z-coordinate of the second point
            # Third point coordinates
            Faces[FaceFlage][7] = a3[0]                         # x-coordinate of the third point
            Faces[FaceFlage][8] = a3[1]                         # y-coordinate of the third point
            Faces[FaceFlage][9] = a3[2]                         # z-coordinate of the third point
            # Center point coordinates (need to be extracted to help locate the surface)
            Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3        # x-coordinate of the center point
            Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3        # y-coordinate of the center point
            Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3        # z-coordinate of the center point
            # After finding each surface, increment the surface number counter:
            FaceFlage += 1                                      # Update the surface number counter
    ############################################# Close the bottom layer: #############################################
    elif abs(layer - (NumLayers-1)) <= 1e-10:
        a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
        a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
        a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]
        # Convert the list coordinates to a tuple [the first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # The penultimate layer has 8 points, connected to the bottom point:
        for Node in range(NumNode4th):
            # Connect point 2
            a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)
            # Connect point 3
            a3[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)
            # After one round of connection, the last connection surface is 192, so the second point needs to be updated to b3:
            if abs(Node - (NumNode4th-1)) <= 1e-10:
                a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the edges of the triangle's three sides through findat, and then generate the surface:
            mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (8 in total), record the center point coordinates of each surface (extract the three vertices of the triangle and then calculate the center point coordinates):
            Faces[FaceFlage][0] = FaceFlage + 1                 # Record the first surface, starting from number 1
            # First point coordinates
            Faces[FaceFlage][1] = a1[0]                         # x-coordinate of the first point
            Faces[FaceFlage][2] = a1[1]                         # y-coordinate of the first point
            Faces[FaceFlage][3] = a1[2]                         # z-coordinate of the first point
            # Second point coordinates
            Faces[FaceFlage][4] = a2[0]                         # x-coordinate of the second point
            Faces[FaceFlage][5] = a2[1]                         # y-coordinate of the second point
            Faces[FaceFlage][6] = a2[2]                         # z-coordinate of the second point
            # Third point coordinates
            Faces[FaceFlage][7] = a3[0]                         # x-coordinate of the third point
            Faces[FaceFlage][8] = a3[1]                         # y-coordinate of the third point
            Faces[FaceFlage][9] = a3[2]                         # z-coordinate of the third point
            # Center point coordinates (need to be extracted to help locate the surface)
            Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3        # x-coordinate of the center point
            Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3        # y-coordinate of the center point
            Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3        # z-coordinate of the center point
            # After finding each surface, increment the surface number counter:
            FaceFlage += 1                                      # Update the surface number counter
    ############################################# Close the middle layers: ########################################################
    elif abs(layer - 3) >= 1e-10:
        ## The middle layers consist of two triangles combined, b and c each control one triangle
        for Node in range(NumMiddleLayer):
            # Generate the triangle controlled by b, starting from 2-10-11
            a1[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a1[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a1[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            b1 = tuple(a1)          # Node 2
            a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)          # Node 10
            a3[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)          # Node 11
            # Generate the triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the single upper +1 node
            a4[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a4[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a4[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            c1 = tuple(a4)          # Node 2
            a5[0] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
            a5[1] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
            a5[2] = RandomNodesCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
            c2 = tuple(a5)          # Node 3
            a6[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a6[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a6[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            c3 = tuple(a6)          # Node 11
            # After one round of connection, update the two triangle surfaces ** Node 9 does not need to be changed, so b1 and c1 do not need to be updated:
            if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                ## Connect 9 - 17 -10 (at this time NODE = 7)
                a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)      # Node 17
                a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                b3 = tuple(a3)      # Node 10
                ## Connect 9 - 10 -2 (at this time NODE = 7)                  
                a5[0] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                a5[1] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                a5[2] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                c2 = tuple(a5)          # Node 2
                a6[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a6[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a6[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                c3 = tuple(a6)          # Node 10
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the edges of the triangle's three sides through findat, and then generate the surface:
            mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
                # For the middle layers, not all connections need to be considered when connecting, but the surfaces are generated as triangles, so the points for b and c need to be generated simultaneously
            mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c2[0]+c3[0])/2, (c2[1]+c3[1])/2, (c2[2]+c3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c3[0]+c1[0])/2, (c3[1]+c1[1])/2, (c3[2]+c1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (16 in total), record the center point coordinates of each surface (extract the three vertices of the triangle and then calculate the center point coordinates):
            Faces[FaceFlage][0] = FaceFlage + 1                 # Record the first surface, starting from number 1
            # Surface controlled by the triangle b
            # First point coordinates
            Faces[FaceFlage][1] = a1[0]                         # x-coordinate of the first point
            Faces[FaceFlage][2] = a1[1]                         # y-coordinate of the first point
            Faces[FaceFlage][3] = a1[2]                         # z-coordinate of the first point
            # Second point coordinates
            Faces[FaceFlage][4] = a2[0]                         # x-coordinate of the second point
            Faces[FaceFlage][5] = a2[1]                         # y-coordinate of the second point
            Faces[FaceFlage][6] = a2[2]                         # z-coordinate of the second point
            # Third point coordinates
            Faces[FaceFlage][7] = a3[0]                         # x-coordinate of the third point
            Faces[FaceFlage][8] = a3[1]                         # y-coordinate of the third point
            Faces[FaceFlage][9] = a3[2]                         # z-coordinate of the third point
            # Center point coordinates (need to be extracted to help locate the surface)
            Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3        # x-coordinate of the center point
            Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3        # y-coordinate of the center point
            Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3        # z-coordinate of the center point
            # After finding each surface, increment the surface number counter:
            FaceFlage += 1                                      # Update the surface number counter
            # Surface controlled by the triangle c
            # First point coordinates
            Faces[FaceFlage][1] = a4[0]                         # x-coordinate of the first point
            Faces[FaceFlage][2] = a4[1]                         # y-coordinate of the first point
            Faces[FaceFlage][3] = a4[2]                         # z-coordinate of the first point
            # Second point coordinates
            Faces[FaceFlage][4] = a5[0]                         # x-coordinate of the second point
            Faces[FaceFlage][5] = a5[1]                         # y-coordinate of the second point
            Faces[FaceFlage][6] = a5[2]                         # z-coordinate of the second point
            # Third point coordinates
            Faces[FaceFlage][7] = a6[0]                         # x-coordinate of the third point
            Faces[FaceFlage][8] = a6[1]                         # y-coordinate of the third point
            Faces[FaceFlage][9] = a6[2]                         # z-coordinate of the third point
            # Center point coordinates (need to be extracted to help locate the surface)
            Faces[FaceFlage][10] = (a4[0]+a5[0]+a6[0])/3        # x-coordinate of the center point
            Faces[FaceFlage][11] = (a4[1]+a5[1]+a6[1])/3        # y-coordinate of the center point
            Faces[FaceFlage][12] = (a4[2]+a5[2]+a6[2])/3        # z-coordinate of the center point
            # After finding each surface, increment the surface number counter:
            FaceFlage += 1                                      # Update the surface number counter
'''**************************************************************************************************************'''
'''******************* Surface -> Solid: Record all center points of surfaces generated for a single aggregate into the facelist collection and find them in ABAQUS *********************************************************'''
# At present, all surfaces have been generated, and the center point coordinates of the surfaces have been saved and recorded in Faces[10](x coordinate); Faces[11](y coordinate); Faces[12](z coordinate)
Facelist = []       # Define an empty set to store all surfaces to be converted into solids
# Store all surfaces into facelist through loop:
for i in range(int(Faces.shape[0])):
    Facelist.append(mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].faces.findAt((Faces[i][10], Faces[i][11],Faces[i][12]), ) )
# Use the cell command in ABAQUS to convert the surface into a solid:
mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].AddCells(faceList=Facelist)
# Place each part into the assembly to become an independent instance
mdb.models[Modelname].rootAssembly.Instance(dependent=ON, name='PolyAgg-'+str(Aggnum), 
    part=mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)])

'''********************************* Prepare coordinate parameters for aggregate placement judgment criteria *****************************************************************************'''
## Use ABAQUS's getMassProperties to read the centroid of each generated aggregate:
AggCentroid[Aggnum][0] = Aggnum+1                                                                # Record aggregate number
AggCentroid[Aggnum][1]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][0]  # x-coordinate of the centroid
AggCentroid[Aggnum][2]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][1]  # y-coordinate of the centroid
AggCentroid[Aggnum][3]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][2]  # z-coordinate of the centroid
'''***************************************ITZ region****************************************************'''
### ITZ region
DatumsPoints = []
for ITZregion in range(NumTotNodes):
    a1 = [0,0,0]
    a2 = [0,0,0]
    b1 = (0,0,0)
    b2 = (0,0,0)
    a1[0] = RandomNodesCar[Aggnum][ITZregion][1]
    a1[1] = RandomNodesCar[Aggnum][ITZregion][2]
    a1[2] = RandomNodesCar[Aggnum][ITZregion][3]
    b1 = tuple(a1)
    a2[0] = AggCentroid[Aggnum][1]
    a2[1] = AggCentroid[Aggnum][2]
    a2[2] = AggCentroid[Aggnum][3]
    b2 = tuple(a2)
    # Use the SEPARATE merge type, otherwise the wire generation will fail
    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=SEPARATE, meshable=
        ON, points=(b2,b1))
    # Create offset point (at position 0.1)
    '''Generate ITZ region (ITZ thickness varies with the radius of the aggregate, but for efficiency of placement, choose to generate inward)'''
    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].DatumPointByEdgeParam(edge=
        mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2)),
        parameter=0.95, isDependent=False)
    '''According to the radius of each aggregate and the fixed ITZ thickness, calculate the parameter'''
    # Delete the generated line segment for later calculation
    del mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].features['Wire-'+str(33)]
# Traverse through the loop to store all the ID numbers into DatumsPoints    
for datum in range(NumTotNodes):   
    # Get all DatumPoints ID numbers 
    DatumsPoints.append(mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].datums[84+datum*2])
# Use the pointOn command to save the left and right point coordinates to the ITZ node region:
for Node in range(NumTotNodes):
    RandomNodesITZCar[Aggnum][Node][0] = RandomNodesPor[Aggnum][Node][0]
    RandomNodesITZCar[Aggnum][Node][1] = DatumsPoints[Node].pointOn[0]
    RandomNodesITZCar[Aggnum][Node][2] = DatumsPoints[Node].pointOn[1]
    RandomNodesITZCar[Aggnum][Node][3] = DatumsPoints[Node].pointOn[2]

'''*********************************************************** Connect 26 points into lines *********************************************************'''
# Define an empty list a, assign points to the list, define an empty tuple b, convert the list into a tuple for easy line connection, c is used for middle layer connection
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
layer = 0
'''***********************************************************Point -> Line: Connect points between the upper layer and the next layer: 123, 145, 156, 167, 178, 189, 192*********************************************************'''
for layer in range(NumLayers):
    ############################################# Top layer (1 point) connected to the next layer #############################################
    if abs(layer - 0) <= 1e-10:
        a1[0] = RandomNodesITZCar[Aggnum][layer][1]
        a1[1] = RandomNodesITZCar[Aggnum][layer][2]
        a1[2] = RandomNodesITZCar[Aggnum][layer][3]
        # Convert list coordinates to tuple [first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # Loop through 8 points in the second layer, connected to the top point
        for Node in range(NumNode2nd):
            # Connect point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]
            a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]
            a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]
            b2 = tuple(a2)
            # Connect point 3
            a3[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][3]
            b3 = tuple(a3)
            # After a round of connection, the last connection is 192, update the second point to b3:
            if abs(Node - (NumNode2nd-1)) <= 1e-10:
                a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesITZCar[Aggnum][NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                ON, points=((b1,b2), (b2, b3), (b3, b1)))
    ############################################# Bottom layer (1 point) connected to the previous layer (second-to-last layer with 8 points) #############################################
    elif abs(layer - (NumLayers-1)) <= 1e-10:
        a1[0] = RandomNodesITZCar[Aggnum][NumTotNodes-1][1]
        a1[1] = RandomNodesITZCar[Aggnum][NumTotNodes-1][2]
        a1[2] = RandomNodesITZCar[Aggnum][NumTotNodes-1][3]
        # Convert list coordinates to tuple [first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # Connect the 8 points in the second-to-last layer to the bottom point:
        for Node in range(NumNode4th):
            # Connect point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)
            # Connect point 3
            a3[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)
            # After a round of connection, the last connection is 192, update the second point to b3:
            if abs(Node - (NumNode4th-1)) <= 1e-10:
                a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                ON, points=((b1,b2), (b2, b3), (b3, b1)))
    ############################################# Middle layers (2,3) -> (1, 2) [excluding 3], 8 points connected to 8 points: ########################################################
    #####********************************* The middle layer must establish a relationship with the layer!!! Otherwise, surface generation will have problems *********************************#####
    elif abs(layer - 3) >= 1e-10:
        ## The middle layer is composed of two triangles, b and c each control one triangle
        for Node in range(NumMiddleLayer):
            # Generate triangle controlled by b, starting from 2-10-11
            a1[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a1[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a1[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            b1 = tuple(a1)          # Point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)          # Point 10
            a3[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)          # Point 11
            # Generate triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the node +1 on the upper layer
            a4[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a4[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a4[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            c1 = tuple(a4)          # Point 2
            a5[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
            a5[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
            a5[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
            c2 = tuple(a5)          # Point 3
            a6[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a6[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a6[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            c3 = tuple(a6)          # Point 11
            # After a round of connection, update the two triangle surfaces ** No need to update point 9, so b1 c1 do not need to be updated:
            if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                ## Connect 9 - 17 -10 (Node = 7)
                a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)      # Point 17
                a3[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                b3 = tuple(a3)      # Point 10
                ## Connect 9 - 10 -2 (Node = 7)                  
                a5[0] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                a5[1] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                a5[2] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                c2 = tuple(a5)          # Point 2
                a6[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a6[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a6[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                c3 = tuple(a6)          # Point 10
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                ON, points=((b1,b2), (b2, b3), (b3, b1)))
'''*************************Line -> Surface: Merge the closed line segments into surfaces using coveredge and findat methods [need to define a new variable to store all generated surfaces for later assignment of solid attributes] *********************************************************'''
for layer in range(NumLayers):
    ############################################# Close the top surface #############################################
    if abs(layer - 0) <= 1e-10:
        a1[0] = RandomNodesITZCar[Aggnum][layer][1]
        a1[1] = RandomNodesITZCar[Aggnum][layer][2]
        a1[2] = RandomNodesITZCar[Aggnum][layer][3]
        # Convert list coordinates to tuple [first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # Loop through 8 points in the second layer, connected to the top point
        for Node in range(NumNode2nd):
            # Connect point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]      # x
            a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]      # y
            a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]      # z
            b2 = tuple(a2)
            # Connect point 3
            a3[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st+1][3]
            b3 = tuple(a3)
            # After a round of connection, the last connection is 192, update the second point to b3:
            if abs(Node - (NumNode2nd-1)) <= 1e-10:
                a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesITZCar[Aggnum][NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the middle points of the three triangle sides using findat, then generate the surface:
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (8 surfaces), record the center coordinates of each surface (extract the three vertices of the triangle and then calculate the midpoint coordinates):
            FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start numbering from 1 for the first surface
            # First point coordinates
            FacesITZ[FaceFlageITZ][1] = a1[0]                         # First point x coordinate
            FacesITZ[FaceFlageITZ][2] = a1[1]                         # First point y coordinate
            FacesITZ[FaceFlageITZ][3] = a1[2]                         # First point z coordinate
            # Second point coordinates
            FacesITZ[FaceFlageITZ][4] = a2[0]                         # Second point x coordinate
            FacesITZ[FaceFlageITZ][5] = a2[1]                         # Second point y coordinate
            FacesITZ[FaceFlageITZ][6] = a2[2]                         # Second point z coordinate
            # Third point coordinates
            FacesITZ[FaceFlageITZ][7] = a3[0]                         # Third point x coordinate
            FacesITZ[FaceFlageITZ][8] = a3[1]                         # Third point y coordinate
            FacesITZ[FaceFlageITZ][9] = a3[2]                         # Third point z coordinate
            # Center coordinates (need to be extracted for surface positioning)
            FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Center point x coordinate
            FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Center point y coordinate
            FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Center point z coordinate
            # After finding each surface, increment the surface counter:
            FaceFlageITZ += 1                                      # Update the surface counter
    ############################################# Close the bottom surface: #############################################
    elif abs(layer - (NumLayers-1)) <= 1e-10:
        a1[0] = RandomNodesITZCar[Aggnum][NumTotNodes-1][1]
        a1[1] = RandomNodesITZCar[Aggnum][NumTotNodes-1][2]
        a1[2] = RandomNodesITZCar[Aggnum][NumTotNodes-1][3]
        # Convert list coordinates to tuple [first vertex, connected to all nodes below]
        b1 = tuple(a1)
        # Connect the 8 points in the second-to-last layer to the bottom point:
        for Node in range(NumNode4th):
            # Connect point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)
            # Connect point 3
            a3[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)
            # After a round of connection, the last connection is 192, update the second point to b3:
            if abs(Node - (NumNode4th-1)) <= 1e-10:
                a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)
                # Update the last connection point b3 to the second node, numbered 2
                a3[0] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                b3 = tuple(a3)
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the middle points of the three triangle sides using findat, then generate the surface:
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (8 surfaces), record the center coordinates of each surface (extract the three vertices of the triangle and then calculate the midpoint coordinates):
            FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start numbering from 1 for the first surface
            # First point coordinates
            FacesITZ[FaceFlageITZ][1] = a1[0]                         # First point x coordinate
            FacesITZ[FaceFlageITZ][2] = a1[1]                         # First point y coordinate
            FacesITZ[FaceFlageITZ][3] = a1[2]                         # First point z coordinate
            # Second point coordinates
            FacesITZ[FaceFlageITZ][4] = a2[0]                         # Second point x coordinate
            FacesITZ[FaceFlageITZ][5] = a2[1]                         # Second point y coordinate
            FacesITZ[FaceFlageITZ][6] = a2[2]                         # Second point z coordinate
            # Third point coordinates
            FacesITZ[FaceFlageITZ][7] = a3[0]                         # Third point x coordinate
            FacesITZ[FaceFlageITZ][8] = a3[1]                         # Third point y coordinate
            FacesITZ[FaceFlageITZ][9] = a3[2]                         # Third point z coordinate
            # Center coordinates (need to be extracted for surface positioning)
            FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Center point x coordinate
            FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Center point y coordinate
            FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Center point z coordinate
            # After finding each surface, increment the surface counter:
            FaceFlageITZ += 1                                      # Update the surface counter
    ############################################# Close the middle surfaces: ########################################################
    elif abs(layer - 3) >= 1e-10:
        ## The middle layer is composed of two triangles, b and c each control one triangle
        for Node in range(NumMiddleLayer):
            # Generate triangle controlled by b, starting from 2-10-11
            a1[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a1[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a1[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            b1 = tuple(a1)          # Point 2
            a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
            a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
            a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
            b2 = tuple(a2)          # Point 10
            a3[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a3[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a3[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            b3 = tuple(a3)          # Point 11
            # Generate triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the node +1 on the upper layer
            a4[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
            a4[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
            a4[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
            c1 = tuple(a4)          # Point 2
            a5[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
            a5[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
            a5[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
            c2 = tuple(a5)          # Point 3
            a6[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
            a6[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
            a6[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
            c3 = tuple(a6)          # Point 11
            # After a round of connection, update the two triangle surfaces ** No need to update point 9, so b1 c1 do not need to be updated:
            if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                ## Connect 9 - 17 -10 (Node = 7)
                a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)      # Point 17
                a3[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a3[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a3[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                b3 = tuple(a3)      # Point 10
                ## Connect 9 - 10 -2 (Node = 7)                  
                a5[0] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                a5[1] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                a5[2] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                c2 = tuple(a5)          # Point 2
                a6[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                a6[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                a6[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                c3 = tuple(a6)          # Point 10
            ############ Use ABAQUS's coveredge method to close the line segments into surfaces; find the middle points of the three triangle sides using findat, then generate the surface:
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                tryAnalytical=True)
                # For the middle layer surfaces, you can skip some connections, but surfaces are generated based on triangles, so both b's points and c's points need to be generated
            mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2), ),   # First point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c2[0]+c3[0])/2, (c2[1]+c3[1])/2, (c2[2]+c3[2])/2), ),   # Second point
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c3[0]+c1[0])/2, (c3[1]+c1[1])/2, (c3[2]+c1[2])/2), )),  # Third point
                tryAnalytical=True)
            ############ After generating the surface (16 surfaces), record the center coordinates of each surface (extract the three vertices of the triangle and then calculate the midpoint coordinates):
            FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start numbering from 1 for the first surface
            # Triangle surface controlled by b
            # First point coordinates
            FacesITZ[FaceFlageITZ][1] = a1[0]                         # First point x coordinate
            FacesITZ[FaceFlageITZ][2] = a1[1]                         # First point y coordinate
            FacesITZ[FaceFlageITZ][3] = a1[2]                         # First point z coordinate
            # Second point coordinates
            FacesITZ[FaceFlageITZ][4] = a2[0]                         # Second point x coordinate
            FacesITZ[FaceFlageITZ][5] = a2[1]                         # Second point y coordinate
            FacesITZ[FaceFlageITZ][6] = a2[2]                         # Second point z coordinate
            # Third point coordinates
            FacesITZ[FaceFlageITZ][7] = a3[0]                         # Third point x coordinate
            FacesITZ[FaceFlageITZ][8] = a3[1]                         # Third point y coordinate
            FacesITZ[FaceFlageITZ][9] = a3[2]                         # Third point z coordinate
            # Center coordinates (need to be extracted for surface positioning)
            FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Center point x coordinate
            FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Center point y coordinate
            FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Center point z coordinate
            # After finding each surface, increment the surface counter:
            FaceFlageITZ += 1                                      # Update the surface counter
            # Triangle surface controlled by c
            # First point coordinates
            FacesITZ[FaceFlageITZ][1] = a4[0]                         # First point x coordinate
            FacesITZ[FaceFlageITZ][2] = a4[1]                         # First point y coordinate
            FacesITZ[FaceFlageITZ][3] = a4[2]                         # First point z coordinate
            # Second point coordinates
            FacesITZ[FaceFlageITZ][4] = a5[0]                         # Second point x coordinate
            FacesITZ[FaceFlageITZ][5] = a5[1]                         # Second point y coordinate
            FacesITZ[FaceFlageITZ][6] = a5[2]                         # Second point z coordinate
            # Third point coordinates
            FacesITZ[FaceFlageITZ][7] = a6[0]                         # Third point x coordinate
            FacesITZ[FaceFlageITZ][8] = a6[1]                         # Third point y coordinate
            FacesITZ[FaceFlageITZ][9] = a6[2]                         # Third point z coordinate
            # Center coordinates (need to be extracted for surface positioning)
            FacesITZ[FaceFlageITZ][10] = (a4[0]+a5[0]+a6[0])/3                        # Center point x coordinate
            FacesITZ[FaceFlageITZ][11] = (a4[1]+a5[1]+a6[1])/3                        # Center point y coordinate
            FacesITZ[FaceFlageITZ][12] = (a4[2]+a5[2]+a6[2])/3                        # Center point z coordinate
            # After finding each surface, increment the surface counter:
            FaceFlageITZ += 1                                      # Update the surface counter
'''*******************Surface -> Solid: Record the center coordinates of all surfaces generated for a single aggregate into the facelist set and find them in abaqus*********************************************************'''
FacelistITZ = []       # Define an empty set to store all surfaces to be converted into solids
# Store all surfaces into facelist through a loop:
for k in range(int(FacesITZ.shape[0])):
    FacelistITZ.append( mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].faces.findAt((FacesITZ[k][10], FacesITZ[k][11],FacesITZ[k][12]), ) )
# Convert the surfaces into solids in ABAQUS using the cell command:
mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].AddCells(faceList=FacelistITZ)
# Place each part into the assembly to become an independent instance
mdb.models[Modelname].rootAssembly.Instance(dependent=ON, name='PolyAggITZregion-'+str(Aggnum), 
    part=mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)])

print('Current Aggregate = {},Total Aggregates = {}'.format(Aggnum,len(AggData)))

'''**********Translate multiple polyhedral aggregates in the ASSEMBLY module****************************************************************************************************'''
b=AggData
for Aggnum in range(len(AggData)):
    # Translate the aggregate region
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAgg-'+str(b[Aggnum][0]-1),),vector=(b[Aggnum][1],0.0, 0.0)) 
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAgg-'+str(b[Aggnum][0]-1),),vector=(0.0, b[Aggnum][2],0.0 ))       
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAgg-'+str(b[Aggnum][0]-1),),vector=(0.0, 0.0, b[Aggnum][3]))   
    # Translate the ITZ region
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAggITZregion-'+str(b[Aggnum][0]-1),),vector=(b[Aggnum][1],0.0, 0.0)) 
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAggITZregion-'+str(b[Aggnum][0]-1),),vector=(0.0, b[Aggnum][2],0.0 ))       
    mdb.models[Modelname].rootAssembly.translate(instanceList=('PolyAggITZregion-'+str(b[Aggnum][0]-1),),vector=(0.0, 0.0, b[Aggnum][3]))   

'''***********************Generate the outer contour of concrete**********************************************************************'''
# Commands to generate the concrete contour
xmin=0
xmax=ConcreteLength+0.2
ymin=0
ymax=ConcreteWidth+0.2
zmin=0
zmax=ConcreteHeight+0.2
zlength=abs(zmax-zmin)
ConcreteSketch = myModel.ConstrainedSketch(name='concretecube',sheetSize=200)
ConcreteSketch.rectangle(point1=(xmin, ymin), point2=(xmax, ymax))
ConcretePart = myModel.Part(dimensionality=THREE_D, name='concretecube', type=DEFORMABLE_BODY)
myPart = ConcretePart.BaseSolidExtrude(depth=zlength, sketch=myModel.sketches['concretecube'])
del myModel.sketches['concretecube']
myModel.rootAssembly.Instance(name='concretecube', part=ConcretePart, dependent=ON)

'''*************Create unified sets for the generated ExtendAgg and Aggregate, merge them, and delete parts to reduce memory usage************'''
## Aggregate operations:
# Create a set for all aggregates:
Aggregate=mdb.models[Modelname].rootAssembly.instances['PolyAgg-'+str(0)].cells[0:0]   
for i in range(len(AggData)): 
    Aggregate=Aggregate+mdb.models[Modelname].rootAssembly.instances['PolyAgg-'+str(AggData[i][0]-1)].cells
mdb.models[Modelname].rootAssembly.Set(cells=Aggregate, name='ExtendAggregates')
# Merge all AGGREGATES into one part:
AggregateInstance = []
for i in range(len(AggData)):
    AggregateInstance.append(mdb.models[Modelname].rootAssembly.instances['PolyAgg-' + str(AggData[i][0]-1)])
AggregateInstance = tuple(AggregateInstance)
mdb.models[Modelname].rootAssembly.InstanceFromBooleanMerge(name='ExtendAggregates', instances=AggregateInstance, keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
## ITZ operations
# Create a set for all ExtendAggregates:
ExtendAggregate=mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-'+str(0)].cells[0:0]   
for i in range(len(AggData)): 
    ExtendAggregate=ExtendAggregate+mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-'+str(AggData[i][0]-1)].cells
mdb.models[Modelname].rootAssembly.Set(cells=ExtendAggregate, name='InnerAggregates')
AggregateInstance = []
# Merge all ExtendAGGREGATES into one part:
ExtentAggregateInstance = []
for i in range(len(AggData)):
    ExtentAggregateInstance.append(mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-' + str(AggData[i][0]-1)])
ExtentAggregateInstance = tuple(ExtentAggregateInstance)
mdb.models[Modelname].rootAssembly.InstanceFromBooleanMerge(name='InnerAggregates', instances=ExtentAggregateInstance, keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
## Delete all individual parts of Aggregate and ExtendAggregate to free memory:
for i in range(len(AggData)):
    del mdb.models[Modelname].parts['PolyAgg-' + str(AggData[i][0]-1)]
    del mdb.models[Modelname].parts['PolyAggITZregion-' + str(AggData[i][0]-1)]

'''************************************Perform Boolean subtraction operation on the generated ExtendAggregate to obtain the outer contour of the ITZ region************************'''
# Perform Boolean operation to obtain the subtraction ITZ region
mdb.models[Modelname].rootAssembly.InstanceFromBooleanCut(name='ITZRegion', 
    instanceToBeCut=mdb.models[Modelname].rootAssembly.instances['ExtendAggregates-1'], 
    cuttingInstances=(mdb.models[Modelname].rootAssembly.instances['InnerAggregates-1'], ), 
    originalInstances=SUPPRESS)
# Restore the aggregate region
mdb.models[Modelname].rootAssembly.features['InnerAggregates-1'].resume()
