#! /user/bin/python
# -*- coding:UTF-8 -*-
# filename：3DSpheres.py
from abaqus import *
from abaqusConstants import *
from caeModules import *
import numpy as np
from visualization import *
from odbAccess import *
from abaqusConstants import *
from caeModules import *
import os
from visualization import *
from odbAccess import *
import math
from ssl import SSLSocket
from abaqusConstants import *
from caeModules import *
import random
import matplotlib.pyplot as plt

session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
Mdb()
myModel=mdb.Model(name='Biochar Concrete')
Modelname = myModel.name
# rectangular region
ConcreteLength=50.0  # mm, length of rectangular concrete
ConcreteWidth=50.0    # mm, width of rectangular concrete
ConcreteHeight=50.0    # mm, height of rectangular concrete


# Modeling based on the percentage of aggregate, not by number
ConcreteVolume=ConcreteLength*ConcreteWidth*ConcreteHeight                 # Concrete volume
AggRatio=0.21                                                   # Aggregate ratio
TargetVolume =ConcreteVolume*AggRatio;                         # Target aggregate volume
TotalAggVolume= 0.0                                            # Accumulated aggregate volume
AggGap = 1.01
BaseGapMax = 0.99
BaseGapMin = 0.01

# AGG iteration method
AggLimite = 10000                                               # Maximum number of AGG iterations
# Adopt Fuller gradation and compare with CCR article
Dmax_1 = 4.75                                               # Maximum aggregate size for the first gradation
Dmin_1 = 2.36                                               # Minimum aggregate size for the first gradation
IterLimite = 10000                                         # Maximum iteration limit


# Determine if the first sphere is inside:
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
    #   points: existing spheres, point: newly generated sphere, judgment between the two
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
        # distance calculate   Judgment criteria: the distance between the centers of the two spheres is greater than the sum of their radii
        distance = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if distance < AggGap *(r1+r2) :
            sign = False
            break
    return sign

# generate spheres
# SphereNum=30        # Number of spherical aggregates
k=0
Aggpoint = []             # Accumulated spherical aggregate data
Aggpoints = []            # Intermediate variable
AggData = []             # Spherical aggregate data

'''Generate spherical aggregates********************************************'''
i=0
for Aggnum in range(AggLimite):
    radius= pow( ( np.random.random((1,1))[0][0]*(pow(Dmax_1,0.5)-pow(Dmin_1,0.5)) ) + pow(Dmin_1,0.5) , 2) /2.0  # Aggregate radius
    '''*********************************Aggregate placement judgment*****************************************************************************'''
    for iter in range (IterLimite):
        if iter < IterLimite - 1 :
            # Generate random centroid coordinates xyz
            x1=np.random.uniform(0+radius,ConcreteLength-radius)     # x-coordinate of the centroid, mm   
            y1=np.random.uniform(0+radius,ConcreteWidth-radius)      # y-coordinate of the centroid, mm  
            z1=np.random.uniform(0+radius,ConcreteHeight-radius)     # z-coordinate of the centroid ,mm    
            # Store the random xyz and the maximum radius of the generated aggregate:
            Aggpoint = (x1,y1,z1,radius)  #x y z r(radius)
            # Sequentially judge the intersection of the generated aggregate:
            if len(Aggpoints) == 0:
                Aggpoints.append(Aggpoint)
                # Store the generated aggregate in the total aggregate data
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
                # Store the qualified aggregate point coordinates and radius in the total aggregate accumulation
                Aggpoints.append(Aggpoint)
                # Store the generated aggregate in the total aggregate data
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




# Biochar Generation
session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
# Define Biochar Library and Biochar model
mdb.Model(name='BiocharLibrary', modelType=STANDARD_EXPLICIT)
mdb.Model(name='BiocharConcrete', modelType=STANDARD_EXPLICIT)
myModel_Lib=mdb.Model(name='BiocharLibrary')
myModel_BC=mdb.Model(name='BiocharConcrete')
Modelname_BC = myModel_BC.name
Modelname_Lib = myModel_Lib.name
# Remove the mask
session.journalOptions.setValues(replayGeometry=INDEX,recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE,recoverGeometry=COORDINATE)
'''*************************************Intersection judgment equation******************************'''
def interact_judgement(points, point):
    #   points: existing aggregates, point: newly generated aggregates, judgment between the two
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
    if (x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength 
    or y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin *  ConcreteWidth 
    or z_z > BaseGapMax *ConcreteHeight or z_f < BaseGapMin* ConcreteHeight) :
        sign = False
    for ii in points:
        x2 = ii[0]
        y2 = ii[1]
        z2 = ii[2]
        r2 = ii[3]
        # distance calculate   Judgment criteria: the distance between the two polyhedra is greater than the sum of their radii
        distance = math.sqrt((x1-x2)**2+(y1-y2)**2+(z1-z2)**2)
        if distance < BioGap *(r1+r2) :
            sign = False
            break
    return sign
'''************************************************Biochar Library geometric definition******************************'''
## Number of items in Biochar Library:
Biochar_Num = 6
## Biochar Library geometric collection definition
BiocharCentroid = np.zeros((Biochar_Num,4))         # Record Biochar number and center coordinates xyz
MaximumDistance_Biochar = np.zeros((Biochar_Num,1)) # Record the maximum distance of Biochar (based on EXCEL)
Length_Biochar = np.zeros((Biochar_Num,1))          # Record the stretch length of Biochar (based on EXCEL)
Volume_Biochar = np.zeros((Biochar_Num,1))          # Record Biochar volume
LDration = np.zeros((Biochar_Num,1))                # Record the aspect ratio of Biochar
## Import geometry from Library (Biochar starts numbering from 1, so str needs to add 1)
for i in range(Biochar_Num):
    iges = mdb.openIges(
        'F:\Papers\Topic 1 Biochar as Fine Aggregate in mortar scale\Biochar Concrete Code\BiocharCode\BiocharGeom/Biochar-'+str(i+1)+'O.sat', msbo=False, 
        trimCurve=DEFAULT, scaleFromFile=OFF)
    mdb.models['BiocharLibrary'].PartFromGeometryFile(name='Biochar-'+str(i+1), 
        geometryFile=iges, combine=False, stitchTolerance=1.0, 
        dimensionality=THREE_D, type=DEFORMABLE_BODY, convertToAnalytical=1, 
        stitchEdges=1)
'''
Take five as example, more geometry can be found in the Biochar Library file.
'''
# Maximum distance
MaximumDistance_Biochar[0][0] = 3.8
MaximumDistance_Biochar[1][0] = 1.7
MaximumDistance_Biochar[2][0] = 2.9
MaximumDistance_Biochar[3][0] = 3.2
MaximumDistance_Biochar[4][0] = 2.8
MaximumDistance_Biochar[5][0] = 3.5
# Stretch length
Length_Biochar[0][0] = 3.8
Length_Biochar[1][0] = 1.7
Length_Biochar[2][0] = 2.9
Length_Biochar[3][0] = 3.2
Length_Biochar[4][0] = 2.8
Length_Biochar[5][0] = 3.5
# Volume
Volume_Biochar[0][0] = 2*3.1416*MaximumDistance_Biochar[0][0]*Length_Biochar[0][0]*Length_Biochar[0][0]/8
Volume_Biochar[1][0] = 2*3.1416*MaximumDistance_Biochar[1][0]*Length_Biochar[1][0]*Length_Biochar[1][0]/8
Volume_Biochar[2][0] = 2*3.1416*MaximumDistance_Biochar[2][0]*Length_Biochar[2][0]*Length_Biochar[2][0]/8
Volume_Biochar[3][0] = 2*3.1416*MaximumDistance_Biochar[3][0]*Length_Biochar[3][0]*Length_Biochar[3][0]/8
Volume_Biochar[4][0] = 2*3.1416*MaximumDistance_Biochar[4][0]*Length_Biochar[4][0]*Length_Biochar[4][0]/8
Volume_Biochar[5][0] = 2*3.1416*MaximumDistance_Biochar[5][0]*Length_Biochar[5][0]*Length_Biochar[5][0]/8

'''********************************************Biochar & Concrete percentage and placement rule definition**********************'''
# Biochar percentage
BiocharFraction = 0.09
# Concrete geometric dimension definition:
ConcreteLength = 50                  
ConcreteWidth = 50
ConcreteHeight = 50
# Biochar volume
TotalAggVolume = 0      # Volume before placement
ConcreteVolume = ConcreteHeight*ConcreteLength*ConcreteWidth
TargetVolume = ConcreteHeight*ConcreteLength*ConcreteWidth*BiocharFraction      # Biochar target volume 
# Preprocess Biochar placement:
TryNumBiochar = 5000          # Maximum number of Biochar placement iterations
IterLimite = 10000        # Maximum number of position iteration
# Biochar placement point definition:
Biopoint = []
Biopoints = []
BioData = []
BioMesh = []
# Position tolerance:
BioGap = 1.03                                           # Distance between aggregates 
BaseGapMax = 0.97                                       # Distance between aggregate and bottom
BaseGapMin = 0.03                                       # Distance between aggregate and top
'''************************************************Biochar placement******************************'''
## Place using a for loop, using random sampling
BioGenerate = 0
for BioNum in range(TryNumBiochar):
    # Define a random number to select the number from the library:
    BioLib = random.randint(1,Biochar_Num)
    mdb.models[Modelname_Lib].rootAssembly.Instance(dependent=ON, name='Biochar-'+str(BioNum), 
        part=mdb.models[Modelname_Lib].parts['Biochar-'+str(BioLib)])
    for iter in range(IterLimite):
        if iter < IterLimite - 1 :
            # Generate random centroid coordinates xyz
            x1=np.random.uniform(0+MaximumDistance_Biochar[BioLib-1][0]/2.0,ConcreteLength-MaximumDistance_Biochar[BioLib-1][0]/2.0)     # x-coordinate of the centroid, mm   
            y1=np.random.uniform(0+MaximumDistance_Biochar[BioLib-1][0]/2.0,ConcreteWidth-MaximumDistance_Biochar[BioLib-1][0]/2.0)      # y-coordinate of the centroid, mm  
            z1=np.random.uniform(0+MaximumDistance_Biochar[BioLib-1][0]/2.0,ConcreteHeight-MaximumDistance_Biochar[BioLib-1][0]/2.0)     # z-coordinate of the centroid ,mm    
            # Store the random xyz and the maximum radius of the generated Biochar:
            Biopoint = (x1,y1,z1,MaximumDistance_Biochar[BioLib-1][0]/2)  #x y z r(radius)
            if len(Biopoints) == 0:
                if interact_judgement(Aggpoints,Biopoint):
                    Biopoints.append(Biopoint)
                    
                    # Record the total generated volume:
                    TotalAggVolume = TotalAggVolume + Volume_Biochar[BioLib-1][0]
                    Process = TotalAggVolume*100/TargetVolume
                    Friction = TotalAggVolume*100/ConcreteVolume
                    print('Success!')
                    print('Current Volume = {:.3f}, Process = {:.3f}%, Friction = {:.3f} %Current Step = {}'.format(TotalAggVolume,Process,Friction,BioNum))
                    print('Current Iter time = {:.3f}'.format(iter))
                    # Random rotation coordinates:
                    axisPoint = (0,0,0)
                    xP=np.random.uniform(-1,1)
                    yP=np.random.uniform(-1,1)
                    zP=np.random.uniform(-1,1)
                    AxisDirection = (xP,yP,zP)
                    # Random rotation angle:
                    RotationAngle = random.randint(0,360)

                    BioData.append([BioNum+1,x1,y1,z1,xP,yP,zP,RotationAngle,MaximumDistance_Biochar[BioLib-1][0]/2,BioLib])
                    BioMesh.append([BioGenerate+1,x1,y1,z1,xP,yP,zP,RotationAngle,MaximumDistance_Biochar[BioLib-1][0]/2,BioLib])
                    AxisDirection = (BioData[BioGenerate][4],BioData[BioGenerate][5],BioData[BioGenerate][6])

                    # Perform random rotation:
                    mdb.models[Modelname_Lib].rootAssembly.rotate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),axisPoint=axisPoint, 
                        axisDirection=AxisDirection, angle=BioData[BioGenerate][7])
                    # Translate after determining the coordinates
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(BioData[BioGenerate][1],0.0, 0.0))
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(0.0, BioData[BioGenerate][2],0.0 ))    
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(0.0, 0.0, BioData[BioGenerate][3]))
                    # Biochar generation counter increases by one:
                    BioGenerate = BioGenerate + 1
                    break
            elif interact_judgement(Biopoints,Biopoint):
                if interact_judgement(Aggpoints,Biopoint):
                    # Store the qualified aggregate point coordinates and radius in the total aggregate accumulation
                    Biopoints.append(Biopoint)
                    # Store the generated aggregate in the total aggregate data
                    
                    # Record the total generated volume:
                    TotalAggVolume = TotalAggVolume + Volume_Biochar[BioLib-1][0]
                    Process = TotalAggVolume*100/TargetVolume
                    Friction = TotalAggVolume*100/ConcreteVolume
                    print('Success!')
                    print('Current Volume = {:.3f}, Process = {:.3f}%, Friction = {:.3f}% Current Step = {}'.format(TotalAggVolume,Process,Friction,BioNum))
                    print('Current Iter time = {:.3f}'.format(iter))
                    ## Perform random rotation before translation, with the origin as the first reference point (centroid), and the subsequent coordinates as a fluctuation range
                    # Random rotation coordinates:
                    axisPoint = (0,0,0)
                    xP=np.random.randint(-1,1)
                    yP=np.random.randint(-1,1)
                    zP=np.random.randint(-1,1)
                    AxisDirection = (xP,yP,zP)
                    # Random rotation angle:
                    RotationAngle = random.randint(0,360)

                    BioData.append([BioNum+1,x1,y1,z1,xP,yP,zP,RotationAngle,MaximumDistance_Biochar[BioLib-1][0]/2,BioLib])
                    BioMesh.append([BioGenerate+1,x1,y1,z1,xP,yP,zP,RotationAngle,MaximumDistance_Biochar[BioLib-1][0]/2,BioLib])
                    AxisDirection = (BioData[BioGenerate][4],BioData[BioGenerate][5],BioData[BioGenerate][6])

                    # Perform random rotation:
                    mdb.models[Modelname_Lib].rootAssembly.rotate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),axisPoint=axisPoint, 
                        axisDirection=AxisDirection, angle=BioData[BioGenerate][7])
                    # Translate after determining the coordinates
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(BioData[BioGenerate][1],0.0, 0.0))
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(0.0, BioData[BioGenerate][2],0.0 ))    
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-'+str(BioData[BioGenerate][0]-1),),vector=(0.0, 0.0, BioData[BioGenerate][3]))
                    # Biochar generation counter increases by one:
                    BioGenerate = BioGenerate + 1
                    break
        else:
            del mdb.models[Modelname_Lib].rootAssembly.features['Biochar-'+str(BioNum)]
            print('Fail!')
            print('Generate Failure')
            #BioNum = BioNum -1
            break
    if (TotalAggVolume - TargetVolume) >= 0.001:
        break

'''Generate aggregates********************************************************************************'''
# create in Abaqus
b=AggData

## Single aggregate definition:
NumLayers = 5
NumMiddleLayer = 8
NumNode1st = 1
NumNode2nd = 8
NumNode3rd = 8
NumNode4th = 8
NumNode5th = 1
AngleGap = 180.0 / (NumLayers-1)                          # The angle size for equal division
NumTotNodes = NumNode1st+NumNode2nd+NumNode3rd+NumNode4th+NumNode5th

### Multiple polygonal aggregates collection (number of aggregates, aggregate index, and node coordinates collection):
NumAgg = len(AggData)                                             # Number of polygonal aggregates
# Information for each aggregate in the collection:
RandomNodesPor = np.zeros((NumAgg,NumTotNodes,4))       # Spherical coordinates (point index, radius, theta, phi)
RandomNodesCar = np.zeros((NumAgg,NumTotNodes,4))       # Cartesian coordinates (point index, x, y, z)
# Face definition when converting SHELL to solid:
NumFaces = NumMiddleLayer*2 + (NumLayers-2-1) * NumMiddleLayer*2
Faces = np.zeros((NumFaces,13))                         # Saved areas, recorded in the order (face index, first point (3), second point (3), third point (3), centroid (3))
FacesITZ = np.zeros((NumFaces,13)) 

AggRadius = np.zeros((NumAgg,1))                          # Aggregate particle size collection
AggCentroid = np.zeros((NumAgg,4))                      # Each generated aggregate's index and centroid coordinates (xyz)


# ITZ region modeling:
RandomNodesITZCar = np.zeros((NumAgg,NumTotNodes,4))       # ITZ region Cartesian coordinates (point index, x, y, z)
AggITZData = []                                            # ITZ data

for Aggnum in range(len(AggData)):
    FaceFlage = 0                    # Record the generated face index, starting from 0, total of 48 faces
    FaceFlageITZ = 0
    # Create the part for the Aggregate
    SizeLength=200
    mdb.models[Modelname].ConstrainedSketch(name='__profile__',sheetSize=1)
    mdb.models[Modelname].sketches['__profile__'].rectangle(point1=(-SizeLength/2, SizeLength/2), 
    point2=(SizeLength/2, -SizeLength/2))
    mdb.models[Modelname].Part(dimensionality=THREE_D, name='PolyAgg-'+str(Aggnum), type=DEFORMABLE_BODY)
    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].BaseSolidExtrude(depth=SizeLength/2, sketch=
    mdb.models[Modelname].sketches['__profile__'])
    del mdb.models[Modelname].sketches['__profile__']
    del mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].features['Solid extrude-1']
    # Create the part for the ITZ region
    mdb.models[Modelname].ConstrainedSketch(name='__profile__',sheetSize=1)
    mdb.models[Modelname].sketches['__profile__'].rectangle(point1=(-SizeLength/2, SizeLength/2), 
    point2=(SizeLength/2, -SizeLength/2))
    mdb.models[Modelname].Part(dimensionality=THREE_D, name='PolyAggITZregion-'+str(Aggnum), type=DEFORMABLE_BODY)
    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].BaseSolidExtrude(depth=SizeLength/2, sketch=
    mdb.models[Modelname].sketches['__profile__'])
    del mdb.models[Modelname].sketches['__profile__']
    del mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].features['Solid extrude-1']
    RadiusVio = 0.15         # Node fluctuation coefficient along the radius (relative to the radius [0,1])
    ## Angle fluctuation can affect convergence, it is recommended to be around 12, or adjust the fluctuation range, otherwise, it may fail to generate
    AngleVio = 9           # Node fluctuation range along the angle (theta & phi)
    RadiusVioAdjust = 1     # Determine whether to fluctuate along the radial direction [0 or 1]
    AngleVioAdjust = 1      # Determine whether to fluctuate along the angle direction [0 or 1]
    radius = AggData[Aggnum][4]  # Radius is pre-generated
    for layer in range(NumLayers):                      # Random distribution of nodes for each layer
        theta = layer * AngleGap/180*math.pi       # Use theta to divide each layer, for 5 layers of aggregate, each layer is 45°
        ############################################## First layer (top layer) node [1 node] #############################################
        if abs(layer-0) <= 1e-10:                           # Convergence judgment
            for Node in range(NumNode1st):
                #Coe = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                phi = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][layer][0] = layer+1    # Top layer 1st point index: 1
                RandomNodesPor[Aggnum][layer][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           # Top layer 1st point control radius
                RandomNodesPor[Aggnum][layer][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Top layer 1st point's Theta
                RandomNodesPor[Aggnum][layer][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Top layer 1st point's Phi
        ############################################## Last layer node [1 node] #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            for Node in range(NumNode5th):
                phi = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][0] = (layer-1)*NumNode2nd+NumNode1st+1                # Bottom layer 26th point index: 1
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])           # Bottom layer 26th point control radius
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Bottom layer 26th point's Theta
                RandomNodesPor[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])    # Bottom layer 26th point's Phi
        ##################################################### Middle three layers [8 nodes] ###############################################
        else:
            for Node in range(NumMiddleLayer):
                # Coe = (-1+2*np.random.random((1,1))[0][0])  # Random fluctuation coefficient
                phi = Node * AngleGap/180 * math.pi         # Rotation angle within each layer
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][0] = Node+(layer-1)*NumNode2nd+NumNode1st+1                # Middle layer 26th point index: 1
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1] = radius + RadiusVioAdjust*RadiusVio*(-1+2*np.random.random((1,1))[0][0])          # Middle layer 26th point control radius
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2] = theta + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])  # Middle layer 26th point's Theta
                RandomNodesPor[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3] = phi + AngleVioAdjust*AngleVio/180.0*math.pi*(-1+2*np.random.random((1,1))[0][0])   # Middle layer 26th point's Phi
    '''***********************************************************Convert spherical coordinates (por) to Cartesian coordinates (Car)*********************************************************'''
    # x = r*sin(theta)*cos(phi);      y = r*sin(theta)*sin(phi);        z = r*cos(theta)
    for Node in range(NumTotNodes):
        RandomNodesCar[Aggnum][Node][0] = RandomNodesPor[Aggnum][Node][0]
        RandomNodesCar[Aggnum][Node][1] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.cos(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][2] = RandomNodesPor[Aggnum][Node][1]*math.sin(RandomNodesPor[Aggnum][Node][2])*np.sin(RandomNodesPor[Aggnum][Node][3])
        RandomNodesCar[Aggnum][Node][3] = RandomNodesPor[Aggnum][Node][1]*math.cos(RandomNodesPor[Aggnum][Node][2])
    '''***********************************************************Connect the 26 points into lines*********************************************************'''
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
    '''***********************************************************Points -> Lines: Connecting the points in the previous layer with the next layer: 123, 145, 156, 167, 178, 189, 192*********************************************************'''
    for layer in range(NumLayers):
        ############################################# Connect the top layer (1 point) with the next layer #############################################
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][layer][1]
            a1[1] = RandomNodesCar[Aggnum][layer][2]
            a1[2] = RandomNodesCar[Aggnum][layer][3]
            # Convert list coordinates to tuple [first vertex, connect with all nodes below]
            b1 = tuple(a1)
            # Second layer of 8 points loop, connect with the top point
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
                # After completing one loop, the last connected face is 192, the second point needs to be updated to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 as the second point, index is 2
                    a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        ############################################# Connect the bottom layer [1 point] with the previous layer (second-to-last layer [8 points]): #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]
            # Convert list coordinates to tuple [first vertex, connect with all nodes above]
            b1 = tuple(a1)
            # Second-to-last layer of 8 points, connect with the bottom point:
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
                # After completing one loop, the last connected face is 192, the second point needs to be updated to b3:
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 as the second point, index is 2
                    a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        #####********************************* The middle layer must be related to the layer!!! Otherwise, there will be issues generating surfaces *********************************#####
        elif abs(layer - 3) >= 1e-10:
            ## The middle layer is composed of two triangles, b and c each controlling one triangle
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
                # Generate the triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the upper layer +1 node
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
                # After completing one loop, update the two triangle faces **Node 9 does not need to be updated, so b1 and c1 do not need to be updated:
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    ## Connect 9 - 17 -10 (at this point NODE = 7)
                    a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      # Node 17
                    a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      # Node 10
                    ## Connect 9 - 10 -2 (at this point NODE = 7)                  
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
                # No need to connect c1c2c3 because after connecting the two triangles, a quadrilateral can be formed, connecting would cause an error
    '''*************************Lines -> Faces: Merge the previously closed line segments into faces, using the coveredge and findat methods [need to define a new variable to store all generated faces for assigning solid attributes later]*********************************************************'''
    for layer in range(NumLayers):
        ############################################# Close the top layer face #############################################
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][layer][1]
            a1[1] = RandomNodesCar[Aggnum][layer][2]
            a1[2] = RandomNodesCar[Aggnum][layer][3]
            # Convert list coordinates to tuple [first vertex, connect with all nodes below]
            b1 = tuple(a1)
            # Second layer of 8 points loop, connect with the top point
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
                # After completing one loop, the last connected face is 192, the second point needs to be updated to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 as the second point, index is 2
                    a3[0] = RandomNodesCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS's coveredge method to close the line segments into a face; find the middle point of the three sides of the triangle through findat, and then generate the face:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the face (8 faces), record the centroid coordinates of each face (extract the three vertices of the triangle, and then calculate the centroid coordinates):
                Faces[FaceFlage][0] = FaceFlage + 1                 # Index starts from 1, record the first face
                # First point coordinates
                Faces[FaceFlage][1] = a1[0]                         # First point x-coordinate
                Faces[FaceFlage][2] = a1[1]                         # First point y-coordinate
                Faces[FaceFlage][3] = a1[2]                         # First point z-coordinate
                # Second point coordinates
                Faces[FaceFlage][4] = a2[0]                         # Second point x-coordinate
                Faces[FaceFlage][5] = a2[1]                         # Second point y-coordinate
                Faces[FaceFlage][6] = a2[2]                         # Second point z-coordinate
                # Third point coordinates
                Faces[FaceFlage][7] = a3[0]                         # Third point x-coordinate
                Faces[FaceFlage][8] = a3[1]                         # Third point y-coordinate
                Faces[FaceFlage][9] = a3[2]                         # Third point z-coordinate
                # Centroid coordinates (need to be extracted to help locate the face)
                Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x-coordinate
                Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y-coordinate
                Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z-coordinate
                # Each time a face is completed, increment the face index counter by 1:
                FaceFlage += 1                                      # Update face index counter
        ############################################# Close the bottom layer face: #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesCar[Aggnum][NumTotNodes-1][3]
            # Convert list coordinates to tuple [first vertex, connect with all nodes above]
            b1 = tuple(a1)
            # Second-to-last layer of 8 points, connect with the bottom point:
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
                # After completing one loop, the last connected face is 192, the second point needs to be updated to b3:
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 as the second point, index is 2
                    a3[0] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS's coveredge method to close the line segments into a face; find the middle point of the three sides of the triangle through findat, and then generate the face:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the face (8 faces), record the centroid coordinates of each face (extract the three vertices of the triangle, and then calculate the centroid coordinates):
                Faces[FaceFlage][0] = FaceFlage + 1                 # Index starts from 1, record the first face
                # First point coordinates
                Faces[FaceFlage][1] = a1[0]                         # First point x-coordinate
                Faces[FaceFlage][2] = a1[1]                         # First point y-coordinate
                Faces[FaceFlage][3] = a1[2]                         # First point z-coordinate
                # Second point coordinates
                Faces[FaceFlage][4] = a2[0]                         # Second point x-coordinate
                Faces[FaceFlage][5] = a2[1]                         # Second point y-coordinate
                Faces[FaceFlage][6] = a2[2]                         # Second point z-coordinate
                # Third point coordinates
                Faces[FaceFlage][7] = a3[0]                         # Third point x-coordinate
                Faces[FaceFlage][8] = a3[1]                         # Third point y-coordinate
                Faces[FaceFlage][9] = a3[2]                         # Third point z-coordinate
                # Centroid coordinates (need to be extracted to help locate the face)
                Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x-coordinate
                Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y-coordinate
                Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z-coordinate
                # Each time a face is completed, increment the face index counter by 1:
                FaceFlage += 1                                      # Update face index counter
        ############################################# Close the middle layer face: ########################################################
        elif abs(layer - 3) >= 1e-10:
            ## The middle layer is composed of two triangles, b and c each controlling one triangle
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
                # Generate the triangle controlled by c, starting from 2-11-3, so c1 = b1; c3 = b3; c2 is the upper layer +1 node
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
                # After completing one loop, update the two triangle faces **Node 9 does not need to be updated, so b1 and c1 do not need to be updated:
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    ## Connect 9 - 17 -10 (at this point NODE = 7)
                    a2[0] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      # Node 17
                    a3[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      # Node 10
                    ## Connect 9 - 10 -2 (at this point NODE = 7)                  
                    a5[0] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          # Node 2
                    a6[0] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          # Node 10
                ############ Use ABAQUS's coveredge method to close the line segments into a face; find the middle point of the three sides of the triangle through findat, and then generate the face:
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                    # For middle layer faces, you may not need to connect all, but when generating the face, it's done as a triangle, so generate for both b and c points
                mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c2[0]+c3[0])/2, (c2[1]+c3[1])/2, (c2[2]+c3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((c3[0]+c1[0])/2, (c3[1]+c1[1])/2, (c3[2]+c1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the face (16 faces), record the centroid coordinates of each face (extract the three vertices of the triangle, and then calculate the centroid coordinates):
                Faces[FaceFlage][0] = FaceFlage + 1                 # Index starts from 1, record the first face
                # The face controlled by b
                # First point coordinates
                Faces[FaceFlage][1] = a1[0]                         # First point x-coordinate
                Faces[FaceFlage][2] = a1[1]                         # First point y-coordinate
                Faces[FaceFlage][3] = a1[2]                         # First point z-coordinate
                # Second point coordinates
                Faces[FaceFlage][4] = a2[0]                         # Second point x-coordinate
                Faces[FaceFlage][5] = a2[1]                         # Second point y-coordinate
                Faces[FaceFlage][6] = a2[2]                         # Second point z-coordinate
                # Third point coordinates
                Faces[FaceFlage][7] = a3[0]                         # Third point x-coordinate
                Faces[FaceFlage][8] = a3[1]                         # Third point y-coordinate
                Faces[FaceFlage][9] = a3[2]                         # Third point z-coordinate
                # Centroid coordinates (need to be extracted to help locate the face)
                Faces[FaceFlage][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x-coordinate
                Faces[FaceFlage][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y-coordinate
                Faces[FaceFlage][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z-coordinate
                # Each time a face is completed, increment the face index counter by 1:
                FaceFlage += 1                                      # Update face index counter
                # The face controlled by c
                # First point coordinates
                Faces[FaceFlage][1] = a4[0]                         # First point x-coordinate
                Faces[FaceFlage][2] = a4[1]                         # First point y-coordinate
                Faces[FaceFlage][3] = a4[2]                         # First point z-coordinate
                # Second point coordinates
                Faces[FaceFlage][4] = a5[0]                         # Second point x-coordinate
                Faces[FaceFlage][5] = a5[1]                         # Second point y-coordinate
                Faces[FaceFlage][6] = a5[2]                         # Second point z-coordinate
                # Third point coordinates
                Faces[FaceFlage][7] = a6[0]                         # Third point x-coordinate
                Faces[FaceFlage][8] = a6[1]                         # Second point y-coordinate
                Faces[FaceFlage][9] = a6[2]                         # Third point z-coordinate
                # Centroid coordinates (need to be extracted to help locate the face)
                Faces[FaceFlage][10] = (a4[0]+a5[0]+a6[0])/3                        # Centroid x-coordinate
                Faces[FaceFlage][11] = (a4[1]+a5[1]+a6[1])/3                        # Centroid y-coordinate
                Faces[FaceFlage][12] = (a4[2]+a5[2]+a6[2])/3                        # Centroid z-coordinate
                # Each time a face is completed, increment the face index counter by 1:
                FaceFlage += 1                                      # Update face index counter
    '''*******************Surface -> Solid: Record the centroid of all surfaces generated for a single aggregate into the facelist collection and find them in ABAQUS*********************************************************'''
    # All surfaces have been generated, and the centroid coordinates of the surfaces have been saved in Faces[10] (x coordinate); Faces[11] (y coordinate); Faces[12] (z coordinate)
    Facelist = []       # Define an empty set to store all the surfaces that need to be turned into solids
    # Loop through all surfaces and store them in facelist:
    for i in range(int(Faces.shape[0])):
        Facelist.append( mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].faces.findAt((Faces[i][10], Faces[i][11],Faces[i][12]), ) )
    # Use the cell command in ABAQUS to convert surfaces into solids:
    mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].AddCells(faceList=Facelist)
    # Place each part into the assembly as an independent instance
    mdb.models[Modelname].rootAssembly.Instance(dependent=ON, name='PolyAgg-'+str(Aggnum), 
        part=mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)])

    '''*********************************Preparation of coordinate parameters for aggregate placement determination*****************************************************************************'''#
    ## Use the getMassProperties command in ABAQUS to read the centroid of each generated aggregate:
    AggCentroid[Aggnum][0] = Aggnum+1                                                                # Record aggregate number
    AggCentroid[Aggnum][1]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][0]  # Centroid x coordinate
    AggCentroid[Aggnum][2]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][1]  # Centroid y coordinate
    AggCentroid[Aggnum][3]=mdb.models[Modelname].rootAssembly.getMassProperties()['volumeCentroid'][2]  # Centroid z coordinate
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
        # Use SEPARATE merge type, otherwise wire will fail to generate
        mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].WirePolyLine(mergeType=SEPARATE, meshable=
            ON, points=(b2,b1))
        # Create an offset point (0.1 position)
        '''Generate ITZ region (ITZ thickness changes with the radius of the aggregate, but for placement efficiency, it is generated inward)'''
        mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].DatumPointByEdgeParam(edge=
            mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2)),
            parameter=0.95, isDependent=False)
        # Delete the generated line for subsequent calculations 
        del mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].features['Wire-'+str(33)]
    # Loop through all ID numbers and store them in DatumsPoints    
    for datum in range(NumTotNodes):   
        # Get all DatumPoints ID numbers 
        DatumsPoints.append(mdb.models[Modelname].parts['PolyAgg-'+str(Aggnum)].datums[84+datum*2])
    # Use the pointOn command to save the coordinates of the left and right points in the ITZ node region:
    for Node in range(NumTotNodes):
        RandomNodesITZCar[Aggnum][Node][0] = RandomNodesPor[Aggnum][Node][0]
        RandomNodesITZCar[Aggnum][Node][1] = DatumsPoints[Node].pointOn[0]
        RandomNodesITZCar[Aggnum][Node][2] = DatumsPoints[Node].pointOn[1]
        RandomNodesITZCar[Aggnum][Node][3] = DatumsPoints[Node].pointOn[2]
    '''***********************************************************Connect 26 points into lines*********************************************************'''
    # Define empty lists to allocate points, then define an empty tuple to convert lists into tuples for connection, c is used for middle layer connection
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
    '''***********************************************************Points -> Lines: Connect the points of the previous layer to the next layer: 123, 145, 156, 167, 178, 189, 192*********************************************************'''
    for layer in range(NumLayers):
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesITZCar[Aggnum][layer][1]
            a1[1] = RandomNodesITZCar[Aggnum][layer][2]
            a1[2] = RandomNodesITZCar[Aggnum][layer][3]
            # Convert list coordinates into tuples [first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # Loop through the 8 points of the second layer, connecting to the vertex
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
                # After one loop, the last connected surface is 192, so update the second point to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 to the second node, numbered 2
                    a3[0] = RandomNodesITZCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        ############################################# Connect the bottom layer point [1 point] to the upper layer (the second to last layer [8 points]): #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesITZCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesITZCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesITZCar[Aggnum][NumTotNodes-1][3]
            # Convert list coordinates into tuples [first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # Connect the 8 points of the second to last layer to the bottom point:
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
                # After one loop, the last connected surface is 192, so update the second point to b3:
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 to the second node, numbered 2
                    a3[0] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use ABAQUS to connect three points into three lines, forming a closed surface:
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
        #####*********************************The middle layer must establish a relationship with the layer!!! Otherwise, the surface generation will have issues *********************************#####
        elif abs(layer - 3) >= 1e-10:
            ## The middle layer is composed of two triangles, with b and c each controlling one triangle
            for Node in range(NumMiddleLayer):
                # Generate the triangle controlled by b, starting with 2-10-11
                a1[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a1[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a1[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                b1 = tuple(a1)          # Node 2
                a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)          # Node 10
                a3[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)          # Node 11
                # Generate the triangle controlled by c, starting with 2-11-3, so c1 = b1; c3 = b3; c2 is a single upper layer+1 node
                a4[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a4[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a4[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                c1 = tuple(a4)          # Node 2
                a5[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
                a5[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
                a5[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
                c2 = tuple(a5)          # Node 3
                a6[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a6[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a6[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                c3 = tuple(a6)          # Node 11
                # After one loop, update the two triangle surfaces **node 9 does not need to be updated, so b1 c1 do not need to be updated:
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    ## Connect 9 - 17 -10 (at this time, NODE = 7)
                    a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      # Node 17
                    a3[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      # Node 10
                    ## Connect 9 - 10 -2 (at this time, NODE = 7)                  
                    a5[0] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          # Node 2
                    a6[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          # Node 10
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].WirePolyLine(mergeType=IMPRINT, meshable=
                    ON, points=((b1,b2), (b2, b3), (b3, b1)))
    '''*************************Line -> Surface: Merge the previously closed line segments into surfaces, using the coveredge and findat methods [a new variable needs to be defined to store all generated surfaces for assigning solid properties later]*********************************************************'''
    for layer in range(NumLayers):
        ############################################# Close the top surface #############################################
        if abs(layer - 0) <= 1e-10:
            a1[0] = RandomNodesITZCar[Aggnum][layer][1]
            a1[1] = RandomNodesITZCar[Aggnum][layer][2]
            a1[2] = RandomNodesITZCar[Aggnum][layer][3]
            # Convert list coordinates into tuples [first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # Loop through the 8 points of the second layer, connecting to the vertex
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
                # After one loop, the last connected surface is 192, so update the second point to b3:
                if abs(Node - (NumNode2nd-1)) <= 1e-10:
                    a2[0] = RandomNodesITZCar[Aggnum][Node+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 to the second node, numbered 2
                    a3[0] = RandomNodesITZCar[Aggnum][NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use the coveredge method in ABAQUS to close the line segments into surfaces; use findat to find the edges where the midpoint of the triangle sides is located, and then generate the surface:
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the surface (8 in total), record the centroid coordinates of each surface (extract the three vertices of the triangle, then calculate the midpoint coordinates):
                FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start recording the first surface with number 1
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
                # Centroid coordinates (need to be extracted to help locate the surface)
                FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x coordinate
                FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y coordinate
                FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z coordinate
                # After finishing each surface, increment the surface number counter by 1:
                FaceFlageITZ += 1                                      # Update the surface number counter
        ############################################# Close the bottom surface: #############################################
        elif abs(layer - (NumLayers-1)) <= 1e-10:
            a1[0] = RandomNodesITZCar[Aggnum][NumTotNodes-1][1]
            a1[1] = RandomNodesITZCar[Aggnum][NumTotNodes-1][2]
            a1[2] = RandomNodesITZCar[Aggnum][NumTotNodes-1][3]
            # Convert list coordinates into tuples [first vertex, connected to all nodes below]
            b1 = tuple(a1)
            # Connect the 8 points of the second to last layer to the bottom point:
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
                # After one loop, the last connected surface is 192, so update the second point to b3:
                if abs(Node - (NumNode4th-1)) <= 1e-10:
                    a2[0] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node + (layer-2)*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)
                    # Update the last connected point b3 to the second node, numbered 2
                    a3[0] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][(layer-2)*NumMiddleLayer+NumNode1st][3]
                    b3 = tuple(a3)
                ############ Use the coveredge method in ABAQUS to close the line segments into surfaces; use findat to find the edges where the midpoint of the triangle sides is located, and then generate the surface:
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the surface (8 in total), record the centroid coordinates of each surface (extract the three vertices of the triangle, then calculate the midpoint coordinates):
                FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start recording the first surface with number 1
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
                # Centroid coordinates (need to be extracted to help locate the surface)
                FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x coordinate
                FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y coordinate
                FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z coordinate
                # After finishing each surface, increment the surface number counter by 1:
                FaceFlageITZ += 1                                      # Update the surface number counter
        ############################################# Close the middle layer surfaces: ########################################################
        elif abs(layer - 3) >= 1e-10:
            ## The middle layer is composed of two triangles, with b and c each controlling one triangle
            for Node in range(NumMiddleLayer):
                # Generate the triangle controlled by b, starting with 2-10-11
                a1[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a1[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a1[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                b1 = tuple(a1)          # Node 2
                a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                b2 = tuple(a2)          # Node 10
                a3[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a3[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a3[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                b3 = tuple(a3)          # Node 11
                # Generate the triangle controlled by c, starting with 2-11-3, so c1 = b1; c3 = b3; c2 is a single upper layer+1 node
                a4[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][1]
                a4[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][2]
                a4[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st][3]
                c1 = tuple(a4)          # Node 2
                a5[0] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][1]
                a5[1] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][2]
                a5[2] = RandomNodesITZCar[Aggnum][Node + (layer-1)*NumMiddleLayer+NumNode1st+1][3]
                c2 = tuple(a5)          # Node 3
                a6[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][1]
                a6[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][2]
                a6[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st+1][3]
                c3 = tuple(a6)          # Node 11
                # After one loop, update the two triangle surfaces **node 9 does not need to be updated, so b1 c1 do not need to be updated:
                if abs(Node - (NumMiddleLayer-1)) <= 1e-10:
                    ## Connect 9 - 17 -10 (at this time, NODE = 7)
                    a2[0] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][1]
                    a2[1] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][2]
                    a2[2] = RandomNodesITZCar[Aggnum][Node + layer*NumMiddleLayer+NumNode1st][3]
                    b2 = tuple(a2)      # Node 17
                    a3[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a3[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a3[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    b3 = tuple(a3)      # Node 10
                    ## Connect 9 - 10 -2 (at this time, NODE = 7)                  
                    a5[0] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][1]
                    a5[1] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][2]
                    a5[2] = RandomNodesITZCar[Aggnum][(layer-1)*NumMiddleLayer+NumNode1st][3]
                    c2 = tuple(a5)          # Node 2
                    a6[0] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][1]
                    a6[1] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][2]
                    a6[2] = RandomNodesITZCar[Aggnum][NumMiddleLayer*layer+NumNode1st][3]
                    c3 = tuple(a6)          # Node 10
                ############ Use the coveredge method in ABAQUS to close the line segments into surfaces; use findat to find the edges where the midpoint of the triangle sides is located, and then generate the surface:
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b1[0]+b2[0])/2, (b1[1]+b2[1])/2, (b1[2]+b2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b2[0]+b3[0])/2, (b2[1]+b3[1])/2, (b2[2]+b3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((b3[0]+b1[0])/2, (b3[1]+b1[1])/2, (b3[2]+b1[2])/2), )),  # Third point
                    tryAnalytical=True)
                    # For middle layer surfaces, the connection does not need to be fully considered, but surface generation follows a triangular pattern, so surfaces must be generated for both b and c points
                mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].CoverEdges(edgeList=(
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c1[0]+c2[0])/2, (c1[1]+c2[1])/2, (c1[2]+c2[2])/2), ),   # First point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c2[0]+c3[0])/2, (c2[1]+c3[1])/2, (c2[2]+c3[2])/2), ),   # Second point
                    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].edges.findAt(((c3[0]+c1[0])/2, (c3[1]+c1[1])/2, (c3[2]+c1[2])/2), )),  # Third point
                    tryAnalytical=True)
                ############ After generating the surface (16 in total), record the centroid coordinates of each surface (extract the three vertices of the triangle, then calculate the midpoint coordinates):
                FacesITZ[FaceFlageITZ][0] = FaceFlageITZ + 1                 # Start recording the first surface with number 1
                # Surface controlled by the triangle b
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
                # Centroid coordinates (need to be extracted to help locate the surface)
                FacesITZ[FaceFlageITZ][10] = (a1[0]+a2[0]+a3[0])/3                        # Centroid x coordinate
                FacesITZ[FaceFlageITZ][11] = (a1[1]+a2[1]+a3[1])/3                        # Centroid y coordinate
                FacesITZ[FaceFlageITZ][12] = (a1[2]+a2[2]+a3[2])/3                        # Centroid z coordinate
                # After finishing each surface, increment the surface number counter by 1:
                FaceFlageITZ += 1                                      # Update the surface number counter
                # Surface controlled by the triangle c
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
                # Centroid coordinates (need to be extracted to help locate the surface)
                FacesITZ[FaceFlageITZ][10] = (a4[0]+a5[0]+a6[0])/3                        # Centroid x coordinate
                FacesITZ[FaceFlageITZ][11] = (a4[1]+a5[1]+a6[1])/3                        # Centroid y coordinate
                FacesITZ[FaceFlageITZ][12] = (a4[2]+a5[2]+a6[2])/3                        # Centroid z coordinate
                # After finishing each surface, increment the surface number counter by 1:
                FaceFlageITZ += 1                                      # Update the surface number counter
    '''*******************Surface -> Solid: Record the centroid of all surfaces generated for a single aggregate into the facelist collection and find them in ABAQUS*********************************************************'''
    FacelistITZ = []       # Define an empty set to store all the surfaces that need to be turned into solids
    # Loop through all surfaces and store them in facelist:
    for k in range(int(FacesITZ.shape[0])):
        FacelistITZ.append( mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].faces.findAt((FacesITZ[k][10], FacesITZ[k][11],FacesITZ[k][12]), ) )
    # Use the cell command in ABAQUS to convert surfaces into solids:
    mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)].AddCells(faceList=FacelistITZ)
    # Place each part into the assembly as an independent instance
    mdb.models[Modelname].rootAssembly.Instance(dependent=ON, name='PolyAggITZregion-'+str(Aggnum), 
        part=mdb.models[Modelname].parts['PolyAggITZregion-'+str(Aggnum)])

    print('Current Aggregate = {},Total Aggregtas = {}'.format(Aggnum,len(AggData)))

'''**********Translate polyhedral aggregates in the ASSEMBLY module****************************************************************************************************'''
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


'''***********************Generate Concrete Outline**********************************************************************'''
# Command to generate the concrete outline
xmin=0
xmax=ConcreteLength
ymin=0
ymax=ConcreteWidth
zmin=0
zmax=ConcreteHeight
zlength=abs(zmax-zmin)
ConcreteSketch = myModel.ConstrainedSketch(name='concretecube',sheetSize=200)
ConcreteSketch.rectangle(point1=(xmin, ymin), point2=(xmax, ymax))
ConcretePart = myModel.Part(dimensionality=THREE_D, name='concretecube', type=DEFORMABLE_BODY)
myPart = ConcretePart.BaseSolidExtrude(depth=zlength, sketch=myModel.sketches['concretecube'])
del myModel.sketches['concretecube']
myModel.rootAssembly.Instance(name='concretecube', part=ConcretePart, dependent=ON)


'''*************Create a unified set for the generated ExtendAgg, Aggregate and merge them, then delete parts to save memory************'''
##Aggregate operation:
# Create a set for all Aggregates:
Aggregate=mdb.models[Modelname].rootAssembly.instances['PolyAgg-'+str(0)].cells[0:0]   
for i in range(len(AggData)): 
    Aggregate=Aggregate+mdb.models[Modelname].rootAssembly.instances['PolyAgg-'+str(AggData[i][0]-1)].cells
mdb.models[Modelname].rootAssembly.Set(cells=Aggregate, name='ExtendAggregates')
# Merge all AGGREGATE into one part:
AggregateInstance = []
for i in range(len(AggData)):
    AggregateInstance.append(mdb.models[Modelname].rootAssembly.instances['PolyAgg-' + str(AggData[i][0]-1)])
AggregateInstance = tuple(AggregateInstance)
mdb.models[Modelname].rootAssembly.InstanceFromBooleanMerge(name='ExtendAggregates', instances=AggregateInstance, keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
##ITZ operation
# Create a set for all ExtenAggregate:
ExtendAggregate=mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-'+str(0)].cells[0:0]   
for i in range(len(AggData)): 
    ExtendAggregate=ExtendAggregate+mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-'+str(AggData[i][0]-1)].cells
mdb.models[Modelname].rootAssembly.Set(cells=ExtendAggregate, name='InnerAggregates')
AggregateInstance = []
# Merge all ExtendAGGREGATE into one part:
ExtentAggregateInstance = []
for i in range(len(AggData)):
    ExtentAggregateInstance.append(mdb.models[Modelname].rootAssembly.instances['PolyAggITZregion-' + str(AggData[i][0]-1)])
ExtentAggregateInstance = tuple(ExtentAggregateInstance)
mdb.models[Modelname].rootAssembly.InstanceFromBooleanMerge(name='InnerAggregates', instances=ExtentAggregateInstance, keepIntersections=ON, originalInstances=SUPPRESS, domain=GEOMETRY)
## Delete all single parts of Aggregate and ExtendAggregate to free memory:
for i in range(len(AggData)):
    del mdb.models[Modelname].parts['PolyAgg-' + str(AggData[i][0]-1)]
    del mdb.models[Modelname].parts['PolyAggITZregion-' + str(AggData[i][0]-1)]


'''************************************Perform Boolean Difference Operation on Generated ExtendAggregate to Obtain ITZ Region of the Outline************************'''
# Perform Boolean operation to obtain the difference ITZ region
mdb.models[Modelname].rootAssembly.InstanceFromBooleanCut(name='ITZRegion', 
    instanceToBeCut=mdb.models[Modelname].rootAssembly.instances['ExtendAggregates-1'], 
    cuttingInstances=(mdb.models[Modelname].rootAssembly.instances['InnerAggregates-1'], ), 
    originalInstances=SUPPRESS)
# Restore the aggregate region
mdb.models[Modelname].rootAssembly.features['InnerAggregates-1'].resume()
