#! /user/bin/python
# -*- coding:UTF-8 -*-
# filename: 3DSpheres.py

from abaqus import *
from abaqusConstants import *
from caeModules import *
import os
import numpy as np
import math
import random
import matplotlib.pyplot as plt
from visualization import *
from odbAccess import *
from ssl import SSLSocket
from driverUtils import executeOnCaeStartup


session.journalOptions.setValues(replayGeometry=INDEX, recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)
Mdb()

myModel = mdb.Model(name='Biochar Concrete')
Modelname = myModel.name

# Rectangular region dimensions
ConcreteLength = 50.0  # mm, length of the concrete block
ConcreteWidth = 50.0   # mm, width of the concrete block
ConcreteHeight = 50.0  # mm, height of the concrete block

# Modeling based on aggregate percentage rather than count
ConcreteVolume = ConcreteLength * ConcreteWidth * ConcreteHeight  # Concrete volume
AggRatio = 0.21  # Aggregate proportion
TargetVolume = ConcreteVolume * AggRatio  # Target volume of aggregates
TotalAggVolume = 0.0  # Accumulated volume of placed aggregates
AggGap = 1.01
BaseGapMax = 0.99
BaseGapMin = 0.01

# AGG iteration method
AggLimit = 3  # Maximum number of AGG iterations
# Use Fuller gradation and CCR article as reference
Dmax_1 = 4.75  # Maximum aggregate size for the first gradation level
Dmin_1 = 2.36  # Minimum aggregate size for the first gradation level
IterLimit = 10000  # Maximum iteration limit

# Function to check if the first sphere is inside the boundary
def first_agg(point):
    x1, y1, z1, r1 = point
    x_z, x_f = x1 + r1, x1 - r1
    y_z, y_f = y1 + r1, y1 - r1
    z_z, z_f = z1 + r1, z1 - r1

    if (x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength or
        y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin * ConcreteWidth or
        z_z > BaseGapMax * ConcreteHeight or z_f < BaseGapMin * ConcreteHeight):
        return False
    return True

# Function to check if the newly generated sphere intersects with existing spheres
def interact_judgement(points, point):
    x1, y1, z1, r1 = point
    x_z, x_f = x1 + r1, x1 - r1
    y_z, y_f = y1 + r1, y1 - r1
    z_z, z_f = z1 + r1, z1 - r1

    if (x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength or
        y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin * ConcreteWidth or
        z_z > BaseGapMax * ConcreteHeight or z_f < BaseGapMin * ConcreteHeight):
        return False

    for ii in points:
        x2, y2, z2, r2 = ii
        # Calculate distance between centers
        distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        if distance < AggGap * (r1 + r2):
            return False
    return True

# Generate spheres
k = 0
Aggpoint = []  # Accumulated data of generated spheres
Aggpoints = []  # Temporary variable
AggData = []  # Data of spherical aggregates

'''Generate spherical aggregates'''
i = 0
for Aggnum in range(AggLimit):
  radius = pow((np.random.random((1, 1))[0][0] * (pow(Dmax_1, 0.5) - pow(Dmin_1, 0.5))) + pow(Dmin_1, 0.5), 2) / 2.0

    '''Aggregate placement judgment'''
    for iter in range(IterLimit):
        if iter < IterLimit - 1:
            # Generate random coordinates for the center
            x1 = np.random.uniform(0 + radius, ConcreteLength - radius)
            y1 = np.random.uniform(0 + radius, ConcreteWidth - radius)
            z1 = np.random.uniform(0 + radius, ConcreteHeight - radius)
            # Store the generated coordinates and radius of the aggregate
            Aggpoint = (x1, y1, z1, radius)
            # Check for intersection with existing aggregates
            if len(Aggpoints) == 0:
                Aggpoints.append(Aggpoint)
                AggData.append([i + 1, x1, y1, z1, radius])
                # Calculate the volume of the generated aggregate
                GenerateVolume = 4 * 3.14 * radius**3 / 3
                TotalAggVolume += GenerateVolume
                Process = TotalAggVolume * 100 / TargetVolume
                Friction = TotalAggVolume * 100 / ConcreteVolume
                i += 1
                print(f'Success! Current Volume = {TotalAggVolume:.3f}, Process = {Process:.3f}%, Friction = {Friction:.3f}% Current Step = {Aggnum}')
                print(f'Current Iteration Time = {iter:.3f}')
                break
            elif interact_judgement(Aggpoints, Aggpoint):
                Aggpoints.append(Aggpoint)
                AggData.append([i + 1, x1, y1, z1, radius])
                GenerateVolume = 4 * 3.14 * radius**3 / 3
                TotalAggVolume += GenerateVolume
                Process = TotalAggVolume * 100 / TargetVolume
                Friction = TotalAggVolume * 100 / ConcreteVolume
                i += 1
                print(f'Success! Current Volume = {TotalAggVolume:.3f}, Process = {Process:.3f}%, Friction = {Friction:.3f}% Current Step = {Aggnum}')
                print(f'Current Iteration Time = {iter:.3f}')
                break
        else:
            print('Fail! Generation Failure')
            break
    if (TotalAggVolume - TargetVolume) >= 0.001:
        break

'''Generate BC'''
mdb.Model(name='BiocharLibrary', modelType=STANDARD_EXPLICIT)
mdb.Model(name='BiocharConcrete', modelType=STANDARD_EXPLICIT)
myModel_Lib = mdb.Model(name='BiocharLibrary')
myModel_BC = mdb.Model(name='BiocharConcrete')
Modelname_BC = myModel_BC.name
Modelname_Lib = myModel_Lib.name

# Remove mask
session.journalOptions.setValues(replayGeometry=INDEX, recoverGeometry=INDEX)
session.journalOptions.setValues(replayGeometry=COORDINATE, recoverGeometry=COORDINATE)

'''Intersection judgment function'''
def interact_judgement(points, point):
    x1, y1, z1, r1 = point
    x_z, x_f = x1 + r1, x1 - r1
    y_z, y_f = y1 + r1, y1 - r1
    z_z, z_f = z1 + r1, z1 - r1

    if (x_z > BaseGapMax * ConcreteLength or x_f < BaseGapMin * ConcreteLength or
        y_z > BaseGapMax * ConcreteWidth or y_f < BaseGapMin * ConcreteWidth or
        z_z > BaseGapMax * ConcreteHeight or z_f < BaseGapMin * ConcreteHeight):
        return False

    for ii in points:
        x2, y2, z2, r2 = ii
        distance = math.sqrt((x1 - x2)**2 + (y1 - y2)**2 + (z1 - z2)**2)
        if distance < BioGap * (r1 + r2):
            return False
    return True

'''Define Biochar Library geometry'''
Biochar_Num = 6
BiocharCentroid = np.zeros((Biochar_Num, 4))  # Record Biochar ID and centroid coordinates xyz
MaximumDistance_Biochar = np.zeros((Biochar_Num, 1))  # Record the maximum distance of Biochar (based on Excel)
Length_Biochar = np.zeros((Biochar_Num, 1))  # Record the length of Biochar (based on Excel)
Volume_Biochar = np.zeros((Biochar_Num, 1))  # Record the volume of Biochar
LDration = np.zeros((Biochar_Num, 1))  # Record the length-to-diameter ratio of Biochar

'''Import geometry from Library'''
for i in range(Biochar_Num):
    iges = mdb.openIges('F:/Papers/Topic 1 Biochar as Fine Aggregate in mortar scale/Biochar Concrete Code/BiocharCode/BiocharGeom/Biochar-' + str(i + 1) + 'O.sat', msbo=False, 
                        trimCurve=DEFAULT, scaleFromFile=OFF)
    mdb.models['BiocharLibrary'].PartFromGeometryFile(name='Biochar-' + str(i + 1), 
                                                      geometryFile=iges, combine=False, stitchTolerance=1.0, 
                                                      dimensionality=THREE_D, type=DEFORMABLE_BODY, 
                                                      convertToAnalytical=1, stitchEdges=1)

# Record biochar geometric size information:
MaximumDistance_Biochar[0][0] = 3.8
MaximumDistance_Biochar[1][0] = 1.7
MaximumDistance_Biochar[2][0] = 2.9
MaximumDistance_Biochar[3][0] = 3.2
MaximumDistance_Biochar[4][0] = 2.8
MaximumDistance_Biochar[5][0] = 3.5
# Take five as example, more geometry can be found in Biochar Library Fiel

Length_Biochar[0][0] = 3.8
Length_Biochar[1][0] = 1.7
Length_Biochar[2][0] = 2.9
Length_Biochar[3][0] = 3.2
Length_Biochar[4][0] = 2.8
Length_Biochar[5][0] = 3.5

Volume_Biochar[0][0] = 2 * 3.1416 * MaximumDistance_Biochar[0][0] * Length_Biochar[0][0]**2 / 8
Volume_Biochar[1][0] = 2 * 3.1416 * MaximumDistance_Biochar[1][0] * Length_Biochar[1][0]**2 / 8
Volume_Biochar[2][0] = 2 * 3.1416 * MaximumDistance_Biochar[2][0] * Length_Biochar[2][0]**2 / 8
Volume_Biochar[3][0] = 2 * 3.1416 * MaximumDistance_Biochar[3][0] * Length_Biochar[3][0]**2 / 8
Volume_Biochar[4][0] = 2 * 3.1416 * MaximumDistance_Biochar[4][0] * Length_Biochar[4][0]**2 / 8
Volume_Biochar[5][0] = 2 * 3.1416 * MaximumDistance_Biochar[5][0] * Length_Biochar[5][0]**2 / 8

'''Define Biochar & Concrete percentage and placement rules'''
BiocharFraction = 0.09  # Biochar percentage
ConcreteLength = 50
ConcreteWidth = 50
ConcreteHeight = 50

TotalAggVolume = 0  # Volume before placement
ConcreteVolume = ConcreteHeight * ConcreteLength * ConcreteWidth
TargetVolume = ConcreteHeight * ConcreteLength * ConcreteWidth * BiocharFraction  # Target volume for Biochar

TryNumBiochar = 5000  # Maximum number of Biochar placement attempts
IterLimit = 10000  # Maximum iteration limit for placement

Biopoint = []
Biopoints = []
BioData = []
BioMesh = []

BioGap = 1.03  # Distance between aggregates (allowing intersection, so it can be set to less than 1)
BaseGapMax = 0.97  # Distance from the bottom
BaseGapMin = 0.03  # Distance from the top

'''Biochar placement'''
BioGenerate = 0
for BioNum in range(TryNumBiochar):
    BioLib = random.randint(1, Biochar_Num)
    mdb.models[Modelname_Lib].rootAssembly.Instance(dependent=ON, name='Biochar-' + str(BioNum), 
                                                    part=mdb.models[Modelname_Lib].parts['Biochar-' + str(BioLib)])
    for iter in range(IterLimit):
        if iter < IterLimit - 1:
            x1 = np.random.uniform(0 + MaximumDistance_Biochar[BioLib-1][0] / 2.0, ConcreteLength - MaximumDistance_Biochar[BioLib-1][0] / 2.0)
            y1 = np.random.uniform(0 + MaximumDistance_Biochar[BioLib-1][0] / 2.0, ConcreteWidth - MaximumDistance_Biochar[BioLib-1][0] / 2.0)
            z1 = np.random.uniform(0 + MaximumDistance_Biochar[BioLib-1][0] / 2.0, ConcreteHeight - MaximumDistance_Biochar[BioLib-1][0] / 2.0)

            Biopoint = (x1, y1, z1, MaximumDistance_Biochar[BioLib-1][0] / 2)

            if len(Biopoints) == 0:
                if interact_judgement(Aggpoints, Biopoint):
                    Biopoints.append(Biopoint)
                    TotalAggVolume += Volume_Biochar[BioLib-1][0]
                    Process = TotalAggVolume * 100 / TargetVolume
                    Friction = TotalAggVolume * 100 / ConcreteVolume
                    print(f'Success! Current Volume = {TotalAggVolume:.3f}, Process = {Process:.3f}%, Friction = {Friction:.3f}% Current Step = {BioNum}')
                    print(f'Current Iteration Time = {iter:.3f}')
                    # Random rotation around the centroid
                    axisPoint = (0, 0, 0)
                    xP, yP, zP = np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)
                    AxisDirection = (xP, yP, zP)
                    RotationAngle = random.randint(0, 360)

                    BioData.append([BioNum+1, x1, y1, z1, xP, yP, zP, RotationAngle, MaximumDistance_Biochar[BioLib-1][0] / 2, BioLib])
                    BioMesh.append([BioGenerate+1, x1, y1, z1, xP, yP, zP, RotationAngle, MaximumDistance_Biochar[BioLib-1][0] / 2, BioLib])

                    AxisDirection = (BioData[BioGenerate][4], BioData[BioGenerate][5], BioData[BioGenerate][6])

                    mdb.models[Modelname_Lib].rootAssembly.rotate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                  axisPoint=axisPoint, 
                                                                  axisDirection=AxisDirection, 
                                                                  angle=BioData[BioGenerate][7])

                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(BioData[BioGenerate][1], 0.0, 0.0))
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(0.0, BioData[BioGenerate][2], 0.0))    
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(0.0, 0.0, BioData[BioGenerate][3]))

                    BioGenerate += 1
                    break
            elif interact_judgement(Biopoints, Biopoint):
                if interact_judgement(Aggpoints, Biopoint):
                    Biopoints.append(Biopoint)
                    TotalAggVolume += Volume_Biochar[BioLib-1][0]
                    Process = TotalAggVolume * 100 / TargetVolume
                    Friction = TotalAggVolume * 100 / ConcreteVolume
                    print(f'Success! Current Volume = {TotalAggVolume:.3f}, Process = {Process:.3f}%, Friction = {Friction:.3f}% Current Step = {BioNum}')
                    print(f'Current Iteration Time = {iter:.3f}')

                    axisPoint = (0, 0, 0)
                    xP, yP, zP = np.random.uniform(-1, 1), np.random.uniform(-1, 1), np.random.uniform(-1, 1)
                    AxisDirection = (xP, yP, zP)
                    RotationAngle = random.randint(0, 360)

                    BioData.append([BioNum+1, x1, y1, z1, xP, yP, zP, RotationAngle, MaximumDistance_Biochar[BioLib-1][0] / 2, BioLib])
                    BioMesh.append([BioGenerate+1, x1, y1, z1, xP, yP, zP, RotationAngle, MaximumDistance_Biochar[BioLib-1][0] / 2, BioLib])

                    AxisDirection = (BioData[BioGenerate][4], BioData[BioGenerate][5], BioData[BioGenerate][6])

                    mdb.models[Modelname_Lib].rootAssembly.rotate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                  axisPoint=axisPoint, 
                                                                  axisDirection=AxisDirection, 
                                                                  angle=BioData[BioGenerate][7])

                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(BioData[BioGenerate][1], 0.0, 0.0))
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(0.0, BioData[BioGenerate][2], 0.0))    
                    mdb.models[Modelname_Lib].rootAssembly.translate(instanceList=('Biochar-' + str(BioData[BioGenerate][0]-1),),
                                                                     vector=(0.0, 0.0, BioData[BioGenerate][3]))

                    BioGenerate += 1
                    break
        else:
            del mdb.models[Modelname_Lib].rootAssembly.features['Biochar-' + str(BioNum)]
            print('Fail! Generation Failure')
            break

    if (TotalAggVolume - TargetVolume) >= 0.001:
        break



