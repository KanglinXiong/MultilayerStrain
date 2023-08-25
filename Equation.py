#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

''' 
    the equilibrium equation of force, moment, interface strain
    m.x == b
    x = [ f1, f2, ..., fn, r]
    where f is force per length, r is inverse of radius of curverture

            ---------------------
       <-- |      layer i+1      | --> external force, mismatch strain
            ---------------------
       <-- |      layer i        | -->
            ---------------------
                    /|\
                     |  R
    interface strain = lattice mismatch * (1 - relax) + thermal mismatch
    lattice mismatch is locked at initial temperature
    
    the model is only valid when R >> total film thickness

    reference: doi: 10.1063/1.323970
'''

import Misc

# balance of force
def forceEq(numOfLayers):
    row = [1.0]*(numOfLayers + 1)
    row[-1] = 0.0
    return [row, 0.0]


# balance of moments
# NOTE: the balanced (neutral) plane is not known, it can be found by minimizing strain energy
def momentEq(numOfLayers, youngList, thickList, neutralPlanePos = 0):
    # the fomula for area moment of inertia of multilayered film from the ref may be wrong
    # think about it: a stack of GaN layers equals to a thick layer of GaN 
    if(neutralPlanePos<0 or neutralPlanePos>sum(thickList)):
        raise Exception("Position of neutral plane should be within 0 and total thickness!")
    row = list(map(lambda i: sum(thickList[0:i])-thickList[i-1]/2.0 - neutralPlanePos,\
                   range(numOfLayers)))
    # by definition, use Integrate[(y)^2, {y, d-x0, d+h-x0}] == ((d+h-x0)^3 - (d-x0)^3)/3
    # x0 is the position of neutral plane
    coeff = sum(map(lambda i:\
        youngList[i]*((sum(thickList[0:i+1]) - neutralPlanePos)**3.0 -\
                      (sum(thickList[0:i])   - neutralPlanePos)**3.0), range(numOfLayers)))/3.0
    row.extend([coeff])
    return [row, 0.0]


# at interface, coherence is somehow mantained
# mismatchStrain + forceInducedStrainDiff == curvertureInducedStrainDiff
# NOTE: i in range(numOfLayers-1)
# double check the signs
def interfaceEq(numOfLayers, i, youngList, thickList, mismatchStrainList):
    row = [0.0]*(numOfLayers + 1)
    row[i] = -1.0/(youngList[i] * thickList[i])         # force strain of i
    row[i+1] = 1.0/(youngList[i+1] * thickList[i+1])    # force strain of i+1
    row[-1] =  -1.0/2.0*(thickList[i] + thickList[i+1]) # curveture strain diff
    return [row, (0.0 - mismatchStrainList[i])]         # mismatch strain diff

# m[rowNum].x = b[rowNum], return row of m and element of b
# only valid when R >> total thickness, otherwise, it turns into nonliear eq for R
def buildEq(numOfLayers, youngList, thickList, mismatchStrainList, neutralPlanePos = 0):
    buf = list(map(lambda i: interfaceEq(numOfLayers, i, youngList, thickList, mismatchStrainList),\
        range(numOfLayers-1)))
    buf.extend( [forceEq(numOfLayers)] )
    buf.extend( [momentEq(numOfLayers, youngList, thickList, neutralPlanePos)] )
    m = [[]]*(numOfLayers + 1)
    b = [0.0]*(numOfLayers +1)
    for i in range(numOfLayers+1):
        m[i] = buf[i][0]
        b[i] = buf[i][1]
    return [m, b]

if __name__ == "__main__":
    numOfLayers = 5
    youngList = [1]*numOfLayers
    thickList = [100]*numOfLayers
    mismatchStrainList = range(numOfLayers)
    [m, b] = buildEq(numOfLayers, youngList, thickList, mismatchStrainList)
    Misc.display(m)
    print(b)
    pass
