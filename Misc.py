#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


import os, os.path, re
import Unit


#########################################################
# manage the locations of files

def listMaterialModuleFiles():
    directory = os.path.dirname(os.path.realpath(__file__))
    directory = os.path.join(directory, "material")
    return os.listdir(directory)

def queryMaterialModule(materialName):
    if(materialName.strip() + ".py" in listMaterialModuleFiles()):
        return True
    return False

def materialModuleName(materialName):
    return("material." + materialName)

def listScripts():
    directory = os.path.dirname(os.path.realpath(__file__))
    directory = os.path.join(directory, "script")
    return list(map(lambda s: (s.split(".txt"))[0], os.listdir(directory)))

def queryScript(script):
    if(script.strip() in listScripts()):
        return True
    return False

def scriptFilename(scriptName, suffix = ".txt"):
    directory = os.path.dirname(os.path.realpath(__file__))
    return(os.path.join(directory, "script", scriptName + suffix))

def outputFilename(dataName, suffix = "_rlt.csv"):
    directory = os.path.dirname(os.path.realpath(__file__))
    return(os.path.join(directory, "output", dataName + suffix))

#########################################################
# save [[temperture, radius, error, stressData, strainData]] to file

def pad2dArray(array, targetLen, pad =" "):
    if(len(array) >= targetLen): return array
    padEle = [pad]*len(array[0])
    array.extend([padEle]*(targetLen - len(array)))
    return array

def pad1dList(array, targetLen, pad ="-"):
    if(len(array) >= targetLen): return array
    array.extend([pad]*(targetLen - len(array)))
    return array

def dataRow2Str(dataList, rowIdx):
    array = []
    for data in dataList:
        array.extend(data[rowIdx])
    return ",".join(map(str, array)) + "\n"

def saveResult(result, scriptName):
    ''' save result of Structure.rampTemperature() into CSV '''
    # prepare data
    rltLen = len(result)
    stressLen = max(map(lambda i: len(result[i][3]), range(rltLen)))
    strainLen = max(map(lambda i: len(result[i][4]), range(rltLen)))
    maxLen = max(rltLen, stressLen, strainLen)
    dataList = []
    tempRadErrList = []
    titleList = [[["T(K)", "R(m)", "neutralPlanePos(um)"]]]
    for rlt in result: 
        tempRadErrList.append([rlt[0], rlt[1]/Unit.length["m"], rlt[2]/Unit.length["um"]])
        for i in range(len(rlt[3])):
            rlt[3][i][0] /= Unit.length["um"]
            rlt[3][i][1] /= Unit.GPa
        dataList.append(pad2dArray(rlt[3], maxLen))
        titleList.append([pad1dList(["x(um)", "stress(GPa)@"+str(rlt[0])], len(rlt[3][0]))])
        for i in range(len(rlt[4])):
            rlt[4][i][0] /= Unit.length["um"]
        dataList.append(pad2dArray(rlt[4], maxLen))
        titleList.append([pad1dList(["x(um)", "strain@"+str(rlt[0])], len(rlt[4][0]))])
    tempRadErrList = pad2dArray(tempRadErrList, maxLen) 
    dataList.insert(0, tempRadErrList) 
    # write title and data
    filename = outputFilename(scriptName)
    file = open(filename, "w")
    file.write(dataRow2Str(titleList, 0))
    for i in range(maxLen):
        file.write(dataRow2Str(dataList, i))
    file.close()

def LoadCSV(scriptName):
    pass

#########################################################
# print in a better way

def display(data):
    if type(data) == type([]):
        for ele in data: print(ele)
    else:
        print(data)
    pass

#########################################################
# given [[x, y], etc], interpolate to get y0 for x0

def linearInterpolate(xyLists, x):
    if(len(xyLists) < 1):
        raise Exception("Invalid data", xyLists)
    # two boundaries
    if(x <= xyLists[0][0]): # lower
        return(xyLists[0][1])
    elif(x >= xyLists[-1][0]): # upper
        return(xyLists[-1][1])
    # find the gap x resides in
    i = 0
    while(i < len(xyLists) and x > xyLists[i][0]):
        i += 1
    # interpolate
    [d1, d2] = [x - xyLists[i-1][0], xyLists[i][0]- x]
    return((xyLists[i-1][1]*d2 + xyLists[i][1]*d1)/(d1+d2))

# Integrate over range x0, x1 with linear interpolate
def linearIntegrate(xyLists, x0, x1):
    numOfPts = max(1000, len(xyLists))
    xList = list(map(lambda i: x0 + i*(x1-x0)/(numOfPts-1), range(numOfPts)))
    yList = list(map(lambda x: linearInterpolate(xyLists, x), xList))
    return sum(map(\
        lambda i: (xList[i+1]-xList[i])*(yList[i]+yList[i+1])/2.0,\
        range(numOfPts-1)))
  
#########################################################

if __name__ == "__main__":
    xyLists =  [[100, 100]]
    xList = range(-100, 200, 20)
    print("xy list", xyLists, "x list", xList)
    print(list(map(lambda i: linearInterpolate(xyLists,i), xList)))
    print(linearIntegrate(xyLists, 0, 100))
    print(linearIntegrate(xyLists, 100, 0))
    xyLists =  [[0, 0], [100, 100]]
    print("xy list", xyLists, "x list", xList)
    print(list(map(lambda i: linearInterpolate(xyLists,i), xList)))
    print(linearIntegrate(xyLists, 0, 100))
    print(linearIntegrate(xyLists, 100, 0))
    print(linearIntegrate(xyLists, 100, 100))
    print(pad2dArray([[1]]*3, 10, 0))
    pass

