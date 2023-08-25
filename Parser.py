#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

import re
import Misc, Unit

class Parser:
    ''' parse script file into a list of list of material name, thickness, relax ratio '''

    def __init__(self):
        self.layerInfoList = []

    def parseThicknessWithUnit(self, thicknessString):
        ''' string => float '''
        #NOTE: no spacing allowed between number and units, 20nm, 1um, 0.1m
        pattern = "((" + ")|(".join(Unit.length.keys()) + "))"
        rltLst = re.split(pattern, thicknessString.strip())
        thickness = float(rltLst[0]) # may fail if invalid symbols exist
        if len(rltLst) > 1 :
            thickness *= Unit.length[rltLst[1]]
        else:
            thickness *= Unit.length["default"]
        return thickness

    @staticmethod
    def discretizeGradedLayer(cmd = "[Al50%05%GaN 120.0 10]"):
        ''' [Ax1%x2BC% thickness numberOflayers] => string for ordinary layers '''
        eleList = cmd.strip("[]").split()
        if(len(eleList) != 3):
            raise Exception("Invalid graded layer " + cmd)
        [name, d, n] = [eleList[0], float(eleList[1]), int(eleList[2])]
        if n <= 2:
            raise Exception("Invalid graded layer " + cmd)
        [x1, x2] = list(map(int, re.findall("[0-9]{2}", name)))
        dx = (x2 - x1)/(n - 1.0)
        xStrList = list(map(lambda i: str(round(x1 + dx*i)).zfill(2) + "%", range(n)))
        labelList = re.split("[0-9]{2}%[0-9]{2}%", name)
        strList = list(map(lambda i: labelList[0] + xStrList[i] + labelList[1], range(n)))
        string = ""
        for i in range(n - 1):
            string = string + strList[i] + " " + str(d/float(n)) + ","
        string = string + strList[-1] + " " + str(d/float(n))
        return(string)

    def parseStringLine(self, line):
        ''' parse a line string '''
        layers = line.partition("{")[2].partition("}")[0]
        repeat = line.partition("*")[2].strip()
        # discretize graded layers like [Al25%05%GaN 120.0 11] if any
        gradedLayerList = re.findall("\[.*?\]", layers)
        for gradedLayer in gradedLayerList:
            layers = layers.replace(gradedLayer, self.discretizeGradedLayer(gradedLayer))
        # parse all the layers, get [[mat, d, relaxIfProvided], ...]
        tmpStack = []
        if(len(repeat)):
            repeat = int(repeat)
        else:
            repeat = 1
        for layer in layers.split(","):
            layer = layer.strip().split(" ")
            layer[1] = self.parseThicknessWithUnit(layer[1])
            # the 3rd param of a layer is lattice relaxation ratio
            if(len(layer) == 2):
                layer.extend([0.0]) # NOTE: default relaxation ratio is 0.0
            elif(len(layer) == 3):
                layer[2] = float(layer[2])
            tmpStack.append(layer)
        for layer in (tmpStack * repeat):
            self.layerInfoList.append(layer)
        pass

    # script syntax is strict
    def run(self, scriptName):
        ''' read script, set layerInfoList as [[matName, d, r], etc] '''
        self.layerInfoList = []
        fileObj = open(Misc.scriptFilename(scriptName), mode="rt", newline='\n')
        for line in fileObj:
            line = line.strip()
            # line with '#' is comment
            if((not len(line)) or ("#" in line)): continue
            line = line.expandtabs(1)
            while(" "*2 in line):
                line = line.replace(" "*2, " ")
            self.parseStringLine(line)
        return self.layerInfoList

    def getResult(self):
        return self.layerInfoList

if __name__ == "__main__":
    parser = Parser()
    gradedLayer = "[Al50%05%GaN 120.0 10]"
    print(gradedLayer, "is expanded as")
    print(parser.discretizeGradedLayer("[Al50%05%GaN 120.0 10]"))
    scriptName = "testParser"
    print("parsing script: ", scriptName)
    print("materialName, thickness, relaxRatio")
    Misc.display(parser.run(scriptName))

