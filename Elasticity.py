#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

''' Read in script, build up structure
    Get parameters for equations to calculate radius and strain under in biaxial strain 
    lattice mismatch is frozen at intial temperature 
    thermal mismatch is caused by temperature deviation from intial temperature
'''

import re, copy
import Misc, Parser, Unit, Equation, Solver, Newton

class Material:
    ''' material elastic parameters '''
    def __init__(self, name):
        self.name = name.strip()
        self.dataModule = None
        self.interpolationFlag = False
        self.boundary = []
        self.importDataModule()
        self.temperature = self.roomTemperature()

    @staticmethod
    def roomTemperature():
        return 300

    @staticmethod
    def getXofAxBC(name = "A25%BC"):
        digits = re.findall("[0-9]{2}", name)
        if len(digits) != 1:
            return -1.0
        return float(digits[0])/100.0

    @staticmethod
    def getABCofAxBC(name = "A25%BC"):
        return re.sub("([0-9]{2}\%)", "", name)

    @staticmethod
    def interpolateSect(x1, y1, x2, y2, x):
        return ((x - x1)*y2 + (x2 - x)*y1)/(x2-x1)

    def setBoundary(self):
        ''' refer to simple materials, otherwise recursive issue '''
        self.boundary = copy.deepcopy(self.dataModule.boundary)
        self.boundary[0][1] = Material(self.dataModule.boundary[0][1])
        self.boundary[1][1] = Material(self.dataModule.boundary[1][1])

    def importDataModule(self):
        ''' parameter file is imported as python module '''
        moduleName = self.getABCofAxBC(self.name)
        if not Misc.queryMaterialModule(moduleName):
            raise Exception("Material " + self.name + " is not supported!")
        fullModuleName = Misc.materialModuleName(moduleName)
        localName = moduleName + "DataModule"
        exec("import " + fullModuleName + " as " + localName)
        self.dataModule =locals()[localName]
        if("boundary" in dir(self.dataModule)):
            self.interpolationFlag = True
            self.setBoundary()

    def setTemperature(self, temperature):
        self.temperature = temperature
        for bd in self.boundary:
            bd[1].setTemperature(temperature)

    def getLattice(self, temp = None):
        t = self.temperature
        if isinstance(temp, float) or isinstance(temp,int): t = temp
        if self.interpolationFlag:
            return self.interpolateSect(\
                self.boundary[0][0], self.boundary[0][1].getLattice(t),\
                self.boundary[1][0], self.boundary[1][1].getLattice(t),\
                self.getXofAxBC(self.name))
        expansionRatio = Misc.linearIntegrate(\
            self.dataModule.thermalExpansionCoefficient, self.roomTemperature(), t) 
        return self.dataModule.lattice300K * (1.0 + expansionRatio)

    def getLattice300K(self):
        return self.getLattice(300)

    def getThermalExpansion(self, tempBegin, tempEnd):
        ''' temperature changes from begin to end '''
        latticeBegin = self.getLattice(tempBegin)
        latticeEnd = self.getLattice(tempEnd)
        expansion = (latticeEnd - latticeBegin)/latticeBegin - 1.0
        return expansion

    def getYoungsModulus(self, temp = None):
        t = self.temperature
        if isinstance(temp, float) or isinstance(temp,int): t = temp
        if self.interpolationFlag:
            return self.interpolateSect(\
                self.boundary[0][0], self.boundary[0][1].getYoungsModulus(t),\
                self.boundary[1][0], self.boundary[1][1].getYoungsModulus(t),\
                self.getXofAxBC(self.name))
        return Misc.linearInterpolate(self.dataModule.youngsModulus, t)*Unit.GPa

    def getPoissonsRatio(self, temp = None):
        t = self.temperature
        if isinstance(temp, float) or isinstance(temp,int): t = temp
        if self.interpolationFlag:
            return self.interpolateSect(\
                self.boundary[0][0], self.boundary[0][1].getPoissonsRatio(t),\
                self.boundary[1][0], self.boundary[1][1].getPoissonsRatio(t),\
                self.getXofAxBC(self.name))
        return Misc.linearInterpolate(self.dataModule.poissonsRatio, t)

    def getGrowthTemperature(self):
        if self.interpolationFlag:
            return self.interpolateSect(\
                self.boundary[0][0], self.boundary[0][1].getGrowthTemperature(),\
                self.boundary[1][0], self.boundary[1][1].getGrowthTemperature(),\
                self.getXofAxBC(self.name))
        return self.dataModule.growthTemperature


class Layer:
    ''' layer properties '''
    def __init__(self, material, thickness = 0.0, relaxationRatio = 0.0):
        self.material = material
        self.thickness = thickness
        self.bottomInterfaceRelax = relaxationRatio
        self.growthTemperature = material.getGrowthTemperature()
        self.name = self.material.name +"("+ str(self.thickness/Unit.length["um"])+")"
        self.force = None
        self.reciprocalOfRadius = None

    def setBottomInterfaceRelax(self, relaxationRatio):
        self.bottomInterfaceRelax = relaxationRatio

    def setGrowthTemperature(self, temp):
        # layer overwrite material
        self.growthTemperature = temp

    def resetGrowthTemperature(self):
        self.growthTemperature = material.getGrowthTemperature()

    def copy(self):
        # do not use copy.deepcopy to a Layer, use this method
        newLayer = Layer(self.material, self.thickness, self.bottomInterfaceRelax)
        newLayer.growthTemperature = self.growthTemperature
        newLayer.force = self.force
        newLayer.reciprocalOfRadius = self.reciprocalOfRadius
        return newLayer

    def setForceAndReciprocalOfRadius(self, f, r):
        ''' must called before getting stress or strain info '''
        self.force = f
        self.reciprocalOfRadius = r

    def getStress(self, x):
        ''' within a layer, biaxial stress(x) = force/thick + young*(x - thick/2)/(2*r) '''
        if(self.force == None or self.reciprocalOfRadius == None):
            raise Exception("Layer not set for stress and strain!")
        if(x<0 or x>self.thickness*(1+1e-8)):
            raise Exception("Position is outside of a layer!")
        biaxialYoung = self.material.getYoungsModulus()/(1.0 - self.material.getPoissonsRatio())
        return self.force/self.thickness +\
            biaxialYoung * (x - self.thickness/2.0)*self.reciprocalOfRadius/2.0

    def getStrain(self, x):
        ''' inplane biaxial strain by definition '''
        biaxialYoung = self.material.getYoungsModulus()/(1.0 - self.material.getPoissonsRatio())
        return self.getStress(x)/biaxialYoung 

    def getStrainEnergy(self):
        ''' total strain energy integrated along thickness, unit N*m/m^2 '''
        # Integrate[strain[x]*(strain[x]*young)/2, {x, 0, thick}]*2
        biaxialYoung = self.material.getYoungsModulus()/(1.0 - self.material.getPoissonsRatio())
        energy = self.force**2 / (self.thickness * biaxialYoung) + \
                 self.thickness**3 * biaxialYoung / 48 * self.reciprocalOfRadius**2
        return energy

class Structure:
    def __init__(self, script):
        self.layerInfoList = []
        self.matDict = None
        self.layerStack = []
        self.script = script
        self.buildStruct()

    def buildStruct(self):
        ''' The struct is a list of Layer instances as [layer1, layer2, etc]  '''
        parser = Parser.Parser()
        self.layerInfoList = parser.run(self.script)
        self.layerInfoList.reverse()
        # unique material names
        materialNameList = []
        for layer in self.layerInfoList:
            if(not layer[0] in materialNameList):
                materialNameList.append(layer[0])
        # unique material instances
        materialDict = []
        for name in materialNameList:
            materialDict.append([name, Material(name)]) # instantialize Material
        self.matDict = dict(materialDict)
        # create struct as Layer instances from self.layerInfoList
        self.layerStack = []
        # instantialize Layers
        for layer in self.layerInfoList:
            self.layerStack.append(Layer(self.matDict[layer[0]], layer[1], layer[2]))
        return self.layerInfoList

    def sampling(self, func = (lambda i,x:[i,x]), maxNumOfStepsPerLayer = 10, minThickStep = 10.0):
        numOfLayers = len(self.layerStack)
        [globalPos, localPos] = [0.0, 0.0]
        posValuePairList = []
        for i in range(numOfLayers):
            num = min(maxNumOfStepsPerLayer, int(self.layerStack[i].thickness/minThickStep)+1)
            if num<3: num = 3
            step = self.layerStack[i].thickness/(num - 1.0)
            localPos = 0.0
            posValuePairList.append([globalPos, func(i,localPos), self.layerStack[i].name])
            for j in range(1, num):
                [globalPos, localPos] = [globalPos+step, localPos+step]
                posValuePairList.append([globalPos, func(i,localPos), self.layerStack[i].name])
        return posValuePairList

    def stress(self):
        func = lambda i, x: self.layerStack[i].getStress(x)
        return self.sampling(func)

    def strain(self):
        func = lambda i, x: self.layerStack[i].getStrain(x)
        return self.sampling(func)

    @staticmethod
    def latticeMismatchStrain(botLayer, topLayer, temp = Material.roomTemperature()):
        ''' top - bot, at given temperature '''
        # assuming balanced lattice is a0, bot lattice a1, top lattice a2
        # (a2 - a0)/a0 - (a1 - a0)/a0 = (a2 - a1)/a0 ~ (a2 - a1)/a2 ~ (a2 - a1)/a1
        # dilemma: AlN/GaN is different from GaN/AlN
        # to resolve the issue, use (a1 + a2)/2 as denominator
        botLattice = botLayer.material.getLattice(temp)
        topLattice = topLayer.material.getLattice(temp)
        denominator = (topLattice + botLattice)/2.0 
        mismatch = (topLattice - botLattice)/denominator * (1.0 - topLayer.bottomInterfaceRelax)
        return mismatch

    @staticmethod
    def thermalMismatchStrain(botLayer, topLayer, tempBegin, tempEnd):
        ''' top - bot, at given temperatures, their difference matters '''
        # thermal mismatch is difference in expansion ratio
        botExpansion = botLayer.material.getThermalExpansion(tempBegin, tempEnd)
        topExpansion = topLayer.material.getThermalExpansion(tempBegin, tempEnd)
        mismatch = topExpansion - botExpansion
        return mismatch

    def getEqParameters(self, tempBegin, tempEnd):
        numOfLayers = len(self.layerStack) 
        # get the values at tempEnd
        youngList = list(map(lambda i: self.layerStack[i].material.getYoungsModulus(tempEnd),\
            range(numOfLayers)))
        poissonList = list(map(lambda i: self.layerStack[i].material.getPoissonsRatio(tempEnd),\
            range(numOfLayers)))
        # biaxial modulus E' = E/(1-v), where E is uniaxial Young's modulus, v is Poisson's ratio
        for i in range(numOfLayers): youngList[i] *= 1.0/(1 - poissonList[i]) 
        # thickness of each layer
        thickList = list(map(lambda i: self.layerStack[i].thickness, range(numOfLayers)))
        # mismatch strain = lattice mismatch * (1 - relax) + thermal mismatch 
        # due to relax, thermal and lattice mismatch need to consider independently
        mismatchStrainList = [0.0]*numOfLayers
        print("\ninterface, bottom, top, lattice Mismatch, thermal mismatch")
        for i in range(numOfLayers - 1):
            lm = self.latticeMismatchStrain(self.layerStack[i], self.layerStack[i+1], tempBegin)
            tm = self.thermalMismatchStrain(self.layerStack[i], self.layerStack[i+1], tempBegin, tempEnd)
            mismatchStrainList[i] = lm + tm # the i-th interface
            print(i, self.layerStack[i].name, self.layerStack[i+1].name, lm, tm)
        return [numOfLayers, youngList, thickList, mismatchStrainList]

    def run(self, eqParams):
        ''' build eq, solve eq, set stack, obtain stress '''
        # try to find neutral plane
        def strainEnergyFunc(neutralPlanePos):
            [m, b] = Equation.buildEq(eqParams[0], eqParams[1], eqParams[2], eqParams[3],\
                                      neutralPlanePos)
            eq = Solver.LinearEq(m, b)
            eq.solve()
            root = eq.getRoot()
            radius = float('inf')
            if(root[-1] != 0.0):  radius = 1.0/root[-1]
            energy = 0.0
            for i in range(len(self.layerStack)):
                self.layerStack[i].setForceAndReciprocalOfRadius(root[i], 1.0/radius)
                energy += self.layerStack[i].getStrainEnergy()
            return energy
        minimizer = Newton.Min(strainEnergyFunc, 0, sum(eqParams[2]), 1e-18)
        print("\ncomputing netrual plane")
        optima = minimizer.run()
        print("neutral plane found")
        print("position, error, energy")
        print(optima, "\n")
        neutralPlanePos = optima[0]
        # solve the equation
        [m, b] = Equation.buildEq(eqParams[0], eqParams[1], eqParams[2], eqParams[3],\
                                  neutralPlanePos) 
        eq = Solver.LinearEq(m, b)
        eq.solve()
        root = eq.getRoot()
        error = eq.error()
        # save root into layers
        radius = float('inf')
        if(root[-1] != 0.0):  radius = 1.0/root[-1]
        for i in range(len(self.layerStack)):
            self.layerStack[i].setForceAndReciprocalOfRadius(root[i], 1.0/radius)
        stressDist = self.stress()
        strainDist = self.strain()
        return [radius, neutralPlanePos, stressDist, strainDist]

    def rampTemperature(self, tempBegin, tempEnd, numOfTempSteps = 10):
        ''' generate a series of parameter set for m.x == b to cooldown or heatup '''
        # ramp to tempEnd from tempBegin, update thermal mismatch at each step
        if(numOfTempSteps<2):
            raise Exception("Error, number of temperture steps less than 2!")
        tempStep = (tempEnd - tempBegin)/(numOfTempSteps - 1.0)
        resultList = [None]*numOfTempSteps
        for i in range(numOfTempSteps):
            currentTemp = tempBegin + float(i)*tempStep
            print("current temperature", currentTemp)
            eqParameters = self.getEqParameters(tempBegin, currentTemp)
            rlt = self.run(eqParameters)
            resultList[i] = [currentTemp]
            resultList[i].extend(rlt)
        #the result is [[temperature, radius, neutralPlanePos, stressDist, strainDist], ...]
        return resultList 

    def statusquo(self, temp):
        eqParameters = self.getEqParameters(temp, temp)
        rlt = self.run(eqParameters)
        resultList = [[]]
        resultList[0] = [temp]
        resultList[0].extend(rlt)
        return resultList 

    def cooldown(self, numOfTempSteps = 100):
        ''' from growth temperature of last layer to RT '''
        return self.rampTemperature(self.layerStack[-1].growthTemperature,\
            Material.roomTemperature(), numOfTempSteps)

    def heatup(self, numOfTempSteps = 100):
        ''' from RT to growth temperature of last layer '''
        return self.rampTemperature(Material.roomTemperature(),\
            self.layerStack[-1].growthTemperature, numOfTempSteps)


if __name__ == "__main__":
    m = Material("GaN")
    print(m.name, "data has been imported as module")
    print(dir(m.dataModule))
    print("lattice at 300, 600, 900K: ", m.getLattice(300), m.getLattice(600), m.getLattice(900))
    m = Material("Al50%GaN")
    print(m.name, "data has been imported as module")
    print(dir(m.dataModule))
    print("lattice at 300, 600, 900K: ", m.getLattice(300), m.getLattice(600), m.getLattice(900))
    print("thermal expansion 900 to 300: ", m.getThermalExpansion(900, 300))
    print(m.boundary)
    print("growth T: ", m.getGrowthTemperature())
    s = Structure("testParser")
    Misc.display(s.layerInfoList)
    Misc.display(s.layerStack)
    print(s.matDict.keys())
    Misc.display(s.getEqParameters(300, 300))
    Misc.display(s.sampling())
    rlt = s.statusquo(300)
    print([rlt[0][0], rlt[0][1]])
    Misc.display(rlt[0][-1])
    Misc.saveResult(rlt, s.script)



