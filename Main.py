#!/mingw64/bin/python

#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################


import sys
import Misc, Elasticity

''' commom computing tasks '''


# given stack at T1, compute r at T2, cooldown or heatup
# lattice mismatch is locked at T1
# thermal mismatch is proportional to T2-T1 for all layers
# sweep T2

helpStr = "Usage: python Main.py <script> T1 [T2] [number of steps]."

def run():
    numOfArgs = len(sys.argv)
    if(numOfArgs < 3 or numOfArgs >5):
        print(helpStr)
        print("Error: number of arguments should be 3 to 5!")
        print("Input:" + " ".join(sys.argv))
        return
    script = sys.argv[1]
    if(not Misc.queryScript(script)):
        print("Cannot find script:", script)
        print("Accessible scripts:", ", ".join(Misc.listScripts()))
        return
    print("parsing", Misc.scriptFilename(script))
    structure = Elasticity.Structure(script)
    print("running")
    if(numOfArgs == 3):
        t1 = float(sys.argv[2])
        result = structure.statusquo(t1)
    if(numOfArgs == 4):
        t1 = float(sys.argv[2])
        t2 = float(sys.argv[3])
        result = structure.rampTemperature(t1, t2)
    if(numOfArgs == 5):
        t1 = float(sys.argv[2])
        t2 = float(sys.argv[3])
        numOfSteps = int(sys.argv[4])
        result = structure.rampTemperature(t1, t2, numOfSteps)
    outName = "_".join(sys.argv[1:])
    print("saving", Misc.outputFilename(outName))
    Misc.saveResult(result, outName)
    print("done")

if __name__ == "__main__":
    run()


