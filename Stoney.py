#!/mingw64/bin/python

#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

''' Stoney's equation '''

import Misc, Unit, Elasticity

class Stoney:
    def __init__(self, waferMaterial, waferThicknessInMicro):
        self.material = Elasticity.Material(waferMaterial)
        self.thickness = waferThicknessInMicro * Unit.length["um"]

    def setTemperature(self, temp):
        self.material.setTemperature(temp)

    def getFilmStress(self, filmThicknessInNano, radiusInMeter, InitRadiusInMeter = None):
        if(float(filmThicknessInNano) == 0.0):
            raise Exception("Thickness of film cannot be zero!")
        if(float(radiusInMeter) == 0.0):
            raise Exception("Radius of curvature cannot be zero!")
        filmThickness = filmThicknessInNano * Unit.length["nm"]
        radius = radiusInMeter * Unit.length["m"]
        # if initially the wafer is curved
        if(isinstance(InitRadiusInMeter, int) or isinstance(InitRadiusInMeter, float)):
            initRadius = InitRadiusInMeter * Unit.length["m"]
        else:
            initRadius = float('inf')* Unit.length["m"]
        young = self.material.getYoungsModulus()
        poisson = self.material.getPoissonsRatio()
        stress = young/(1.0 - poisson) * (self.thickness**2)/(6 * filmThickness)*(1/radius - 1/initRadius)
        return stress / Unit.GPa

    def getRadiusOfCurvature(self, filmThicknessInNano, stressInGPa, InitRadiusInMeter = None):
        if(float(filmThicknessInNano) == 0.0):
            raise Exception("Thickness of film cannot be zero!")
        if(float(stressInGPa) == 0.0):
            raise Exception("stress cannot be zero!")
        filmThickness = filmThicknessInNano * Unit.length["nm"]
        stress = stressInGPa * Unit.GPa
        # if initially the wafer is curved
        if(isinstance(InitRadiusInMeter, int) or isinstance(InitRadiusInMeter, float)):
            initRadius = InitRadiusInMeter * Unit.length["m"]
        else:
            initRadius = float('inf')* Unit.length["m"]
        young = self.material.getYoungsModulus()
        poisson = self.material.getPoissonsRatio()
        effRadius = young/(1.0 - poisson) * (self.thickness**2)/(6 * stress * filmThickness)
        radius  = 1/(1/effRadius - 1/initRadius)
        return radius / Unit.length["m"]
      


if __name__ == "__main__":
    stoney = Stoney("Si111", 500)
    thickness = 2000
    print("please input radius in meter:")
    r = float(input())
    stress = stoney.getFilmStress(thickness, r)
    print("stress of film of", thickness, "nm:", stress,"GPa")
    print("please input stress in GPa:")
    stress = float(input())
    print("radius of curvature:", stoney.getRadiusOfCurvature(thickness, stress))

