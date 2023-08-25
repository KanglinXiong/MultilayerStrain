#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

# AlxGa1-xN alloy, for 0 < x < 1

# name of the alloy
name = "AlGaN"

# two materials to interpolate, must be simple material as GaN, Si111
# otherwise, there will be recursive loop
boundary = [
    [1.0, "AlN"],
    [0.0, "GaN"]
]

# typical temperature for growth, minus value means not set
growthTemperature  = 1050
