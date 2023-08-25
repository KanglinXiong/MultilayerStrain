#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

# in-plane parameter file of material use the python syntax

# name of the material
name = "GaN"

# lattice constant at room temperature, Angstrom
lattice300K = 3.189

# thermal expansion coefficient at different temperatures, [[T0,v0], [T1,v1], etc]
thermalExpansionCoefficient = [
    [300, 5.59e-6]
]

# Young's modulus at different temperatures, [[T0,v0], [T1,v1], etc], GPa
youngsModulus = [
    [300, 320]
]

# Poisson's ratio at different temperatures, [[T0,v0], [T1,v1], etc]
poissonsRatio = [
    [300, 0.25]
]

# spontaneous polarization in x, y, z direction, C/cm2
spontaneousPolarization = [0.0, 0.0, 0.0]

# piezoelectricity under biaxial-strain uxx = uyy, pz = coeff * uxx, C/cm2
piezoElectricStrainCoefficient= 0.0

# typical temperature for growth, minus value means not set
growthTemperature  = 850
