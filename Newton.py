#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

import random, math

''' find max or min value for a func using Newton's method '''

def grad(f, x):
    #NOTE: set dx depending on range of x
    dx = 1e-5*max(1, abs(x))
    return (f(x + dx/2) -f(x - dx/2))/dx

class Min:
    def __init__(self, func, x1, x2, tol = 1e-10):
        self.func = func
        self.diff = lambda x:grad(func, x)
        self.xmin = x1
        self.xmax = x2
        self.tolerance = tol

    def run(self, maxIter = 1e2):
        x = (self.xmin + self.xmax)/2.0
        count = 0
        print("iter, x, error, func(x)")
        while(count < maxIter):
            count += 1
            # root found
            if(abs(self.diff(x)) <= self.tolerance):
                break
            # update x
            if(grad(self.diff, x) == 0.0):
                x = random.uniform(self.xmin, self.xmax)
            else:
                dx = - self.diff(x)/grad(self.diff, x)
                if(x + dx < self.xmin or x + dx > self.xmax):
                    x = random.uniform(self.xmin, self.xmax)
                else:
                    x += dx
            # print process
            print(count, x, self.diff(x), self.func(x))
        if(abs(self.diff(x)) > self.tolerance):
            raise Exception("Newton failed to find local optima!")
        return [x, self.diff(x), self.func(x)]

Max = Min

if __name__ == "__main__":
    m = Min(lambda x:(x-1)**2 + 0.2, 0, 18)
    print(m.run())
    m = Min(lambda x:0 - math.cos(x), -2, 3)
    print(m.run())
    m = Max(lambda x: math.cos(x)**2, -2, 3)
    print(m.run())




