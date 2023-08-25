#!/mingw64/bin/python

#########################################################
# This file is part of the MultilayerStrain package.    #
# Version 0.1.0                                         #
# Copyright (c) 2016 and later, Kanglin Xiong.          #
#########################################################

''' simply draw line charts '''


# the [0, 0] is at the center of the canvas
# direction: x is right, y is up

import turtle, math
import Color

def showClickXYonTitle(x, y):
    turtle.penup()
    turtle.goto(x, y)
    turtle.pendown()
    turtle.write("test")
    turtle.title("[x, y] = [" + str(x) + ", " + str(y) + "]" )
    turtle.penup()

def config():
    turtle.hideturtle() # hide turtle to speed up plot
    turtle.speed(10) # moving speed
    turtle.delay(0)  # animation speed

def transpose(xyList):
    return list(map(lambda xy: [xy[1], xy[0]], xyList))

def translate(xyList, tx, ty):
    return list(map(lambda xy: [xy[0]+tx, xy[1]+ty], xyList))

def scale(xyList, sx, sy):
    return list(map(lambda xy: [xy[0]*sx, xy[1]*sy], xyList))

def drawAxes(x1, x2, y1, y2):
    pass

def drawLine(xyList, label = "", colorString = "black", penWidth = 2):
    if(len(xyList) < 2):
        raise Exception("Cannot draw a line for less than 2 points!")
    # setup pen
    turtle.width (penWidth) 
    turtle.pencolor(colorString)
    # draw
    turtle.penup()
    [x, y] = xyList[0]
    turtle.goto(x, y)
    turtle.pendown()
    for [x, y] in xyList:
        turtle.goto(x, y)
        #turtle.dot()
    turtle.penup()
    turtle.forward(10)
    turtle.pendown()
    turtle.write(label, "center")
    turtle.penup()

def drawMultiLine(xyListArray, labelList, colorStringList):
    chartWidth = 600
    chartHeight = 500
    xmin = min(map(lambda xyList: min(map(lambda xy: xy[0], xyList)), xyListArray)) 
    xmax = max(map(lambda xyList: max(map(lambda xy: xy[0], xyList)), xyListArray)) 
    ymin = min(map(lambda xyList: min(map(lambda xy: xy[1], xyList)), xyListArray)) 
    ymax = max(map(lambda xyList: max(map(lambda xy: xy[1], xyList)), xyListArray)) 
    # translate
    tx = 0 - (xmin + xmax)/2
    ty = 0 - (ymin + ymax)/2
    # scale
    sx = 1
    if(xmax - xmin != 0.0): sx = chartWidth/(xmax - xmin)
    sy = 1
    if(ymax - ymin != 0.0): sy = chartWidth/(ymax - ymin)
    # process and plot, line by line
    for i in range(len(xyListArray)):
        xyList = scale(translate(syListArray[i], tx, ty), sx, sy)
        drawLine(xyList, labelList[i], colorStringList[i])
    return

def done():
    turtle.done()

def mainloop():
    turtle.mainloop()

if __name__ == "__main__":
    config()
    xyList = list(map(lambda x:[x, 100*math.sin(x/20)], range(100)))
    drawLine(transpose(xyList),"sin",  "blue")
    turtle.onclick(showClickXYonTitle)
    turtle.done()
