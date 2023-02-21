# -*- coding: utf-8 -*-
"""
Created on Mon Feb 20 11:46:50 2023

@author: howling
"""

#from https://www.femm.info/wiki/HalbachExample
#Modified Halbach array by Jose B. Almeida, November 2017

#Feedback from E. Tsesmelis
# Assume a wiggler system based on permanent magnets. 
# Permanent magnets are typically used for synchrotron light sources as they allow for the possibility of smaller gaps 
# and shorter periods compared to coil dominated electromagnets.
# Referring to the lecture in ‘Undulators and Wigglers’ from last term and the 
# MAD-X file from Exercise 1, start modelling in FEMM a normal-conducting magnetic
#  wiggler to be used for the pDR. Assume each wiggler has a length of 6.64 m and 
#  a field of 1.66 T. Further reading: Chapter 8 (Wigglers and Undulators) in K. Wille, 
#  The Physics of Particle Accelerators – An Introduction, Oxford University Press.

#currently using 50mm diameter beampipe

import femm
import matplotlib.pyplot as plt
import numpy as np

# The package must be initialized with the openfemm command.
femm.openfemm()

# modify these values to suit your needs
units = 'millimeters'	# should be one of "inches","millimeters","centimeters","mils","meters","micrometers" 
forceUnits = 'Newtons'		# should be one of "lbf", "Newtons", "kgf"
magnetWidth = 200 #vertically
magnetLength = 50.00191431900077 # mm horizontally
Ldrift=30 #mm
magnetSep = 25  # one half of separation between magnet layers vertically
magnetType= 'N52' #'SmCo32'  # has to be the EXACT name of a material in the materials library
steelThick = 50  # thickness of magnetic circuit steel
interMagnet = 0  # thickness of air between magnets
steelType = '1020 Steel'  # has to be the EXACT name of a material in the materials library
numrepeats= 41#41 gives pDR spec
mess_size=0 #set to zero for automatic
margin=1
boundhoriz=700
oddmagnetLength = Ldrift-2*interMagnet
print('oddmagnetlength is', oddmagnetLength, 'mm')
#left_edge=-1.5*magnetLength-2*interMagnet-margin
og_thickness=2*magnetLength+2*oddmagnetLength+4*interMagnet
lengthtotal=og_thickness*(numrepeats+1)
print('Total length is', lengthtotal/1000, 'm')
#-----------------------------------------------------------------------------------------------
#Wiggler details from MADX file 

#//Wiggler

# ANGLEB2 = 0.01750089335228726
# LB2 = 0.05000191431900077
# ANGLEB3 = -0.01750089335228726
# LB3 = 0.05000191431900077
# ANGLEB4 = -0.01750089335228726
# LB4 = 0.05000191431900077
# ANGLEB5 = 0.01750089335228726
# LB5 = 0.05000191431900077


#wdrift : drift, L := 0.03;


# wstart: marker;
# wend: marker;

# B2 : RBEND, ANGLE:= ANGLEB2, E1:= -.5 * ANGLEB2, E2:= .5 * ANGLEB2, L:= LB2;
# B3 : RBEND, ANGLE:= ANGLEB3, E1:= .5 * ANGLEB3, E2:= -.5 * ANGLEB3, L:= LB3;
# B4 : RBEND, ANGLE:= ANGLEB4, E1:= -.5 * ANGLEB4, E2:= .5 * ANGLEB4, L:= LB4;
# B5 : RBEND, ANGLE:= ANGLEB5, E1:= .5 * ANGLEB5, E2:= -.5 * ANGLEB5, L:= LB5;

# //WIGGLER: LINE=(wstart, 21*(B2,driftw,B3,driftw,B4,driftw,B5,driftw), wend);

# WIGGLER_CELL: LINE=(B2,WDRIFT,B3,WDRIFT,B4,WDRIFT,B5,WDRIFT);

# WIGGLER: LINE=(21*WIGGLER_CELL);

#-----------------------------------------------------------------------------------------------

# We need to create a new Magnetostatics document to work on.
femm.newdocument(0)

#addpath('C:\femm42\mfiles')

femm.mi_hidegrid()


## define the problem type, equivalent to the top level menu    Problem
#Function attributes: mi_probdef(freq,units,type,precision,depth,minangle,(acsolver))
femm.mi_probdef(0, units, 'planar', 1E-8, magnetWidth)

## adds these materials from the Material Library to the project
femm.mi_getmaterial(magnetType)	
femm.mi_getmaterial(steelType)
femm.mi_getmaterial('Air')


## this bit is a bit of a pain. Although we have set the problem units above,
## that has no effect on the calculation of c0 for the Asymptotic Boundary Condition
## We therefore need to create a scaling factor for c0, dependant on the units
c0_scale=1.0
if units=='micrometers':
    c0_scale= 1000000.0
elif units=='millimeters':
    c0_scale= 1000.0
elif units=='centimeters':
    c0_scale= 100.0
elif units=='meters':
    c0_scale= 1.0
elif units=='inches':
    c0_scale= 1.0/0.0254
elif units=='mils':
    c0_scale= 1000.0/0.0254


forceScale = 1.0
if forceUnits=='Newtons':
    forceScale= 1.0	## the natural units of the weighted stress tensor
elif forceUnits=='lbf':
    forceScale= 0.2248089
elif forceUnits=='kgf':
    forceScale= 0.1019716


## define the boundary
uo = 4 * np.pi * 10**-7  ## value of magnetic constant
femm.mi_zoom(-2*magnetLength,-1.5*magnetWidth+magnetSep,numrepeats*2*magnetLength,1.5*magnetWidth+magnetSep)	  ## set the window to a nice size for the problem
c0 = c0_scale/(2*uo*magnetLength)
femm.mi_addboundprop("Asymptotic",0,0,0,0,0,0,c0,0,2)		## create the Asymptotic Boundary Condition
femm.mi_addboundprop("Periodic",0,0,0,0,0,0,0,0,4)		## create the Periodic Boundary Condition
femm.mi_addboundprop("Antiperiodic",0,0,0,0,0,0,0,0,5)		## create the anti-Periodic Boundary Condition

femm.mi_addnode(-boundhoriz,4*magnetWidth)
femm.mi_addnode(-boundhoriz,-4*magnetWidth)
femm.mi_addnode(og_thickness*(numrepeats+1)+2*margin+boundhoriz,4*magnetWidth)
femm.mi_addnode(og_thickness*(numrepeats+1)+2*margin+boundhoriz,-4*magnetWidth)
femm.mi_drawline(-boundhoriz,4*magnetWidth,-boundhoriz,-4*magnetWidth) #left boundary
femm.mi_drawline(og_thickness*(numrepeats+1)+2*margin+boundhoriz,4*magnetWidth,
                 og_thickness*(numrepeats+1)+2*margin+boundhoriz,-4*magnetWidth) #right boundary
femm.mi_drawline(-boundhoriz,4*magnetWidth,og_thickness*(numrepeats+1)+2*margin+boundhoriz,4*magnetWidth) #top boundary
femm.mi_drawline(-boundhoriz,-4*magnetWidth,og_thickness*(numrepeats+1)+2*margin+boundhoriz,-4*magnetWidth) #bottom boundary

#adding nodes to analyse through centre
femm.mi_addnode(0,0)
femm.mi_addnode(og_thickness*(numrepeats+1)+2*margin+oddmagnetLength,0)


def drawMagnet(centerx,centery,magnetLength, magnetWidth):
    femm.mi_addnode(centerx-magnetLength/2,centery+magnetWidth/2)
    femm.mi_addnode(centerx+magnetLength/2,centery+magnetWidth/2)
    femm.mi_addnode(centerx-magnetLength/2,centery-magnetWidth/2)
    femm.mi_addnode(centerx+magnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx-magnetLength/2,centery+magnetWidth/2,centerx+magnetLength/2,centery+magnetWidth/2)
    femm.mi_drawline(centerx+magnetLength/2,centery+magnetWidth/2,centerx+magnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx+magnetLength/2,centery-magnetWidth/2,centerx-magnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx-magnetLength/2,centery-magnetWidth/2,centerx-magnetLength/2,centery+magnetWidth/2)
    
def drawoddMagnet(centerx,centery,oddmagnetLength, magnetWidth):
    femm.mi_addnode(centerx-oddmagnetLength/2,centery+magnetWidth/2)
    femm.mi_addnode(centerx+oddmagnetLength/2,centery+magnetWidth/2)
    femm.mi_addnode(centerx-oddmagnetLength/2,centery-magnetWidth/2)
    femm.mi_addnode(centerx+oddmagnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx-oddmagnetLength/2,centery+magnetWidth/2,centerx+oddmagnetLength/2,centery+magnetWidth/2)
    femm.mi_drawline(centerx+oddmagnetLength/2,centery+magnetWidth/2,centerx+oddmagnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx+oddmagnetLength/2,centery-magnetWidth/2,centerx-oddmagnetLength/2,centery-magnetWidth/2)
    femm.mi_drawline(centerx-oddmagnetLength/2,centery-magnetWidth/2,centerx-oddmagnetLength/2,centery+magnetWidth/2)

drawoddMagnet(oddmagnetLength/2+interMagnet+margin,magnetWidth/2+magnetSep,oddmagnetLength, magnetWidth)
drawMagnet(magnetLength/2+oddmagnetLength+2*interMagnet+margin,magnetWidth/2+magnetSep,magnetLength, magnetWidth)
drawoddMagnet(oddmagnetLength/2+magnetLength+oddmagnetLength+3*interMagnet+margin,
            magnetWidth/2+magnetSep,oddmagnetLength, magnetWidth)
drawMagnet(magnetLength/2+oddmagnetLength+magnetLength+oddmagnetLength+4*interMagnet+margin,
            magnetWidth/2+magnetSep,magnetLength, magnetWidth)

drawoddMagnet(oddmagnetLength/2+interMagnet+margin,-(magnetWidth/2+magnetSep),oddmagnetLength, magnetWidth)
drawMagnet(magnetLength/2+oddmagnetLength+2*interMagnet+margin,-(magnetWidth/2+magnetSep),magnetLength, magnetWidth)
drawoddMagnet(oddmagnetLength/2+magnetLength+oddmagnetLength+3*interMagnet+margin,
            -(magnetWidth/2+magnetSep),oddmagnetLength, magnetWidth)
drawMagnet(magnetLength/2+oddmagnetLength+magnetLength+oddmagnetLength+4*interMagnet+margin,
            -(magnetWidth/2+magnetSep),magnetLength, magnetWidth)


## set magnetic properties
#Function attributes: mi_setblockprop('blockname', automesh, meshsize, 'incircuit', magdir, group, turns)
#Top row 
femm.mi_clearselected()
femm.mi_addblocklabel(1.5*magnetLength+2*oddmagnetLength+4*interMagnet+margin,magnetSep+magnetWidth/2)#4th
femm.mi_selectlabel(1.5*magnetLength+2*oddmagnetLength+4*interMagnet+margin,magnetSep+magnetWidth/2)
femm.mi_setblockprop(magnetType, 0, mess_size, "", 90, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(1*magnetLength+1.5*oddmagnetLength+3*interMagnet+margin,magnetSep+magnetWidth/2)#3rd
femm.mi_selectlabel(1*magnetLength+1.5*oddmagnetLength+3*interMagnet+margin,magnetSep+magnetWidth/2)
femm.mi_setblockprop(magnetType, 0, mess_size, "", 180, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(0.5*magnetLength+1*oddmagnetLength+2*interMagnet+margin,magnetSep+magnetWidth/2) #2nd
femm.mi_selectlabel(0.5*magnetLength+1*oddmagnetLength+2*interMagnet+margin,magnetSep+magnetWidth/2)
femm.mi_setblockprop(magnetType, 0, mess_size, "", -90, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(0.5*oddmagnetLength+interMagnet+margin,magnetSep+magnetWidth/2) #1st
femm.mi_selectlabel(0.5*oddmagnetLength+interMagnet+margin,magnetSep+magnetWidth/2)
femm.mi_setblockprop(magnetType, 0, mess_size, "", 0, 0, 0)
femm.mi_clearselected()

#Bottom row
femm.mi_addblocklabel(1.5*magnetLength+2*oddmagnetLength+4*interMagnet+margin,-(magnetSep+magnetWidth/2))#4th
femm.mi_selectlabel(1.5*magnetLength+2*oddmagnetLength+4*interMagnet+margin,-(magnetSep+magnetWidth/2))
femm.mi_setblockprop(magnetType, 0, mess_size, "", 90, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(1*magnetLength+1.5*oddmagnetLength+3*interMagnet+margin,-(magnetSep+magnetWidth/2))#3rd
femm.mi_selectlabel(1*magnetLength+1.5*oddmagnetLength+3*interMagnet+margin,-(magnetSep+magnetWidth/2))
femm.mi_setblockprop(magnetType, 0, mess_size, "", 0, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(0.5*magnetLength+1*oddmagnetLength+2*interMagnet+margin,-(magnetSep+magnetWidth/2)) #2nd
femm.mi_selectlabel(0.5*magnetLength+1*oddmagnetLength+2*interMagnet+margin,-(magnetSep+magnetWidth/2))
femm.mi_setblockprop(magnetType, 0, mess_size, "", -90, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(0.5*oddmagnetLength+interMagnet+margin,-(magnetSep+magnetWidth/2)) #1st
femm.mi_selectlabel(0.5*oddmagnetLength+interMagnet+margin,-(magnetSep+magnetWidth/2))
femm.mi_setblockprop(magnetType, 0, mess_size, "", 180, 0, 0)
femm.mi_clearselected()
  
#-----------------------------------------------------------------------------------------------

#Repeating section
#Function attributes: mi_copytranslate(dx, dy, copies)

femm.mi_selectrectangle(0+margin,magnetSep+magnetWidth+steelThick,
                        2*magnetLength+2*oddmagnetLength+4*interMagnet+margin,-magnetSep-magnetWidth-steelThick, 4)
femm.mi_copytranslate2(og_thickness, 0, numrepeats, 4)
femm.mi_clearselected()

#adding final drift blocks 

femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                    +4*interMagnet)+margin, magnetSep+magnetWidth)
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                    +4*interMagnet)+margin+oddmagnetLength, magnetSep)
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, magnetSep+magnetWidth,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                                +4*interMagnet)+margin+oddmagnetLength, magnetSep+magnetWidth)
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin+oddmagnetLength, magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                                    +4*interMagnet)+margin+oddmagnetLength, magnetSep+magnetWidth)

femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, -magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                    +4*interMagnet)+margin, -(magnetSep+magnetWidth))
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, -magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                    +4*interMagnet)+margin+oddmagnetLength, -magnetSep)
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin, -(magnetSep+magnetWidth),(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                                +4*interMagnet)+margin+oddmagnetLength, -(magnetSep+magnetWidth))
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)
                 +margin+oddmagnetLength, -magnetSep,(numrepeats+1)*(2*magnetLength+2*oddmagnetLength
                                                                    +4*interMagnet)+margin+oddmagnetLength, -(magnetSep+magnetWidth))
femm.mi_clearselected()
femm.mi_addblocklabel((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+margin+oddmagnetLength/2,magnetSep+magnetWidth/2)#4th
femm.mi_selectlabel((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+margin+oddmagnetLength/2,magnetSep+magnetWidth/2)
femm.mi_setblockprop(magnetType, 0, mess_size, "", 0, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+margin+oddmagnetLength/2,-(magnetSep+magnetWidth/2))#4th
femm.mi_selectlabel((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+margin+oddmagnetLength/2,-(magnetSep+magnetWidth/2))
femm.mi_setblockprop(magnetType, 0, mess_size, "", 180, 0, 0)
femm.mi_clearselected()

# draw steel sheet

femm.mi_drawline(margin,magnetSep+magnetWidth+steelThick,
(numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,magnetSep+magnetWidth+steelThick)#outer long line
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,magnetSep+magnetWidth+steelThick,
  (numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,magnetSep+magnetWidth)

femm.mi_drawline(margin,magnetSep+magnetWidth,
(numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,magnetSep+magnetWidth)#inner long line
femm.mi_drawline(margin,-magnetSep-magnetWidth,
(numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,-magnetSep-magnetWidth)  

femm.mi_drawline(margin,-magnetSep-magnetWidth-steelThick,
(numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,-magnetSep-magnetWidth-steelThick)  
femm.mi_drawline((numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,-magnetSep-magnetWidth-steelThick,
(numrepeats+1)*(2*magnetLength+2*oddmagnetLength+4*interMagnet)+2*margin+oddmagnetLength,-magnetSep-magnetWidth)

femm.mi_drawline(margin,magnetSep+magnetWidth+steelThick,margin,magnetSep+magnetWidth)
femm.mi_drawline(margin,-magnetSep-magnetWidth-steelThick,margin,-magnetSep-magnetWidth)

if steelThick > 0:
  femm.mi_addblocklabel(margin+interMagnet+magnetLength,magnetSep + magnetWidth + 0.5*steelThick)
  femm.mi_selectlabel(margin+interMagnet+magnetLength,magnetSep + magnetWidth + 0.5*steelThick)
  femm.mi_addblocklabel(margin+interMagnet+magnetLength,-magnetSep - magnetWidth - 0.5*steelThick)
  femm.mi_selectlabel(margin+interMagnet+magnetLength,-magnetSep - magnetWidth - 0.5*steelThick)
  femm.mi_setblockprop(steelType, 0, mess_size, "", 0, 1, 0)
  femm.mi_clearselected()
  


# set surrounding properties
femm.mi_addblocklabel(10, 0)
femm.mi_selectlabel(10 ,0)
femm.mi_setblockprop("Air", 0, mess_size, "", 0, 1, 0)
femm.mi_clearselected()

femm.mi_setblockprop("Air", 1, 0, "", 0, 1, 0)

femm.mi_refreshview()

##----------------------------------------------------------------------------------------------
femm.mi_saveas('wiggle_attempt4.FEM')
