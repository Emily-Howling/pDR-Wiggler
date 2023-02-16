#from https://www.femm.info/wiki/HalbachExample
#Modified Halbach array by Jose B. Almeida, November 2017

import femm
import matplotlib.pyplot as plt
import numpy as np

# The package must be initialized with the openfemm command.
femm.openfemm()

# modify these values to suit your needs
units = 'millimeters'	# should be one of "inches","millimeters","centimeters","mils","meters","micrometers" 
forceUnits = 'Newtons'		# should be one of "lbf", "Newtons", "kgf"
magnetLength = 20.0 #horizontal
magnetThick = 10.0 #vertical
magnetSep = 7  # one half of separation between magnet layers vertically
magnetType= 'N52'  # has to be the EXACT name of a material in the materials library
steelThick = 1  # thickness of magnetic circuit steel
interMagnet = 0.7  # thickness of steel inter magnet pieces.
steelType = '1020 Steel'  # has to be the EXACT name of a material in the materials library
numrepeats=3
#-----------------------------------------------------------------------------------------------


# We need to create a new Magnetostatics document to work on.
femm.newdocument(0)

#addpath('C:\femm42\mfiles')

femm.mi_hidegrid()


## define the problem type, equivalent to the top level menu    Problem
#Function attributes: mi_probdef(freq,units,type,precision,depth,minangle,(acsolver))
femm.mi_probdef(0, units, 'planar', 1E-8, magnetLength)

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
femm.mi_zoom(-2*magnetThick,-1.5*magnetThick,2*magnetThick,1.5*magnetThick)	  ## set the window to a nice size for the problem
c0 = c0_scale/(2*uo*magnetThick)
femm.mi_addboundprop("Asymptotic",0,0,0,0,0,0,c0,0,2)		## create the Asymptotic Boundary Condition
femm.mi_addboundprop("Periodic",0,0,0,0,0,0,0,0,4)		## create the Periodic Boundary Condition
femm.mi_addboundprop("Antiperiodic",0,0,0,0,0,0,0,0,5)		## create the anti-Periodic Boundary Condition

left_edge=-1.5*magnetThick-2*interMagnet-0.1
og_thickness=4*magnetThick+4*interMagnet+0.2

femm.mi_addnode(-1.5*magnetThick-2*interMagnet-0.1,3*magnetThick)
femm.mi_addnode(-1.5*magnetThick-2*interMagnet-0.1,-3*magnetThick)
femm.mi_addnode(og_thickness*(numrepeats+1)+left_edge,3*magnetThick)
femm.mi_addnode(og_thickness*(numrepeats+1)+left_edge,-3*magnetThick)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet-0.1,3*magnetThick,-1.5*magnetThick-2*interMagnet-0.1,-3*magnetThick)
femm.mi_drawline(og_thickness*(numrepeats+1)+left_edge,3*magnetThick,og_thickness*(numrepeats+1)+left_edge,-3*magnetThick)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet-0.1,3*magnetThick,og_thickness*(numrepeats+1)+left_edge,3*magnetThick)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet-0.1,-3*magnetThick,og_thickness*(numrepeats+1)+left_edge,-3*magnetThick)

#adding nodes to analyse through centre
femm.mi_addnode(-1.5*magnetThick-2*interMagnet-0.1,0)
femm.mi_addnode(og_thickness*(numrepeats+1)+left_edge,0)

#NOT SURE WHAT BOUNDARIES ARE NEEDED SO OMITTED FOR NOW
## set boundary properties

# femm.mi_selectsegment(0,2*magnetThick)
# femm.mi_selectsegment(0,-2*magnetThick)
# femm.mi_setsegmentprop("Asymptotic",0,0,0,0)
# femm.mi_clearselected()
# femm.mi_selectsegment(-1.5*magnetThick-2*interMagnet-0.1,0)
# femm.mi_selectsegment(2.5*magnetThick+2*interMagnet,0)
# femm.mi_setsegmentprop("Periodic",0,0,0,0)
# femm.mi_clearselected()

## draw magnets

def drawMagnet(centerx,centery,magnetThick):
    femm.mi_addnode(centerx-magnetThick/2,centery+magnetThick/2)
    femm.mi_addnode(centerx+magnetThick/2,centery+magnetThick/2)
    femm.mi_addnode(centerx-magnetThick/2,centery-magnetThick/2)
    femm.mi_addnode(centerx+magnetThick/2,centery-magnetThick/2)
    femm.mi_drawline(centerx-magnetThick/2,centery+magnetThick/2,centerx+magnetThick/2,centery+magnetThick/2)
    femm.mi_drawline(centerx+magnetThick/2,centery+magnetThick/2,centerx+magnetThick/2,centery-magnetThick/2)
    femm.mi_drawline(centerx+magnetThick/2,centery-magnetThick/2,centerx-magnetThick/2,centery-magnetThick/2)
    femm.mi_drawline(centerx-magnetThick/2,centery-magnetThick/2,centerx-magnetThick/2,centery+magnetThick/2)
#endfunction

drawMagnet(-magnetThick-interMagnet,magnetThick/2+magnetSep,magnetThick)
drawMagnet(0,magnetThick/2+magnetSep,magnetThick)
drawMagnet(magnetThick+interMagnet,magnetThick/2+magnetSep,magnetThick)
drawMagnet(2*magnetThick+2*interMagnet,magnetThick/2+magnetSep,magnetThick)
femm.mi_addnode(-1.5*magnetThick - 2*interMagnet,magnetSep + magnetThick + steelThick)
femm.mi_addnode(-1.5*magnetThick - 2*interMagnet,magnetSep)
femm.mi_drawline(-1.5*magnetThick - 2*interMagnet,magnetSep + magnetThick + steelThick,-1.5*magnetThick - 2*interMagnet,magnetSep)


drawMagnet(-magnetThick-interMagnet,-magnetThick/2-magnetSep,magnetThick)
drawMagnet(0,-magnetThick/2-magnetSep,magnetThick)
drawMagnet(magnetThick+interMagnet,-magnetThick/2-magnetSep,magnetThick)
drawMagnet(2*magnetThick+2*interMagnet,-magnetThick/2-magnetSep,magnetThick)
femm.mi_addnode(-1.5*magnetThick - 2*interMagnet,-magnetSep - magnetThick - steelThick)
femm.mi_addnode(-1.5*magnetThick - 2*interMagnet,-magnetSep)
femm.mi_drawline(-1.5*magnetThick - 2*interMagnet,-magnetSep - magnetThick - steelThick,-1.5*magnetThick - 2*interMagnet,-magnetSep)

# now close intermagnet spaces
femm.mi_drawline(-1.5*magnetThick - 2*interMagnet,magnetSep + magnetThick,
  2.5*magnetThick+2*interMagnet,magnetThick+magnetSep)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet,magnetSep,2.5*magnetThick+2*interMagnet,magnetSep)

femm.mi_drawline(-1.5*magnetThick-2*interMagnet,-magnetThick-magnetSep,
  2.5*magnetThick+2*interMagnet,-magnetThick-magnetSep)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet,-magnetSep,2.5*magnetThick+2*interMagnet,-magnetSep)

# draw steel sheet

femm.mi_drawline(-1.5*magnetThick-2*interMagnet,magnetSep+magnetThick+steelThick,
2.5*magnetThick+2*interMagnet,magnetSep+magnetThick+steelThick)
  
  
femm.mi_drawline(2.5*magnetThick+2*interMagnet,magnetSep+magnetThick+steelThick,
 2.5*magnetThick+2*interMagnet,magnetSep+magnetThick)
femm.mi_drawline(-1.5*magnetThick-2*interMagnet,-magnetSep-magnetThick-steelThick,
2.5*magnetThick+2*interMagnet,-magnetSep-magnetThick-steelThick)  
    
femm.mi_drawline(2.5*magnetThick+2*interMagnet,-magnetSep-magnetThick-steelThick,
2.5*magnetThick+2*interMagnet,-magnetSep-magnetThick)

## set magnetic properties
#Function attributes: mi_setblockprop('blockname', automesh, meshsize, 'incircuit', magdir, group, turns)
femm.mi_clearselected()
femm.mi_addblocklabel(2*magnetThick+2*interMagnet,magnetSep+magnetThick/2)
femm.mi_selectlabel(2*magnetThick+2*magnetSep,magnetSep+magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 90, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(-magnetThick-interMagnet,magnetSep+magnetThick/2)
femm.mi_selectlabel(-magnetThick-interMagnet,magnetSep+magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 0, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(magnetThick+interMagnet,magnetSep+magnetThick/2)
femm.mi_selectlabel(magnetThick+interMagnet,magnetSep+magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 180, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(-0.5,magnetSep+magnetThick/2)
femm.mi_selectlabel(-0.5,magnetSep+magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", -90, 0, 0)
femm.mi_clearselected()

femm.mi_addblocklabel(2*magnetThick+2*interMagnet,-magnetSep-magnetThick/2)#bottom far right
femm.mi_selectlabel(2*magnetThick+2*magnetSep,-magnetSep-magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 90, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(-magnetThick-interMagnet,-magnetSep-magnetThick/2)#bottom far left 
femm.mi_selectlabel(-magnetThick-interMagnet,-magnetSep-magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 180, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(magnetThick+interMagnet,-magnetSep-magnetThick/2)#bottom left
femm.mi_selectlabel(magnetThick+interMagnet,-magnetSep-magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", 0, 0, 0)
femm.mi_clearselected()
femm.mi_addblocklabel(0,-magnetSep-magnetThick/2)#bottom right
femm.mi_selectlabel(0,-magnetSep-magnetThick/2)
femm.mi_setblockprop(magnetType, 0, 0.5, "", -90, 0, 0)
femm.mi_clearselected()

# set block properties

femm.mi_addblocklabel(0,interMagnet/2 + 1.8*magnetThick)
femm.mi_selectlabel(0,interMagnet/2 + 1.8*magnetThick)
femm.mi_setblockprop("Air", 1, 1, "", 0, 1, 0)
femm.mi_clearselected()

if interMagnet > 0:
  femm.mi_addblocklabel(-1.5*magnetThick-1.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_addblocklabel(-0.5*magnetThick-0.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_addblocklabel(1.5*magnetThick+1.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_addblocklabel(0.5*magnetThick+0.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_selectlabel(-1.5*magnetThick-1.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_selectlabel(-0.5*magnetThick-0.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_selectlabel(1.5*magnetThick+1.5*interMagnet,magnetSep+magnetThick/2)
  femm.mi_selectlabel(0.5*magnetThick+0.5*interMagnet,magnetSep+magnetThick/2)
  
  femm.mi_addblocklabel(-1.5*magnetThick-1.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_addblocklabel(-0.5*magnetThick-0.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_addblocklabel(1.5*magnetThick+1.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_addblocklabel(0.5*magnetThick+0.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_selectlabel(-1.5*magnetThick-1.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_selectlabel(-0.5*magnetThick-0.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_selectlabel(1.5*magnetThick+1.5*interMagnet,-magnetSep-magnetThick/2)
  femm.mi_selectlabel(0.5*magnetThick+0.5*interMagnet,-magnetSep-magnetThick/2)
  
  femm.mi_setblockprop("Air", 1, 1, "", 0, 1, 0)
  femm.mi_clearselected()


if steelThick > 0:
  femm.mi_addblocklabel(0,magnetSep + magnetThick + 0.5*steelThick)
  femm.mi_selectlabel(0,magnetSep + magnetThick + 0.5*steelThick)
  femm.mi_addblocklabel(0,-magnetSep - magnetThick - 0.5*steelThick)
  femm.mi_selectlabel(0,-magnetSep - magnetThick - 0.5*steelThick)
  femm.mi_setblockprop(steelType, 1, 0.1, "", 0, 1, 0)
  femm.mi_clearselected()
  
  
#Repeating section
#could do it using groups but trying with rectangle
#Function attributes: mi_copytranslate(dx, dy, copies)

#femm.mi_selectgroup()
femm.mi_selectrectangle(-1.5*magnetThick-2*interMagnet-10,magnetSep+magnetThick+steelThick+10,
                        2.5*magnetThick+2*interMagnet+10,-magnetSep-magnetThick-steelThick-10, 4)
femm.mi_copytranslate2(og_thickness, 0, numrepeats, 4)
femm.mi_clearselected()

femm.mi_refreshview()


##----------------------------------------------------------------------------------------------
femm.mi_saveas('wiggle_attempt1.FEM')



##----------------------------------------------------------------------------------------------

# start processing

femm.mi_analyze
femm.mi_loadsolution # opens the solution window
#femm.mo_showdensityplot(1,0,2.0,0.0,'mag') ##Shows the flux density
# femm.mo_addcontour(-1.5*magnetThick-2*interMagnet,0)
# femm.mo_addcontour(2.5*magnetThick+2*interMagnet,0)
# femm.mo_makeplot(2,200)


