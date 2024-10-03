from processes import (ProfileZ, 
                       UniformGrid1D, 
                       Erosion, 
                       SeaLevel, 
                       VerticalDisp, 
                       WaterHeight, 
                       InitTopo, 
                       SLRise, 
                       SLFile,
                       Construction,
                       MyHorizontalFactor, 
                       MyVerticalFactor, 
                       ErosiveMemory,
                       SedimClastics,
                       InitTopoTerr,
                      )
import xsimlab as xs

#*********************************************************************
#*********************************************************************

reef = xs.Model({
    "vertical"  : VerticalDisp,
    "grid"      : UniformGrid1D,
    "water"     : WaterHeight,
    "profile"   : ProfileZ,
    "init"      : InitTopo,
    "sealevel"  : SeaLevel,
    "SLstory"   : SLFile,
    "construct" : Construction,
    "hfactor"   : MyHorizontalFactor,
    "vfactor"   : MyVerticalFactor,
    "eros"      : Erosion,
    "erosmem"   : ErosiveMemory,
    "depot"     : SedimClastics
})

#*********************************************************************

reef_platform = xs.Model({
    "vertical"  : VerticalDisp,
    "grid"      : UniformGrid1D,
    "water"     : WaterHeight,
    "profile"   : ProfileZ,
    "init"      : InitTopoTerr,
    "sealevel"  : SeaLevel,
    "SLstory"   : SLFile,
    "construct" : Construction,
    "hfactor"   : MyHorizontalFactor,
    "vfactor"   : MyVerticalFactor,
    "eros"      : Erosion,
    "erosmem"   : ErosiveMemory,
    "depot"     : SedimClastics
})