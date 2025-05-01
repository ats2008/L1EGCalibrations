#!/usr/bin/env python
import sys,argparse
import ROOT
from GBR2LUTEmulator import GBR2LUTEmulator

parser = argparse.ArgumentParser()
parser.add_argument('-i',"--inputFile", help="Input File")
args = parser.parse_args()

inputFile =  args.inputFile
print(inputFile)

version = "v17.04.04"
sortedShapesFile = "compressedSortedShapes.txt"
regresionName ="pfetReg"
outputDir="./"
suffix="_{}.txt".format(version)

print( "Reading the Input Result file File as   :  ", inputFile)
print( "Reading the Sorted Shapes File as       :  ", sortedShapesFile)
print( "Setting the output file name as          :  ", outputDir+regresionName+suffix )


header = """\
# Calibration vs |ieta|,shape,E. Derived from Run 283478 data, with semi-parametric regression
# The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits
# Index is compressedShape+compressedE<<4+compressedIeta<<8. 
# Compression version is v4
#anything after # is ignored with the exception of the header
#the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header>
#<header> V8 12 10 </header>
"""

ietapoints = []
shapepoints = []
Epoints = []
for i in range(0,16):
    ietapoints.append(i)

for i in range(0,16):
    shapepoints.append(i)

for i in range(0,16):
    Epoints.append(i)

gbr2luts = []

gbr2luts.append(GBR2LUTEmulator())
gbr2luts[-1].name = regresionName
gbr2luts[-1].inputFileName = inputFile
gbr2luts[-1].outputFileName = outputDir+gbr2luts[-1].name+ suffix
gbr2luts[-1].sortedShapes = sortedShapesFile
gbr2luts[-1].variablePoints.append(("abs(compressedieta)",ietapoints))
gbr2luts[-1].variablePoints.append(("compressedE",Epoints))
gbr2luts[-1].variablePoints.append(("compressedsortedshape",shapepoints))
gbr2luts[-1].retrieveForest()
gbr2luts[-1].createLUT(header)

