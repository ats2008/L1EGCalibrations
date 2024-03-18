import ROOT
import array
import math
import operator

compressedIetaFile  = "data/egCompressEtaLUT_4bit_v4.txt"
compressedEFile     = "data/egCompressELUT_4bit_v4.txt"
compressedShapeFile = "data/egCompressShapesLUT_calibr_4bit_v4.txt"
inputFileName       = "./regressionTrainerFile.root"
outputFileName      = "compressed_regressionTrainerFile.root"
treeName = "eIDSimpleTree"


def readLUT(lutFileName):
    lut = {}
    with open(lutFileName, 'r') as f:
        lines = f.readlines()
        for line in lines:
            if line[0]=='#': continue
            tokens = line.split()
            if len(tokens)!=2: continue
            lut[int(tokens[0])] = int(tokens[1])
    return lut

def sortShapes(shapeHisto):
    numbers = {}
    nbins = shapeHisto.GetNbinsX()
    for b in range(1,nbins+1):
        numbers[b-1] = shapeHisto.GetBinContent(b)
    #
    sortedShapes = sorted(numbers.items(), key=operator.itemgetter(1))
    #
    zeros = []
    nonzeros = []
    for shape,n in sortedShapes:
        if n==0:
            zeros.append((shape,n))
        else:
            nonzeros.append((shape,n))
    #
    nonzeros.reverse()
    sortedShapes = nonzeros
    sortedShapes.extend(zeros)
    #
    lut = {}
    sortedShape = 0
    for shape,n in sortedShapes:
        lut[shape] = sortedShape
        sortedShape += 1
    #
    with open("data/compressedSortedShapes.txt", 'w') as f:
        for shape in range(0,128):
            sortedShape = lut[shape]
            f.write(str(shape) + ' ' + str(sortedShape) +'\n' )
    return lut


## Read Ieta,E,shape compression mapping
compressedIeta  = readLUT(compressedIetaFile)
compressedE     = readLUT(compressedEFile)
compressedShape = readLUT(compressedShapeFile)

## Reading input trees 
## + filling compressed shape histo
inputFile = ROOT.TFile.Open(inputFileName)
inputTree = inputFile.Get(treeName)
inputTree.__class__ = ROOT.TTree
shapeHisto = ROOT.TH1F("compressedShapeHisto", "compressedShapeHisto", 128, -0.5, 127.5)


data = {"Run"    :array.array('i',[0]),
        "Event"  :array.array('l',[0]),
        #"Weight"  :array.array('f',[0.]),
        #"group"   :array.array('i',[0]),
        "ieta"   :array.array('i',[0]),
        #"iphi"    :array.array('i',[0]),
        "E"      :array.array('i',[0]),
        "shape"  :array.array('i',[0]),
        #"side"    :array.array('i',[0]),
        "target" :array.array('f',[0.]),
        #"ptoff"   :array.array('f',[0.]),
        #"etl1"    :array.array('f',[0.]),
        "Run2IdLevel"  :array.array('i',[0]),
        "compressedieta" :array.array('i',[0]),
        "compressedE" :array.array('i',[0]),
        "compressedshape" :array.array('i',[0]),
        "compressedsortedshape" :array.array('i',[0]),
       }


print( "First pass: reading tree to build compressed shape histo")
nentries = inputTree.GetEntriesFast()
for e in range(nentries):
    if e%2000==0:
        print( "\t entry = ",e," / ",nentries,"  [ ",100.0*e/nentries ," ] ")
    inputTree.GetEntry(e)
    data["Run"][0]    = int(inputTree.Run)
    data["Event"][0]  =  int(inputTree.Event)
    #data["Weight"][0] =  inputTree.Weight
    #data["group"][0]  =  int(inputTree.group)
    data["ieta"][0]   =  int(inputTree.ieta)
    #data["iphi"][0]   =  int(inputTree.iphi)
    data["E"][0]      =  int(inputTree.E)
    data["shape"][0]  =  int(inputTree.shape)
    #data["side"][0]   =  int(inputTree.side)
    data["target"][0] =  inputTree.target
    #data["ptoff"][0]  =  inputTree.ptoff
    #data["etl1"][0]   =  inputTree.etl1
    data["compressedieta"][0]  = int(math.copysign(compressedIeta[abs(data["ieta"][0])], data["ieta"][0]))
    data["compressedE"][0]     = compressedE[min(data["E"][0],255)]
    data["compressedshape"][0] = compressedShape[data["shape"][0]]
    shapeHisto.Fill(data["compressedshape"][0])
    #outputTree.Fill()
## Sort compressed shapes and write in file
compressedSortedShape = sortShapes(shapeHisto)

## TODO
## Reading and filling tree with compressed and sorted values
outputFile = ROOT.TFile.Open(outputFileName, "RECREATE")
outputTree = ROOT.TTree(treeName, treeName)
for name, a in data.items():
    outputTree.Branch(name, a, "{0}/{1}".format(name, a.typecode.upper()))
shapeHisto.Write()

print( "Second pass: reading tree for filling output tree")
for e in range(nentries):
    if e%2000==0:
        print( "\t entry = ",e," / ",nentries,"  [ ",100.0*e/nentries ," ] ")
    inputTree.GetEntry(e)
    data["Run"][0]    = int(inputTree.Run)
    data["Event"][0]  =  int(inputTree.Event)
    #data["Weight"][0] =  inputTree.Weight
    #data["group"][0]  =  int(inputTree.group)
    data["ieta"][0]   =  int(inputTree.ieta)
    #data["iphi"][0]   =  int(inputTree.iphi)
    data["E"][0]      =  int(inputTree.E)
    data["shape"][0]  =  int(inputTree.shape)
    #data["side"][0]   =  int(inputTree.side)
    data["target"][0] =  inputTree.target
    #data["ptoff"][0]  =  inputTree.ptoff
    #data["etl1"][0]   =  inputTree.etl1
    data["Run2IdLevel"][0] = inputTree.Run2IdLevel

    data["compressedieta"][0]  = int(math.copysign(compressedIeta[abs(data["ieta"][0])], data["ieta"][0]))
    data["compressedE"][0]     = compressedE[min(data["E"][0],255)]
    data["compressedshape"][0] = compressedShape[data["shape"][0]]
    data["compressedsortedshape"][0] = compressedSortedShape[data["compressedshape"][0]]
    outputTree.Fill()


outputFile.cd()
outputTree.Write()
outputFile.Close()
inputFile.Close()

