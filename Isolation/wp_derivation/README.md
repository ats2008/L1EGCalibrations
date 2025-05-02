# Isolation WP Derivations

[TOC] 

## Introduction

This module derives the Isolation working points for the EG objects from TagAndProbe ntuple . Need to add more information here !

The module can also be used to emulate a give LUT ( in cmssw-complacent txt format) over TagAndProbe Ntuples as well as L1Ntuples.

## Usage

### Derive Isolation LUTs

```bash
[]$ python python/makeIsolationLUTs.py -h
usage: makeIsolationLUTs.py [-h] [-i INPUTFILE] [-o DEST] [-p]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        Input File
  -o DEST, --dest DEST  destination To Use
  -p, --doSanityPlots   export the basic sanity plots

```
- `Input File` is the TagAndProbe Ntuple with re-emulated branches filled.
- `DEST` is the destination folder where the derived LUTs are going to be exported. Plots will also be saved here.
- `doSanityPlots` will make a set of plots for the input varable distributions

### Export the LUTs for Efficiencies
```bash
[]$ python python/ApplyIsolationLUTForEG_Rates.py -h
usage: ApplyIsolationLUTForEG_Rates.py [-h] [-i INPUTFILE] [-l INPUTLUT] [-o DEST] [-t TAG] [--iso {loose,tight}] [-p] [-e]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        Input File
  -l INPUTLUT, --inputLUT INPUTLUT
                        Input Isolation LUT File
  -o DEST, --dest DEST  destination To Use
  -t TAG, --tag TAG     tag To the made files
  --iso {loose,tight}   choose isolation bit
  -p, --doSanityPlots   export the basic sanity plots
  -e, --doEmulated      Use emulated branches

```
- `Input File` is the L1 Ntuple
- `DEST` is the destination folder where the derived LUTs are going to be exported. Plots will also be saved here.
- `INPUTLUT`  path to the LUT in the txt format. If the LUT is not provided the emulated ( or unpacked ) branches are used to make the efficincies. This can help 1-to-1  comparisons.
- `TAG` : additional prefix to the saved files
- `iso` : in case no `INPUTLUT` is provided, choose the loose or tight isolation bit
- `doEmulated` : Use emulated branches for the studies. Note : need the this option if proper isolation calcluateion has to be done for the `INPUTLUT`
- `doSanityPlots` will make a set of plots for the input varable distributions

Sample Usage
```bash
    python3 python/ApplyIsolationLUTForEG_Rates.py -l IsoDerivation/loose_LUT.txt -i workarea/zs_studies/L1Ntuple_default.root -e -o results/loose_iso/rates/ -p
```


### Export the LUTs for Efficiencies
```bash
[] $python python/ApplyIsolationLUTForEG_Efficiency.py -h
usage: ApplyIsolationLUTForEG_Efficiency.py [-h] [-i INPUTFILE] [-l INPUTLUT] [-o DEST] [-t TAG] [--iso {loose,tight}] [-p] [-e]

optional arguments:
  -h, --help            show this help message and exit
  -i INPUTFILE, --inputFile INPUTFILE
                        Input File
  -l INPUTLUT, --inputLUT INPUTLUT
                        Input Isolation LUT File
  -o DEST, --dest DEST  destination To Use
  -t TAG, --tag TAG     tag To the made files
  --iso {loose,tight}   choose isolation bit
  -e, --doEmulated      Use emulated branches
```
- `Input File` is the TagAndProbe Ntuple with re-emulated branches filled.
- `DEST` is the destination folder where the derived LUTs are going to be exported. Plots will also be saved here.
- `INPUTLUT`  path to the LUT in the txt format. If the LUT is not provided the emulated ( or unpacked ) branches are used to make the efficincies. This can help 1-to-1  comparisons.
- `TAG` : additional prefix to the saved files
- `iso` : in case no `INPUTLUT` is provided, choose the loose or tight isolation bit
- `doEmulated` : Use emulated branches for the studies. Note : need the this option if proper isolation calcluateion has to be done for the `INPUTLUT`
- `doSanityPlots` will make a set of plots for the input varable distributions

Sample usage

```bash
    python3 python/ApplyIsolationLUTForEG_Efficiency.py -l IsoDerivation/loose_LUT.txt -i workarea/zs_studies/TandP_default.root -e -o results/loose_iso/effs/
```

