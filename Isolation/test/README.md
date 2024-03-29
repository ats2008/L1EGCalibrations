# Deriving Isolation LUTs
## Basic Setup
```
cmsrel CMSSW_X_Y_Z   # Any X > 12 will do
cd CMSSW_X_Y_Z
git clone https://github.com/ats2008/L1EGCalibrations.git
cd L1EGCalibrations/Isolation/test/
```
## Introduction
#### 1. Computation of the isolation cut for a given flat efficiency.
The isolation energy (E_iso) is defined as the difference between the energy deposit in a 6x9 TT window (E_6x9) and the raw energy deposited by the tau/EG (ET_raw) . E_iso depends of ET_raw, iEta and nTT (being this last one a pileup estimator for taus and EG). We reject the background by requiring that the tau is isolated, e.g., this E_iso is smaller than a certain value (E_isocut). The value of the isolation cut is chosen accordingly to a target signal efficiency, that is constant in ieta and nTT but can be varied as a function of ET. It is computed for a given flat efficiency with simulated events in exclusive intervals of the isolation input variables ET_raw, iEta and nTT.
The script Build_Isolation_WPs_MC_NewCompression_Thomas_nTT_OlivierFlatWP_With2017Layer1.C computes this isolation cuts and gives as output a root file (named Iso_LUTs_Distributions_MC or something similar) which contains histograms labelled by (i,j,k) corresponding to the value of the isocut for each bin (ieta, iEt supercompr., nTT). (For the taus we previously do a supercompression of the Et bins to have more statistics in each bin.)
#### 2. Relaxation of the isocut with pT
After that, we relax the isolation cut with pT to reach a 100% plateau efficiency while controlling the trigger rate.  We do the relaxation linearly, starting from a flat efficiency A up to a certain pT B, then relaxing linearly up to a certain pT C, when the 100% efficiency is reached (see picture). We test different isolation options corresponding to different combinations of the values A/B/C (this is what we call Option 1, 2…). 
To derive these new isocuts for each of the options you can use the script Fill_Isolation_TH3_MC_2017Layer1Calibration.C. 
This gives as ouput another root file (Iso_LUTs_Options_MC or similar) that gives you histograms with the value of the isocuts for each Option for (x,y,z = ieta, iEt supercompr., nTT).

#### 3. Rate estimation:
Now one should check the rates of the different isolation options as a function of the ET threshold. For this we use ZeroBias samples (which you should previously calibrate). 
Make sure to modify the correct number of bunches for the run you are considering (see here for example https://cmswbm.cern.ch/cmsdb/servlet/FillReport?FILL=6311), as well as the luminosity. 
The rate is computed as:
    N_tot = number of zero bias events
    N_pass = number of zero bias events that pass the trigger that we are testing (L1_DoubleIsoTau_30)
    Rate = (N_pass/PS)/(N_total/PS) x scale
    Scale = frequency at which you see a zero bias event at the LHC (Hz) 
          = (bunch revolution frequency) x (#bunches)
          = (c/27km) x #bunches = 11245.6 x #bunches (1900-2600)    (multiply by 0.001 to get kHz)

#### 4. Apply the isolation:
The last step is to apply the isolation (relaxed) cuts onto your signal Tag and Probe  samples, and get the turnons.

#### 5. Make the isolation LUT:
You can produce the isolation LUT (text file) with the script **TO BE FILLED** .

## Step 1 and 2 : Making the Working points and filling the Isolations maps
  * For compiling :
  ```
  make isolationAnalysis
  ```
  * Execution
    The `step 1` and `step 2` can be performed separtely by passsing `do_1` or `do_2` or `do_all` as the second argument to `isolationAnalysis.exe`.
  ```
  ./isolationAnalysis.exe ParList_IsoAna.dat do_1
  ```

  Since `step 1` takes a long time to run, is a covininat methord to do `step 1` once and reuse its output for defining the LUT histogram maps.

  The Configuration file is described as below :
  ```
    NtupleFileName=<file path to the list of Tag and Probe Ntuples.>
    OutputFileName=<Name of the output filename>
    EtLUTFileName=data/compressionLuts/egCompressELUT_4bit_v4.txt
    EtaLUTFileName=data/compressionLuts/egCompressEtaLUT_4bit_v4.txt
    NTTLUTFileName=data/compressionLuts/tauCompressnTTLUT_5bit_v8.txt
    SCEtaLUTFileName=<Doesnot Care>
    SCNTTLUTFileName=<Doesnot Care>
    OutputWPStepFileName=<Path to the output file of step 1 , to be used if you invoke do_2 >
    MaxEntries=-1
    LUTProgressionOptions=<lut description string>
    ReportEvery=5000
  ```
  **A <lut description string> for the grid search can be obtained by using the custon script available at `misc/makeOptionsGrid.py`**
  ```
  python misc/makeOptionsGrid.py
  ```
  Will print out a custon LUT description string as described in the script.

##### FIT Diagonstics for step 1 :

- Plots the fits for a given eta and et bins for different efficiencies
Uasge :
```
python python/makeFitDiagonosticPlots_forEtEta.py -i /home/aravind/Documents/dumpX/EraGDefault_Step1.root  -t defaultV1 --et 5 --eta 9
```

- Plots the fits for each eta amd et bins for a given efficiency
Uasge :
```
 python python/makeFitDiagonosticPlots.py -i /home/aravind/Documents/dumpX/EraGDefault_Step1.root  -t defaultV1 -e 50
```

#### Step 2 In Condor
```
     python misc/makeOptionsGrid.py -c misc/Par_Step1Step2.dat.tpl.cfg --opt CalibFiles/EraGDefault_Step1.root  -t defV4 -r misc/runStep2.tpl.sh --doStep2
```

## Step 3 and 4 : Obtaining Rate and Efficincies
  * For compiling :
  ```
  make applyIsolation
  ```
  * Execution [not recomended if you have large LUTs under consideration , or have big TandP set to measure]
### Using condor for grid search
 -  Customize the `misc/makeOptionsGrid.py` to point to the correct files
    - Need to customize `misc/Par_ApplyIsolation.dat.tpl.cfg`
        - Need to customize `config/Ntuple_turnon.txt`
        - Need to customize `config/Ntuple_rate.txt`
 - Make the Condor job files 
   ```
   # For default scheme
   python misc/makeOptionsGrid.py -c misc/Par_ApplyIsolation.dat.tpl.cfg  --opt CalibFiles/EraG_step2Output.root-t eraGExtrapolated  
   # For extrapolation scheme
   python misc/makeOptionsGrid.py -c misc/Par_ApplyIsolation.dat.tpl.cfg  --opt CalibFiles/EraG_step2Output.root-t eraGExtrapolated  --doV2
   ```
## Step 7 : Deriving the optimal LUTs
Use the optimization scripts available under python directory
*Please setup the baseline files insside the scripts properly !*
```
python python/doubleEGOptimization.py -f eraGExtrapolated_v4_upd1/eraGExtrapolated_v4_upd1.fls  -o results/optimization/
```
```
python python/singleEGOptimization.py -f eraGExtrapolated_v4_upd1/eraGExtrapolated_v4_upd1.fls  -o results/optimization/
```
  * One can also use the notebooks in `analyzer` directory
  * TODO : should move to a json based customization scheme

##  Export the Isolation results for a set of options

```
make exportIsolation
./exportIsolation.exe config/ParList_Export.dat
```

## Making the visualization of the LUTs exported to txt files 
```
# For EG LUTs
python python/getLUTImages.py -t LUT_TAG_HERE -i <input file path >

# For Tau LUTs
python python/getTauLUTImages.py -t LUT_TAG_HERE -i <input file path >

```
# Augmenting to make dataset flat in nTT

 - See the `analyzers/Augumentation.ipynb` for more details
 - use the script `python/augumentDataset_flatInNTT.py` for making the augmented dataset
    ```bash
     python python/augumentDataset_flatInNTT.py -i workarea/DYToLL_TandP_calo_v6_ZS0p0_v1.root -o workarea/results/dataAugmentation_mc_v1/ --export
    ```
