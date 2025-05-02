import uproot,json,os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
import itertools,argparse
import awkward as ak
import util as utl

BPTX_RATE=26.5
hep.style.use("CMS")
plt.figure()
plt.close()
matplotlib.rcParams['figure.figsize'] = (6, 5.5)
matplotlib.rcParams['figure.dpi'] = 100

BINS=np.arange(0.0,100.1,1.0)-0.5
BIN_X=0.5*(BINS[1:]+BINS[:-1])
BINS=[float(i) for i in BINS]
BIN_X=[float(i) for i in BIN_X]

def main():
    
    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--inputFile", help="Input File")
    parser.add_argument('-l',"--inputLUT", help="Input Isolation LUT File")
    parser.add_argument("-o","--dest", help="destination To Use", default='results/' )
    parser.add_argument("-t","--tag" , help="tag To the made files", default=None )
    parser.add_argument("--iso", help="choose isolation bit ", choices=['loose','tight'] )
    parser.add_argument("-p","--doSanityPlots", help="export the basic sanity plots",default=False,action='store_true')
    parser.add_argument("-e","--doEmulated", help="Use emulated branches",default=False,action='store_true')
    args = parser.parse_args()

    prefix=args.dest
    fname=args.inputFile
    
    print("Input Ntuple : ",fname)
    print("Input LUT : ",args.inputLUT)
    print("Destination : ",prefix)
    
    cmd=f'mkdir -p {prefix}'
    print(cmd); os.system(cmd)
    
    doEmulated=args.doEmulated
    
    fileIn= uproot.open(fname)

    data=fileIn['Ntuplizer/TagAndProbe']
    
    dataStore=data.arrays(['l1tEmuPt','l1tEmuNTT','l1tEmuRawEt','RunNumber','l1tEmuTowerIEta','l1tEmuEta',
                           'l1tPt','l1tEta','l1tIso','Nvtx','l1tEmuIsoEt','l1tEmuIso',
                           'eleProbeSclEt','isProbeLoose','eleProbePt',
                           'eleProbeEta','eleProbePhi',
                           'eleTagEta','eleTagPhi'])
    
    print("Number of events : ",len(dataStore))
    
    dataStore['dR']= np.sqrt( 
                             ( dataStore['eleProbeEta'] - dataStore['eleTagEta'] )**2 + 
                             ( dataStore['eleProbePhi'] - dataStore['eleTagPhi'] )**2 
                            ) 
    
    preselection_mask = np.abs(dataStore['eleProbeEta']) < 2.5
    preselection_mask = np.logical_and( preselection_mask , dataStore['dR'] > 0.6 )
    preselection_mask = np.logical_and( preselection_mask , dataStore['isProbeLoose'] > 0.5 ) 
    print("Preliminary selection : ",np.sum(preselection_mask)," / ",len(preselection_mask))

    count_all,binx_=np.histogram(dataStore['eleProbePt'][preselection_mask],bins=BINS)
    COUNT_ALL=[ int(i) for i in count_all]
    
    l1EGPt=dataStore['l1tPt']
    if args.doEmulated:
        print("Doing Emulated Branches ! ")
        if args.tag is None:
            args.tag='emulated'
        l1EGPt=dataStore['l1tEmuPt']
        if args.inputLUT is not None:
            print("Doing Isolation LUT Emulation ! ")
            intNTT=np.array(dataStore['l1tEmuNTT'],dtype=int)
            dataStore['compressedNTT']=utl.compressionNTTMap[intNTT]
            intRawEt=np.array(dataStore['l1tEmuRawEt'],dtype=int)
            intRawEt[intRawEt>255]=255
            intRawEt[intRawEt<0]=0
            dataStore['compressedRawEt']=utl.compressionRawEtMap[intRawEt]
            intIEta=np.abs(np.array(dataStore['l1tEmuTowerIEta'],dtype=int))
            intIEta[intIEta>31]=0
            dataStore['compressedIEta']=utl.compressionEtaMap[intIEta]
            
            ieta=dataStore['compressedIEta']
            iet =dataStore['compressedRawEt']
            intt=dataStore['compressedNTT']
            lut_index = ieta<<9 | iet <<5 | intt   
                
            isoLUT=utl.getIsoLUT(args.inputLUT)
            isolation_thr = isoLUT[lut_index]    
            isolation_bit = (dataStore['l1tEmuIsoEt'] < isolation_thr) & ( dataStore.l1tIso >=0)
        else:
            if args.isolation=='loose':
                args.tag+='_loose'
                isolation_bit = (dataStore['l1tIso'] ==2 ) | ( dataStore['l1tIso'] ==3 )
            if args.isolation=='tight':
                args.tag+='_tight'
                isolation_bit = (dataStore['l1tIso'] ==1 ) | ( dataStore['l1tIso'] ==3 )
    else:
        print("Doing Unpacked !")
        if args.tag is None:
            args.tag='unpacked'
        if args.isolation=='loose':
            print(" setting isolation  bit to loose")
            args.tag+='_loose'
            isolation_bit = (dataStore['l1tIso'] ==2 ) | ( dataStore['l1tIso'] ==3 )
        if args.isolation=='tight':
            print(" setting isolation  bit to tight")
            args.tag+='_tight'
            isolation_bit = (dataStore['l1tIso'] ==1 ) | ( dataStore['l1tIso'] ==3 )

    efficiencyMap={ }
    for et in range(81):
        effDict={}
        eff_name=f"SingleEG_{et}"
        print(f"\rProcessing {eff_name}   !",end="")
        triggerSeedMask = l1EGPt >= et
        triggerSeedMask = np.logical_and(triggerSeedMask,preselection_mask)
        
        count_trig,binx_= np.histogram(dataStore['eleProbePt'][triggerSeedMask],bins=BINS)
        effDict['pass'] = [int(i) for i in count_trig ]
        effDict['all']  = COUNT_ALL
        eff=count_trig/(count_all+1e-9)
        effDict['efficiency']  = [ float(i) for i in eff ]
        effDict['efficiency_err']  = [ 0.02*i for i in eff ]
        effDict['bin']  = BINS
        effDict['bin_center'] = BIN_X
        
        efficiencyMap[eff_name] = effDict
 
    for et in range(81):
        effDict={}
        eff_name=f"SingleEGIso_{et}"
        print(f"\rProcessing {eff_name}   !",end="")
        triggerSeedMask = l1EGPt >= et
        triggerSeedMask = np.logical_and(triggerSeedMask,preselection_mask)
        triggerSeedMask = np.logical_and(triggerSeedMask,isolation_bit)
   
        count_trig,binx_= np.histogram(dataStore['eleProbePt'][triggerSeedMask],bins=BINS)
        effDict['pass'] = [int(i) for i in count_trig ]
        effDict['all']  = COUNT_ALL
        eff=count_trig/(count_all+1e-9)
        effDict['efficiency']  = [float(i) for  i in eff]
        effDict['efficiency_err']  = [ float(0.02*i) for i in eff ]
        effDict['bin']  = BINS
        effDict['bin_center'] = BIN_X
   
        efficiencyMap[eff_name] = effDict

    ## Eff Plotting    

    f=plt.figure(figsize=(12,6))
    
    for i in range(10):
        print(f"\rPlotting {et} GeV Turnon",end="")
        et = 10+3*i
        eff_name=f"SingleEG_{et}"
        x=efficiencyMap[eff_name]['bin_center']
        y=efficiencyMap[eff_name]['efficiency']
        yerr=efficiencyMap[eff_name]['efficiency_err']
        plt.errorbar(x,y,yerr,
                     ecolor='r',
                     label=eff_name
                    )
    
        eff_name=f"SingleEGIso_{et}"
        x=efficiencyMap[eff_name]['bin_center']
        y=efficiencyMap[eff_name]['efficiency']
        yerr=efficiencyMap[eff_name]['efficiency_err']
        plt.errorbar(x,y,yerr,
                     ecolor='r',
                     label=eff_name
                    )
    
    
    plt.legend(fontsize=10,ncol=2)
    plt.xlim([-10,50])
    plt.ylabel("Efficiency")
    plt.xlabel("Offline E$_{T}$ [ GeV]")
    plt.hlines(1.0,10.0,100.0,color='k',alpha=0.5) 
    hep.cms.label("Internal",data=True,com=13.6,fontsize=15)
    print()

    foutname=f"{prefix}/{args.tag}_Efficiencies.png"
    print("Exporting ",foutname)
    f.savefig(foutname,bbox_inches='tight',dpi=150)
    foutname=f"{prefix}/{args.tag}_Efficiencies.pdf"
    print("Exporting ",foutname)
    f.savefig(foutname,bbox_inches='tight',dpi=150)
    
    ## Rate Export
    foutname=f"{prefix}/{args.tag}_Efficiencies.json"
    print("Exporting turnons to : ",foutname)
    with open(foutname,'w') as f:
        json.dump(efficiencyMap,f,indent=4)

if __name__ == "__main__":
    main()













