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

def singleIEtSeed(data,iet):
    seed_cond=  data['egEt'] >= iet
    return ak.sum(seed_cond,axis=1) > 0

def singleIEtIsoSeed(data,iet):
    seed_cond = data['is_iso'] & (data['egEt'] >= iet)
    return ak.sum(seed_cond,axis=1) > 0


def doubleIEtSeed(data,iet1,iet2):
    leg1 = data['egEt'] >= iet1
    leg2 = data['egEt'] >= iet2
    nSeedPass= ( ak.sum(leg1,axis=1) > 0) & ( ak.sum(leg2,axis=1) > 1)
    return nSeedPass

def doubleIEtSingleIsoSeed(data,iet1,iet2):
    leg1 = (data['egEt'] >= iet1) & data['is_iso']
    leg2 = (data['egEt'] >= iet2) 
    nSeedPass= ( ak.sum(leg1,axis=1) > 0) & ( ak.sum(leg2,axis=1) > 1)
    return nSeedPass

def doubleIEtDoubelIsoSeed(data,iet1,iet2):
    leg1 = (data['egEt'] >= iet1) & data['is_iso']
    leg2 = (data['egEt'] >= iet2) & data['is_iso']
    nSeedPass= ( ak.sum(leg1,axis=1) > 0) & ( ak.sum(leg2,axis=1) > 1)
    return nSeedPass

selectorStore={}
selectorStore['SingleEG_Inclusive'] = singleIEtSeed
selectorStore['SingleIsoEG'] = singleIEtIsoSeed
selectorStore['DoubleEG_Inclusive'] = lambda d,e1: doubleIEtSeed(d,e1,e1)
selectorStore['DoubleIsoEG'] = lambda d,e1: doubleIEtDoubelIsoSeed(d,e1,e1) 
selectorStore['DoubleEGXandXp10_SingleIsoXplus10'] = lambda d,e1: doubleIEtSingleIsoSeed(d,e1+10,e1)

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
    isolationFname=args.inputLUT
    
    print("Input Ntuple : ",fname)
    print("Input LUT : ",isolationFname)
    print("Destination : ",prefix)
    
    cmd=f'mkdir -p {prefix}'
    print(cmd); os.system(cmd)
    
    isoLUT=utl.getIsoLUT(isolationFname)
    doEmulated=args.doEmulated
    
    fileIn= uproot.open(fname)
    if doEmulated:
        tree=fileIn['l1UpgradeEmuTree/L1UpgradeTree/L1Upgrade']
    else:
        tree=fileIn['l1UpgradeTree/L1UpgradeTree/L1Upgrade']
    
    branches=['egEt', 'egEta', 'egPhi', 'egIEt', 'egIEta', 'egIPhi', 'egIso', 
              'egBx', 'egTowerIPhi', 'egTowerIEta', 'egRawEt', 'egIsoEt', 
              'egFootprintEt', 'egNTT', 'egShape', 'egTowerHoE', 'egHwQual']
    
    data=tree.arrays(branches)
    N_EVTS=len(data)

    if doEmulated:  
        print("Processing the new LUT bits")
        if args.tag is None:
            args.tag='exported'
        dataStore={}
        
        counts=ak.num(data.egRawEt)
        
        ## DEFINING COMPRESSED VARIABLES and flttened object vars
        
        print("Making the compressed variables")

        intRawEt = np.array(ak.flatten( data.egRawEt))
        intRawEt[ intRawEt > 255 ] =255
        dataStore['l1tEmuRawEt'] = intRawEt
        dataStore['compressedRawEt']=utl.compressionRawEtMap[intRawEt]
        data['compressedRawEt']=ak.unflatten(dataStore['compressedRawEt'],counts)
        
        counts=ak.num(data.egNTT)
        
        intNTT = np.array(ak.flatten( data.egNTT))
        intNTT[ intNTT > 1023 ] =1023
        dataStore['counts']=counts
        dataStore['l1tEmuNTT']=intNTT
        dataStore['compressedNTT']=utl.compressionNTTMap[intNTT]
        data['compressedNTT']=ak.unflatten(dataStore['compressedNTT'],counts)
        
        counts=ak.num(data.egIEta)
        
        intIEta = np.array(ak.flatten( data.egIEta))
        dataStore['l1tEmuTowerIEta']=intIEta
        cbis=np.abs(intIEta)
        cbis[cbis>31]=31
        dataStore['compressedIEta']=utl.compressionEtaMap[cbis]
        data['compressedIEta']=ak.unflatten(dataStore['compressedIEta'],counts)
        
        dataStore['l1tEmuIsoEt']=ak.flatten(data['egIsoEt'])
        
        saveFig=True
        
        if args.doSanityPlots:
            print("MAKE")
            utl.makeSanityPlotsDicts(dataStore,prefix,saveFig)
        
        print("Making the LUT index")
        ieta=dataStore['compressedIEta']
        iet =dataStore['compressedRawEt']
        intt=dataStore['compressedNTT']
        lut_index = ieta<<9 | iet <<5 | intt
        
        print("Making the Iso Bit")
        isolation_thr = isoLUT[lut_index]
        dataStore['isolation_bit'] = dataStore['l1tEmuIsoEt'] < isolation_thr 
        data['is_iso'] = ak.unflatten(dataStore['isolation_bit'],dataStore['counts'])

    else:
        if args.tag is None:
            args.tag='unpacked'

        if args.iso == 'loose':
            print("Processing unpacked loose isolation")
            data['is_iso'] = (data.egIso==2) | (data.egIso== 3)
            args.tag+='_loose'
        elif args.isolation == 'tight':
            print("Processing unpacked tight isolation")
            data['is_iso'] = (data.egIso==1) | (data.egIso== 3)
            args.tag+='_tight'
        

    
    ## RATE EVALUATION ! 
    rate_scale_factor= BPTX_RATE*1e3/N_EVTS

    rateStore={}
    for selector_ky in selectorStore:
        print("Making rates for  ",selector_ky)
        selector=selectorStore[selector_ky]
        rateDict={'threashold':[],'rate':[],'rate_err':[],'pass':[],'n_events':int(N_EVTS)}
        et_thr=[]
        rates=[]
        for i in range(81):
            print(f"\rDoing Et = {i}   ",end="")
            seed_events = selector(data,i)
            pass_count=int(ak.sum(seed_events))
            rateDict['pass'].append(pass_count)
            rateDict['rate'].append(pass_count*rate_scale_factor)
            rateDict['rate_err'].append(pass_count**0.5*rate_scale_factor)
            rateDict['threashold'].append(i)
        rateStore[selector_ky]=rateDict
        print()
    
    ## Rate Plotting    
    f=plt.figure(figsize=(6,5.5))
    for ky in rateStore:
        print("Plotting : ",ky)
        plt.errorbar(rateStore[ky]['threashold'],
                     rateStore[ky]['rate'],
                     rateStore[ky]['rate_err'],
                     ecolor='r',
                     label=ky
                    )
    plt.semilogy()
    plt.legend(fontsize=10)
    plt.ylim([5e-1,1e5])
    plt.axhline(18.0,color='k',alpha=0.4)
    plt.text(50.0,25,"18 kHz",fontsize=10)
    plt.axhline(4.0,color='k',alpha=0.4)
    plt.text(3.0,5,"4 kHz",fontsize=9)
    plt.ylabel('L1 Rate [ kHz ]')
    plt.xlabel('L1 E$_{T}$ [ GeV ]')
    hep.cms.label('Internal',data=True,fontsize=15,com=13.6)
    
    foutname=f"{prefix}/{args.tag}_Rates.png"
    print("Exporting ",foutname)
    f.savefig(foutname,bbox_inches='tight',dpi=150)
    foutname=f"{prefix}/{args.tag}_Rates.pdf"
    print("Exporting ",foutname)
    f.savefig(foutname,bbox_inches='tight',dpi=150)
    
    ## Rate Export
    foutname=f"{prefix}/{args.tag}_Rates.json"
    print("Exporting rate to : ",foutname)
    with open(foutname,'w') as f:
        json.dump(rateStore,f,indent=4)

if __name__ == "__main__":
    main()
