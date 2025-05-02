import uproot,json,os
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
from scipy import stats
import itertools,argparse
import util as utl

from util import *


# NOTE : iet is in compressed units , see the comression LUT here : https://github.com/ats2008/L1EGCalibrations/blob/regredssionUpd/Calibration/RegressionTraining/data/egCompressELUT_4bit_v4.txt 

##  For tight Isolation
def getEfficiency(ieta,iet):
    if iet <6.0:  # 6 --> 20 GeV
        return 0.7
    if iet <10.0: # 10--> 40 GeV
        return 0.7+0.3*(iet-6)/4
    return 1.0

##  For loose Isolation
def getEfficiency(ieta,iet):
    if iet <2.0:  # 2 --> 10 GeV
        return 0.7
    if iet <10.0: # 10--> 40 GeV
        return 0.7+0.3*(iet-2)/8

    return 1.0
    


hep.style.use("CMS")
plt.figure()
plt.close()
matplotlib.rcParams['figure.figsize'] = (6, 5.5)
matplotlib.rcParams['figure.dpi'] = 100

def main():

    parser = argparse.ArgumentParser()
    parser.add_argument('-i',"--inputFile", help="Input File",default='workarea/TagAndProbe_ReEmul_EGEraG_Remul_caloParam_v0_6_recalib_v0p2.root')
    parser.add_argument("-o","--dest", help="destination To Use", default='workarea/results/dataAugmentation/' )
    parser.add_argument("-p","--doSanityPlots", help="export the basic sanity plots",default=False,action='store_true')
    args = parser.parse_args()

    prefix=args.dest
    input_fname=args.inputFile
    
    # DEFINE THE BIN BOUNDARIES FOR DERIVATION OF FITS
    bin_splits={
        'compressedIEta' : [[0,3],[4,8],[9,12],[13,15]],
        'compressedRawEt' : [[i,i] for i in range(16)],
    }

    print("Opening the tag and probe file : ",input_fname)
    fileIn= uproot.open(input_fname)
    data=fileIn['Ntuplizer/TagAndProbe']
    dataStore=data.arrays(['l1tEmuNTT','l1tEmuRawEt','l1tEmuTowerIEta','l1tEmuEta','Nvtx','l1tEmuIsoEt'])
    
    dataMask=np.logical_and(dataStore['Nvtx'] > 20,dataStore['l1tEmuRawEt'] > 0)
    dataStore=dataStore[dataMask]
    cmd=f'mkdir -p {prefix}'
    print(cmd); os.system(cmd)
    
    ## ADDING EXTRA VARIABLES
    intNTT=np.array(dataStore['l1tEmuNTT'],dtype=int)
    dataStore['compressedNTT']=compressionNTTMap[intNTT]
    
    intRawEt=np.array(dataStore['l1tEmuRawEt'],dtype=int)
    intRawEt[intRawEt>255]=255
    dataStore['compressedRawEt']=compressionRawEtMap[intRawEt]
    
    intIEta=np.abs(np.array(dataStore['l1tEmuTowerIEta'],dtype=int))
    dataStore['compressedIEta']=compressionEtaMap[intIEta]
    
    saveFig=True
    
    if args.doSanityPlots:
        utl.makeSanityPlots(dataStore,prefix,saveFig)

    ## The bins here are inclusive to the boundary values provided. these variables need to be there in the data
    print("Making the dataset split in bins ")   
    print("The bin boundaries are : ")
    for ky in bin_splits:
        print(f"Bins along {ky} : {bin_splits[ky]}")
    bin_dim_names=list(bin_splits.keys())
    print()
    bins_bounds_=list(itertools.product(*[bin_splits[ky] for ky in bin_dim_names]))
    bins_bounds={getBinName(bin_dim_names,bb) : bb for bb in bins_bounds_}
    print(f"Number of bins = {len(bins_bounds)}")
    

    ## Splitting the data in bins , for derivation
    data_split_inBins={}
    nbins=len(bins_bounds)
    bin_mask=np.ones(len(dataStore['Nvtx']),dtype=bool)
    print("Splitting the datset into proposed bins ! ")
    for i,bn_name in enumerate(bins_bounds):
        print(f"\r  > Processing bin {i+1}/{nbins}",end="")
        bin_mask=dataStore['Nvtx'] > 0  ## A TRIVIAL MASK , all True
        for bn,bb in zip(bin_dim_names,bins_bounds[bn_name]):
            bin_mask=np.logical_and(bin_mask,dataStore[bn]>=bb[0])
            bin_mask=np.logical_and(bin_mask,dataStore[bn]<=bb[1])
        n=np.sum(bin_mask)
        if n < 1000:
            print("   > NEvents : ",n,"  (  < 1000 )")
            print()
        data_split_inBins[bn_name]=dataStore[bin_mask]
    
    print()    
    
    ## COMPRESSION SCHEME
    compressedEta_nbits=4
    compressedEt_nbits=4
    compressednTT_nbits=5
    nEtaBins=1<<compressedEta_nbits
    nNTTBins=1<<compressednTT_nbits
    nEtBins =1<<compressedEt_nbits
    
    ## DOING THE NTT vs IsoEt fits in each bins. The fits are done so that we can have a robust extrapolation to higher / lower NTT values
    k=-1
    extracted_data={}
    print("Making the Fits in each proposed bins !")
    for i,bin_name in enumerate(data_split_inBins):
        print(f"\r  > Processing bin {i+1}/{nbins}",end="")
        ldta=data_split_inBins[bin_name]
        x,y=np.array(ldta['compressedNTT']),np.array(ldta['l1tEmuIsoEt'])
        y[y>255]=255
    
        ## GET THE EFFICINCY FOR A GIVEN ieta and iet , NOTE : these are in compressed unit
        ieta=0.5*(bins_bounds[bin_name][0][0]+bins_bounds[bin_name][0][0])
        iet =0.5*(bins_bounds[bin_name][1][0]+bins_bounds[bin_name][1][0])
        qtile=getEfficiency(ieta,iet)*100 ## this function shold be customized according to need
        
        # getting the quantile of the IsoEt variable in each ntt bins 
    
        projection_default=getXProjectionsAtPercentiles( x,y,qtile)
        extracted_data[bin_name]={}
        extracted_data[bin_name]['x']=x
        extracted_data[bin_name]['y']=x
        extracted_data[bin_name]['c_ntt']  =projection_default['x']
        extracted_data[bin_name]['c_isoet']=projection_default['projection']
        extracted_data[bin_name]['c_isoet_std_u']=projection_default['stdev_u']
        extracted_data[bin_name]['c_isoet_std_d']=projection_default['stdev_d']
        
        ## Liner fit to the quantile points
        fit_range_mask = projection_default['x']  > -100 # TRIVIAL MASK
        fit_range_mask = np.logical_and(fit_range_mask , projection_default['x']  >= 0 ) 
        fit_range_mask = np.logical_and(fit_range_mask , projection_default['x']  < 5 ) 
    
        extracted_data[bin_name]['fit_mask']= fit_range_mask
        extracted_data[bin_name]['quantile']= qtile
        res = stats.linregress(projection_default['x'][fit_range_mask],
                               projection_default['projection'][fit_range_mask]
                              )
        
        extracted_data[bin_name]['linear_fit']={}
        extracted_data[bin_name]['linear_fit']["Remarks"] = "m*x+c"
        extracted_data[bin_name]['linear_fit']["c"] = res.intercept
        extracted_data[bin_name]['linear_fit']["m"] = res.slope
        ##  just setting upper IsoEt Cut 
        if qtile>99.9:
            extracted_data[bin_name]['linear_fit']["c"]=255
            extracted_data[bin_name]['linear_fit']["m"]=255
    

    print()    
    print("   Completed the fits !")
    
    if args.doSanityPlots:
        print("Exporting the fit validation plots ! ")
        axidx=-1
        f=None
        nPic=0
        for bin_name in extracted_data:
            if axidx < 0:
                if f is not None:
                    nPic+=1
                    foname=f"{prefix}/ntt_isoet_fitPlot_{nPic}.png"
                    print(f"Saving {foname}")
                    f.savefig(foname,bbox_inches='tight')
                f,axs=plt.subplots(4,4,figsize=(20,20),dpi=100)
                axs=np.ndarray.flatten(axs)
                axidx=0
            ax_=axs[axidx]
            _=ax_.hist2d(x,y,bins=[np.arange(0.0,32,2),np.arange(0.0,40,2)],norm='log')
            ax_.scatter(extracted_data[bin_name]['c_ntt'] , extracted_data[bin_name]['c_isoet'],color='magenta')
            mask=extracted_data[bin_name]['fit_mask']
            ax_.scatter(extracted_data[bin_name]['c_ntt'][mask] , extracted_data[bin_name]['c_isoet'][mask],color='green',marker='x')
            ax_.plot(np.arange(0,32),
                     np.arange(0,32)*extracted_data[bin_name]['linear_fit']["m"] + extracted_data[bin_name]['linear_fit']["c"] ,
                     color='r')
            ax_.set_xlim([0.0,32.0])
            ax_.set_ylim([0.0,40.0])
            txt=bin_name.replace('_compressed',' ').replace("_","").strip()
            quatile=extracted_data[bin_name]['quantile']
            txt+="\n"+f"Eff = {quatile}"
            ax_.text(15.0,30.0,txt,fontsize=10,bbox=dict(facecolor='cyan',alpha=0.2))
            axidx+=1
            if axidx==len(axs):
                axidx=-1
    
    
        if f is not None:
            nPic+=1
            foname=f"{prefix}/ntt_isoet_fitPlot_{nPic}.png"
            print(f"Saving {foname}")
            f.savefig(foname,bbox_inches='tight')            

    ## Evaluating the IsoEt thresholds from the fits in each bins : LUT
    print("Making the LUT ! ")
    LUT_MAP={'index':[],'threashold':[],'map':{}}
    for ieta in range(nEtaBins):
        if ieta not in LUT_MAP['map']:
            LUT_MAP['map'][ieta]={}
        for iet in range(nEtBins):
            print(f"\r Processing Bin {ieta=:>2} {iet=:>2}",end="")
            if iet not in LUT_MAP['map'][ieta]:
                LUT_MAP['map'][ieta][iet]={}
            
            bin_values={'compressedIEta':ieta,'compressedRawEt':iet}
            bin_bounds=getBinForValues(bin_values,bin_splits,bin_dim_names)
            bin_name=getBinName(bin_dim_names,bin_bounds)           
            edata=extracted_data[bin_name]
            print(bin_name,edata['linear_fit']['m'],edata['linear_fit']['c'])
            for intt in range(nNTTBins):
                index= (ieta<<9) | (iet<<5) | intt
                threashold=edata['linear_fit']['m']*intt+edata['linear_fit']['c']
                if threashold<0:
                    threashold=0
                if threashold>127:
                    threashold=127
                    
                threashold=int(threashold+1)
                if intt not in LUT_MAP['map'][ieta][iet]:
                    LUT_MAP['map'][ieta][iet][intt]={}
                LUT_MAP['map'][ieta][iet][intt]['threashold']=threashold
                LUT_MAP['map'][ieta][iet][intt]['index']     =index
                
                LUT_MAP['index'].append(index)
                LUT_MAP['threashold'].append(threashold)
    
    print()    
    print("   LUT is made !")
    
    ## Validation plot for the evaluated fits
    print("Making the Threashold vs NTT plots")       
    f,ax=plt.subplots(4,4,figsize=(16,16),dpi=150)
    ax=np.ndarray.flatten(ax)
    for ietas in bin_splits['compressedIEta']:
        ieta=ietas[0]
        for iet in LUT_MAP['map'][ieta]:
            x=[]
            y=[]
            for intt in LUT_MAP['map'][ieta][iet]:
                x.append(intt)
                thr=LUT_MAP['map'][ieta][iet][intt]['threashold']
                if thr>39: thr=39
                y.append(thr)
            ax[iet].scatter(x,y,s=5,label=f"ieta={ieta}")
    
    for ieta in LUT_MAP['map']:
        for iet in LUT_MAP['map'][ieta]:
            ax[iet].legend(ncol=2,fontsize=12)
            ax[iet].set_ylim([0,42])
            ax[iet].text(4,20,f"comp. iET : {iet}",fontsize=12)
            
            if (iet%4)==0:
                ax[iet].set_ylabel("IsoEt thr.")
            if  iet>11:
                ax[iet].set_xlabel("Compressed NTT")
        break
    fname_out=f"{prefix}/thresholdVsNTT.png"
    #if saveFig:
    print(f"   exporting {fname_out}")
    f.savefig(fname_out,bbox_inches='tight',dpi=150)
    print()    
    
    ## Making the standard LUT pictures 
    print("Making the LUT plots")       
    LUT=np.zeros((16,16,32),int)
    for i in LUT_MAP['map']:
        for j in LUT_MAP['map'][i]:
            for k in LUT_MAP['map'][i][j]:
                LUT[i][j][k]=LUT_MAP['map'][i][j][k]['threashold']
    
    utl.plotLUT(LUT,prefix)
    
    # Exporting the LUT into CMSSW compacent format
    print("Exporting the LUT !")
    k=0
    LUT_TEXT=""
    LUT_TEXT+="# EG isolation LUT with WP = ISOLATION WP\n"
    LUT_TEXT+="# iso LUT structure is ieta --> iEt -->  nTT\n"
    LUT_TEXT+="# Compr bits: ieta: 4 iEt: 4 nTT: 5\n"
    LUT_TEXT+="#<header> V10.0 13 9 </header>\n"
    LUT_TEXT+="\n"
    for ieta in range(nEtaBins):
        for iet in range(nEtBins):
            for intt in range(nNTTBins):
                index=LUT_MAP['index'][k]
                thr=LUT_MAP['threashold'][k]
                
                ostr=f"{index} {thr} # ieta : iEt : inTT = {ieta} : {iet} : {intt}"
                LUT_TEXT+=ostr
                LUT_TEXT+="\n"
                k+=1
    foutname=f"{prefix}/LUT.txt"
    with open(foutname,'w') as f:
        print("LUT exported to ",foutname)
        f.write(LUT_TEXT)

if __name__ == "__main__":
    main()

