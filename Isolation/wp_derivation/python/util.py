import uproot,json,os
import ROOT
import numpy as np
import matplotlib
import matplotlib.pyplot as plt
import mplhep as hep
from scipy import stats
import itertools

def getBinName(bin_dim_names,bin_bounds):
    name="_"
    for bn,bb in zip(bin_dim_names,bin_bounds):
        name+=f"{bn}IN{bb}_"
    return name.replace(" ","")

def getBinForValues(bin_values,bin_splits,bin_dim_names):
    bins=[]
    for i in bin_dim_names:
        if i not in bin_splits:
            print(f"Bin splits does not have {i}")
            raise KeyError
        if i not in bin_values:
            print(f"Bin values does not have {i}")
            raise KeyError
        v=bin_values[i]
        gotBin=False
        for bn in bin_splits[i]:
            if (v<=bn[1]) and (v>=bn[0]):
                bins.append(bn)
                gotBin=True
        if not gotBin:
            print(f"For dim {i}, value {v} did not find a bin match !")
            raise KeyError
    return bins
            
def getXProjectionsAtPercentiles(x_arr,y_arr,quantile=50,nMinProj=20):
    xU = np.unique(x_arr)
    xVal     =[]     
    xProj    =[]
    xProj_err_up=[]
    xProj_err_dn=[]
    for x in xU:
        mask= x_arr==x
        if np.sum(mask) < nMinProj:
            continue
        val=np.percentile(y_arr[mask],quantile)
        n=np.sum(mask)
        var=y_arr[mask]-val
        varU=var[var>=0]
        varD=var[var<0]
        std_u=np.sqrt(np.sum(varU*varU)/n)
        std_d=np.sqrt(np.sum(varD*varD)/n)
        xProj_err_up.append(std_u)
        xProj_err_dn.append(std_d)
        
        xProj.append(val)
        xVal.append(x)
        
    return {'x':np.array(xVal) , 'projection' : np.array(xProj) ,
            'stdev_u' : np.array(xProj_err_up) ,'stdev_d' : np.array(xProj_err_dn) }




def getXProjections(x_arr,y_arr,nMinProj=-1):
    xU = np.unique(x_arr)
    xVal     =[]     
    xProj    =[]
    xProj_err=[]
    for x in xU:
        mask= x_arr==x
        if np.sum(mask) < nMinProj:
            continue
        avg = np.average( y_arr[mask] )
        std = np.std( y_arr[mask] )
        xProj.append(avg)
        xProj_err.append(std)
        xVal.append(x)
    return {'x':np.array(xVal) , 'projection' : np.array(xProj) , 'stdev' : np.array(xProj_err) }



### OBJECTS
compressionNTTMap=np.array([0,0,0,0,0,0,1,1,1,1,1,2,2,2,2,2,3,3,3,3,3,
                         4,4,4,4,4,5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,8
                         ,8,8,8,8,9,9,9,9,9,10,10,10,10,10,
                         11,11,11,11,11,12,12,12,12,12,13,13,13,13,13,
                         14,14,14,14,14,15,15,15,15,15,16,16,16,16,16,
                         17,17,17,17,17,18,18,18,18,18,19,19,19,19,19,
                         20,20,20,20,20,21,21,21,21,21,22,22,22,22,22,
                         23,23,23,23,23,24,24,24,24,24,25,25,25,25,25,
                         26,26,26,26,26,27,27,27,27,27,28,28,28,28,28,
                         29,29,29,29,29,30,30,30,30,30,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,
                         31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,31,])

compressionRawEtMap=np.array([0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,0,1,1,2,2,3,3,3,3,3,3,4,4,4,4,
                              5,5,5,5,5,6,6,6,6,6,7,7,7,7,7,7,7,7,7,7,8,8,8,8,8,8,8,8,8,8,8,
                              9,9,9,9,9,9,9,9,9,9,10,10,10,10,10,10,10,10,11,11,11,11,11,11,
                              12,12,12,12,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,13,
                              13,13,13,13,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
                              14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,14,
                              14,14,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,
                              15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15,15]
                            )


compressionEtaMap=np.array([0,0,0,0,0,
                            1,
                            2,2,2,
                            3,
                            4,4,
                            5,6,
                            7,7,7,
                            8,9,10,
                            11,11,11,
                            12,13,14,
                            15,15,15,
                            15,15,15])


def makeSanityPlots(dataStore,prefix='',saveFig=False):
    print("Making the sanity check plots  !")
    if 'Nvtx' in dataStore.fields:
        print("  >> N-Reco-Vertex plot")
        plt.figure()
        _=plt.hist(dataStore['Nvtx'],np.arange(-0.5,80.0,1),alpha=0.6)
        mu=np.average(dataStore['Nvtx'])
        stdev=np.std(dataStore['Nvtx'])
        plt.text(4,max(_[0])*0.5,f"<N Vertex> : {mu:.1f}\n RMS<N Vertex> {stdev:.1f}",fontsize=13)
        plt.xlabel('N Vertex',fontsize=14)
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024')
        if saveFig:
            fname_out=f'{prefix}/nvertex.png'
            print(f"   exporting {fname_out}")
            plt.savefig(fname_out,bbox_inches='tight')
        
    if 'l1tEmuIsoEt' in dataStore.fields:
        print("  >> IsoEt  plot")
        plt.figure()
        _=plt.hist(dataStore['l1tEmuIsoEt'],np.arange(-0.5,40.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['l1tEmuIsoEt'])
        stdev=np.std(dataStore['l1tEmuIsoEt'])
        plt.text(20,max(_[0])*0.7,f"<IsoEt> : {mu:.1f}\n RMS<IsoEt> {stdev:.1f}",fontsize=13)
        plt.xlabel('IsoEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024')
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        if saveFig:
            fname_out=f'{prefix}/IsoEt.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/IsoEt.png',bbox_inches='tight')
        
    if ('l1tEmuRawEt' in dataStore.fields) and ('compressedRawEt' in dataStore.fields):
        print("  >> Raw Et plot")
        f,ax=plt.subplots(1,2,figsize=(12,4))
        _=ax[0].hist(dataStore['l1tEmuRawEt'],np.arange(-0.5,180.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['l1tEmuRawEt'])
        stdev=np.std(dataStore['l1tEmuRawEt'])
        ax[0].text(110,max(_[0])*0.7,f"<RawEt> : {mu:.1f}\n RMS<RawEt> {stdev:.1f}",fontsize=13)
        ax[0].set_xlabel('RawEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        
        _=ax[1].hist(dataStore['compressedRawEt'],np.arange(-0.5,16.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['compressedRawEt'])
        stdev=np.std(dataStore['compressedRawEt'])
        ax[1].text(0,max(_[0])*0.7,f"<CompRawEt> : {mu:.1f}\n RMS<CompRawEt> {stdev:.1f}",fontsize=13)
        ax[1].set_xlabel('compressedRawEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        if saveFig:
            fname_out=f'{prefix}/compressedRawEt.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/compressedRawEt.png',bbox_inches='tight')    
        
    if ('l1tEmuTowerIEta' in dataStore.fields) and ('compressedIEta' in dataStore.fields):
        print("  >> IEta plot")
        f,ax=plt.subplots(2,1,figsize=(16,12))
        _=ax[0].hist(dataStore['l1tEmuTowerIEta'],np.arange(-30.5,31.0,1),alpha=0.6,color='g')
        # mu=np.average(dataStore['l1tEmuEta'])
        # stdev=np.std(dataStore['l1tEmuEta'])
        # plt.text(110,max(_[0])*0.7,f"<l1tEmuEta> : {mu:.1f}\n RMS<l1tEmuEta> {stdev:.1f}",fontsize=13)
        ax[0].set_xlabel('TowerIEta',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        
        _=ax[1].hist(dataStore['compressedIEta'],np.arange(-0.5,17.0,1),alpha=0.6,color='g')
        # mu=np.average(dataStore['l1tEmuEta'])
        # stdev=np.std(dataStore['l1tEmuEta'])
        # plt.text(110,max(_[0])*0.7,f"<l1tEmuEta> : {mu:.1f}\n RMS<l1tEmuEta> {stdev:.1f}",fontsize=13)
        ax[1].set_xlabel('compressedIEta',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[1])
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        if saveFig:
            fname_out=f'{prefix}/compressedIEta.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/compressedIEta.png',bbox_inches='tight')

def plotLUT(LUT,prefix='./'):
    f,axlist=plt.subplots(8,2,figsize=(10,18))
    axlist=np.ndarray.flatten(axlist)
    for et in range(16):
        ax=axlist[et]
        lut_proj=LUT[:,et,:]
    #     f,ax=plt.subplots(figsize=(32/4,16/4))
        maxIsoEt=40
        if et > 9:
            maxIsoEt=120
        lut_proj[lut_proj>maxIsoEt]=maxIsoEt
        c=ax.imshow(lut_proj,cmap='tab20c',origin='lower',vmin=0, vmax=maxIsoEt)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4),minor=True,fontsize=12)
        t=ax.set_yticks(np.arange(0.0,16,2)-0.5,np.arange(0.0,16,2),minor=True,fontsize=12)
        t=ax.set_xticks(np.arange(0.0,32,4)-0.5,np.arange(0.0,32,4),fontsize=10)
        t=ax.set_yticks(np.arange(0.0,16,4)-0.5,np.arange(0.0,16,4),color='w')
        ax2 = f.add_axes([ax.get_position().x1+0.005,ax.get_position().y0,0.02,ax.get_position().height])
        ax.grid(color='k',alpha=1,which='both')
        t=ax.text(25.2,12.3,"E$_{T}$ = "+str(et) ,fontweight='bold',c='r',alpha=1,fontsize=10)
        t.set_bbox(dict(facecolor='w',alpha=1.0, edgecolor='r'))
        cbar=plt.colorbar(c,cax=ax2)
        cbar.ax.tick_params(labelsize=10)
        if (et%2)==0:
            ax.set_ylabel("IsoEt thr.",fontsize=15)
        if et>13:
            ax.set_xlabel("comp. NTT",fontsize=15)
    
    fname_out=f"{prefix}/isoEtLUTs.png"
    #if saveFig:
    print(f"   exporting {fname_out}")
    f.savefig(fname_out,bbox_inches='tight',dpi=150)



def makeSanityPlotsDicts(dataStore,prefix='',saveFig=False):
    print("Making the sanity check plots  !")
    if 'Nvtx' in dataStore:
        print("  >> N-Reco-Vertex plot")
        plt.figure()
        _=plt.hist(dataStore['Nvtx'],np.arange(-0.5,80.0,1),alpha=0.6)
        mu=np.average(dataStore['Nvtx'])
        stdev=np.std(dataStore['Nvtx'])
        plt.text(4,max(_[0])*0.5,f"<N Vertex> : {mu:.1f}\n RMS<N Vertex> {stdev:.1f}",fontsize=13)
        plt.xlabel('N Vertex',fontsize=14)
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024')
        if saveFig:
            fname_out=f'{prefix}/nvertex.png'
            print(f"   exporting {fname_out}")
            plt.savefig(fname_out,bbox_inches='tight')
        
    if 'l1tEmuIsoEt' in dataStore:
        print("  >> IsoEt  plot")
        plt.figure()
        _=plt.hist(dataStore['l1tEmuIsoEt'],np.arange(-0.5,40.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['l1tEmuIsoEt'])
        stdev=np.std(dataStore['l1tEmuIsoEt'])
        plt.text(20,max(_[0])*0.7,f"<IsoEt> : {mu:.1f}\n RMS<IsoEt> {stdev:.1f}",fontsize=13)
        plt.xlabel('IsoEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024')
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        #plt.semilogy()
        if saveFig:
            fname_out=f'{prefix}/IsoEt.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/IsoEt.png',bbox_inches='tight')
        
    if ('l1tEmuRawEt' in dataStore) and ('compressedRawEt' in dataStore):
        print("  >> Raw Et plot")
        f,ax=plt.subplots(1,2,figsize=(12,4))
        _=ax[0].hist(dataStore['l1tEmuRawEt'],np.arange(-0.5,180.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['l1tEmuRawEt'])
        stdev=np.std(dataStore['l1tEmuRawEt'])
        ax[0].text(110,max(_[0])*0.01,f"<RawEt> : {mu:.1f}\n RMS<RawEt> {stdev:.1f}",fontsize=13)
        ax[0].set_xlabel('RawEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        ax[0].semilogy()
        
        _=ax[1].hist(dataStore['compressedRawEt'],np.arange(-0.5,16.0,1),alpha=0.6,color='g')
        mu=np.average(dataStore['compressedRawEt'])
        stdev=np.std(dataStore['compressedRawEt'])
        ax[1].text(4,max(_[0])*0.01,f"<CompRawEt> : {mu:.1f}\n RMS<CompRawEt> {stdev:.1f}",fontsize=13)
        ax[1].set_xlabel('compressedRawEt',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        ax[1].semilogy()
        if saveFig:
            fname_out=f'{prefix}/compressedRawEt.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/compressedRawEt.png',bbox_inches='tight')    
        
    if ('l1tEmuTowerIEta' in dataStore) and ('compressedIEta' in dataStore):
        print("  >> IEta plot")
        f,ax=plt.subplots(2,1,figsize=(16,12))
        _=ax[0].hist(dataStore['l1tEmuTowerIEta'],np.arange(-30.5,31.0,1),alpha=0.6,color='g')
        # mu=np.average(dataStore['l1tEmuEta'])
        # stdev=np.std(dataStore['l1tEmuEta'])
        # plt.text(110,max(_[0])*0.7,f"<l1tEmuEta> : {mu:.1f}\n RMS<l1tEmuEta> {stdev:.1f}",fontsize=13)
        ax[0].set_xlabel('TowerIEta',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[0])
        
        _=ax[1].hist(dataStore['compressedIEta'],np.arange(-0.5,17.0,1),alpha=0.6,color='g')
        # mu=np.average(dataStore['l1tEmuEta'])
        # stdev=np.std(dataStore['l1tEmuEta'])
        # plt.text(110,max(_[0])*0.7,f"<l1tEmuEta> : {mu:.1f}\n RMS<l1tEmuEta> {stdev:.1f}",fontsize=13)
        ax[1].set_xlabel('compressedIEta',fontsize=14)
        hep.cms.label('Internal',fontsize=12,data=False,year='2024',ax=ax[1])
        # hep.cms.label('Internal',fontsize=18,data=True,year='2023 C',lumi=9.7)
        if saveFig:
            fname_out=f'{prefix}/compressedIEta.png'
            print(f"   exporting {fname_out}")
            plt.savefig(f'{prefix}/compressedIEta.png',bbox_inches='tight')




def exportIsoLUT():
    pass

def getIsoLUT(isolationFname):
    with open(isolationFname) as f:
        txt=f.readlines()
    idx=0
    isoLUT=[]
    for l in txt:
        if len(l) <2 : continue
        if l.startswith('#') : continue
        l=l.split('#')[0].strip().split()
        index=int(l[0])
        if index!=idx:
            print(f"LUT not proper ! for expected index {idx} , the value obtained is {index}")
            raise ValueError
        isoLUT.append(int(l[1]))
        idx+=1
    isoLUT=np.array(isoLUT)        
    return isoLUT

