//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Dec  1 07:58:18 2021 by ROOT version 6.24/06
// from TTree L1UpgradeTree/L1UpgradeTree
// found on file: L1Ntuple.root
//////////////////////////////////////////////////////////

#ifndef ApplyIsolation_h
#define ApplyIsolation_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <map>
#include <TFile.h>
#include <TH1F.h>
#include <TH1D.h>
#include <TH3F.h>
#include <vector>
#include <map>
#include <TGraphAsymmErrors.h>


// Header file for the classes stored in the TTree if any.
#define nBins_fine 100
#define nBins 38
class ApplyIsolation {
public :
    TChain          *fChain;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent; //!current Tree number in a TChain

    TChain          *fChain_1;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent_1; //!current Tree number in a TChain


    TChain          *fChain1;   //!pointer to the analyzed TTree or TChain
    Int_t           fCurrent1; //!current Tree number in a TChain


    // Declaration of leaf types

    //Rate
    UShort_t        nEGs;
    Int_t           nPV_True;
    UInt_t run;
    std::vector<float>   egEt;
    std::vector<float>   egEta;
    std::vector<short>   egBx;
    std::vector<short>   egIsoEt;
    std::vector<short>   egNTT;
    std::vector<short> egRawEt;
    std::vector<short> egTowerIEta;


    //Turn-On
    Float_t           l1tEmuPt;
    Float_t           l1tEmuTowerIEta;
    Float_t           l1tEmuNTT;
    Float_t           l1tEmuRawEt;
    Float_t           eleProbeSclEt;
    Float_t           l1tEmuIsoEt;
    Float_t           l1tEmuIso;

    Int_t   hasL1Emu_tightiso22;
    Int_t   hasL1Emu_tightiso24;
    Int_t   hasL1Emu_tightiso26;
    Int_t   hasL1Emu_looseiso22;
    Int_t   hasL1Emu_looseiso24;
    Int_t   hasL1Emu_looseiso26;
    Int_t   hasL1Emu_22;
    Int_t   hasL1Emu_24;
    Int_t   hasL1Emu_26;

    Float_t           Nvtx;

    Float_t   isProbeLoose;
    Float_t   eleTagPhi;
    Float_t   eleTagEta;
    Float_t   eleProbePhi;
    Float_t   eleProbeEta;


    ApplyIsolation(std::string& inputFileName);
    ~ApplyIsolation();

    void accessTree(std::string & input_filelist_rate, std::string & input_filelist_turn_on);
    void readParameters(const std::string jfile);
    void readTree();
    void loops();
    void bookHistograms();
    void saveHistograms();
    void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
    void scaleHistograms();
    void createEfficiencyHistograms(TH1* h_P, TH1* h_F, TH1* h_E);
    double BetaInverse(double x,double p, double q);
    double FindRate(std::string it, UInt_t e);
    double FindRateD(std::string it, UInt_t e);

    short getEtaBin(short eta);
    short getEtBin(short et);
    short getnTTBin(short ntt);


    Long64_t nEntries_;
    Long64_t nEntries1_;
    Long64_t maxEntriesForEfficiency_;
    Long64_t maxEntriesForRate_;
    Long64_t reportEvery_;
    UInt_t etMin_;
    UInt_t etMax_;
    float  Et_Double;
    std::string ntupleFileNameRate_;
    std::string ntupleFileNameTurnOn_;
    std::string outputFileName_;
    std::string weightFileName_;
    std::string optionsFileName_;
    std::vector<std::string> lutProgOptVec_;


    TFile* optionsFile_;
    TFile* outputFile_;
    TFile* weightFile_;
    bool bookedHistograms_;
    bool nTTRange_;

    void readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map);

    UInt_t nBinsIEta;
    UInt_t nBinsIEt;
    UInt_t nBinsNTT;
    UInt_t ET_MAX = 255;
    UInt_t bunches =2540;

    std::map<short, short> lutMapEt;
    std::map<short, short> lutMapEta;
    std::map<short, short> lutMapNTT;

    std::map<TString,TH1F*> ptMap_;
    std::map<TString,TH1F*> rateMap_;
    std::map<TString,TH1F*> pt_pass_Map_;
    std::map<TString,TGraphAsymmErrors*> turnOn_Map_;
    std::map<TString,TH1F*> pt_pass_Map_Run2_;
    std::map<TString,TGraphAsymmErrors*> turnOn_Map_Run2_;
    std::map<TString,TH1F*> th1fStore;
    std::map<TString,TH1D*>  Nvtx_pass_Map_;
    std::map<TString,TH1D*> Nvtx_fail_Map_;
    std::map<TString,TH1D*>  Nvtx_Eff_Map_;
    std::map<TString,TH1D*>  eta_pass_Map_;
    std::map<TString,TH1D*> eta_fail_Map_;
    std::map<TString,TH1D*>  eta_Eff_Map_;
    std::map<TString,TH1D*>  nTT_pass_Map_;
    std::map<TString,TH1D*> nTT_fail_Map_;
    std::map<TString,TH1D*>  nTT_Eff_Map_;

    //   int etMax_;
    bool doRateOptimization_;
    std::vector<std::string>  histoLabelVec_;
    std::vector<TH3F*>optionsHistoVec_;
    short in_compressediEta;
    short in_compressediEt;
    short in_compressedNTT;


    TGraphAsymmErrors*   turnon_inclusive;
    TH1F* PtPass_inclusive;
    std::vector<UInt_t> et_option;
    std::vector<UInt_t> et_optionD;

    TH1F* pT_all;
    std::vector<UInt_t> Et_fr;


    bool check_pt_rate_dir = false;
    bool check_rate_dir = false;
    bool check_pt_turn_on_dir = false;
    bool check_turn_on_dir = false;
    bool check_turn_on_double_dir = false;
    bool check_pt_rate_double_dir =false;
    bool check_rate_double_dir =false;
    bool check_turn_on_dir_Run2 =false;
    bool check_pt_turn_on_fr_dir =false;
    bool check_pt_turn_on_fr_double_dir =false ;
    bool check_turn_on_fr_dir =false;
    bool check_turn_on_fr_double_dir =false ;
    bool check_Nvtx =false;
    bool check_eta =false;
    bool check_nTT=false;
    bool check_Nvtx_Iso =false;
    bool check_eta_Iso =false;
    bool check_nTT_Iso = false;
    bool check_pt_rateI_dir =false;
    bool check_pt_rateI_double_dir =false;
    bool check_rateI_dir =false;
    bool check_rateI_double_dir =false;
    TDirectory* td;
    TDirectory* td1;
    TDirectory* tdD;
    TDirectory* td1D;
    TDirectory* td2;
    TDirectory* td3;
    TDirectory* td4;
    TDirectory* td7;
    TDirectory* td8;
    TDirectory* td5;
    TDirectory* td6;
    TDirectory* td9;
    TDirectory* td10;
    TDirectory* td11;
    TDirectory* td12;
    TDirectory* td13;
    Double_t xEdges[39] = {1., 3., 5., 7., 9.,  10., 12., 15., 18., 20., 22., 24., 26., 28.,
                            29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41.,
                            42., 43., 45., 50., 60., 70., 100., 200., 300., 400., 600., 1000.
                           };


    Double_t xEdges_fine[nBins_fine+1];




};

#endif












