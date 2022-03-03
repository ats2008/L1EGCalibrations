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

class ApplyIsolation {
 public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  TChain          *fChain1;   //!pointer to the analyzed TTree or TChain                                             
  Int_t           fCurrent1; //!current Tree number in a TChain                                                                                                                                            
  
  
  // Declaration of leaf types
  
  //Rate
  UShort_t        nEGs;
  std::vector<float>   egEt;
  std::vector<short>   egBx;
  std::vector<short>   egIsoEt;
  std::vector<short>   egNTT;
  std::vector<short> egRawEt;  
  std::vector<short> egTowerIEta;


  //Turn-On
   Float_t         l1tEmuPt;
   Int_t           l1tEmuTowerIEta;
   Int_t           l1tEmuNTT;
   Int_t           l1tEmuRawEt;    
   Float_t         eleProbeSclEt;
   Int_t           l1tEmuIsoEt;

   Int_t   hasL1Emu_tightiso24;
   Int_t   hasL1Emu_tightiso26;
   Int_t   hasL1Emu_looseiso24;
   Int_t   hasL1Emu_looseiso26;
   Int_t   hasL1Emu_24;
   Int_t   hasL1Emu_26;
   
   
   Int_t   isProbeLoose;
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
   
   
   Long64_t nEntries_;
   Long64_t nEntries1_;
   Long64_t maxEntriesForEfficiency_;
   Long64_t maxEntriesForRate_;
   Long64_t reportEvery_;
   float frate_;
   
   std::string ntupleFileNameRate_;
   std::string ntupleFileNameTurnOn_;
   std::string outputFileName_;
   std::string optionsFileName_;    
   std::vector<std::string> lutProgOptVec_;

   
   TFile* optionsFile_;
   TFile* outputFile_;

   bool bookedHistograms_;
   bool nTTRange_;

   void readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map);

   UInt_t nBinsIEta;
   UInt_t nBinsIEt;
   UInt_t nBinsNTT;
   UInt_t ET_MAX = 255;
   UInt_t bunches =2544;
   
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

   TH1F* pT_all;
   std::vector<UInt_t> et_option;

   bool check_pt_rate_dir = false;
   bool check_rate_dir = false;
   bool check_pt_turn_on_dir = false;
   bool check_turn_on_dir = false;                                                                                             
   bool check_turn_on_dir_Run2 =false;
   TDirectory* td;
   TDirectory* td1;
   TDirectory* td2;
   TDirectory* td3;
   TDirectory* td4;

   Double_t binning[39] ={1., 3., 5., 7., 9.,  10., 12., 15., 18., 20., 22., 24., 26., 28.,
                          29., 30., 31., 32., 33., 34., 35., 36., 37., 38., 39., 40., 41.,
                          42., 43., 45., 50., 60., 70., 100., 200., 300., 400., 600., 1000.};

};

#endif












