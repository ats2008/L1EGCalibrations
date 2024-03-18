//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Tue Oct 19 09:59:39 2021 by ROOT version 6.24/06
// from TTree L1UpgradeTree/L1UpgradeTree
// found on file: L1Ntuple.root
//////////////////////////////////////////////////////////

#ifndef IsolationAnalysis_h
#define IsolationAnalysis_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <TProfile.h>
#include <TString.h>
#include <chrono>
// Header file for the classes stored in the TTree if any.

class IsolationAnalysis {
public :
  TChain          *fChain;   //!pointer to the analyzed TTree or TChain
  Int_t           fCurrent; //!current Tree number in a TChain
  
  // Fixed size dimensions of array or collections stored in the TTree if any.
  
  // Declaration of leaf types

  // Int_t          et;
  // Int_t          eta;
  // Int_t          ntt;
  // Int_t          iso;
  float et;
  float eta;
  float ntt;
  float iso;
  float eleProbeEta	;
  float eleProbePhi	;
  float eleTagEta   ;
  float eleTagPhi	;
  float Nvtx	;
  //Int_t isProbeLoose;
  float isProbeLoose;

  Long64_t       maxEntries;
  Int_t          reportEvery;
  IsolationAnalysis(const std::string& inputFileName);
  ~IsolationAnalysis();

  void accessTree(std::string & input_filelist);
  void readParameters(const std::string jfile);
  
  void readTree();
  void analyse();
  
  void bookHistograms(const std::string option);
  void saveHistograms();
  
  
  void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters);
  
  TFile* outputFile_;
  
  Long64_t nEntries_;
  std::string ntupleFileName_;
  std::string outputFileName_;
  std::string outputWPFileName_;
  int etMax_;
  bool bookedHistograms_;
  
  void readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map);
  bool getHistoBin(TString str, short& et_bin, short& eta_bin, short& ntt_bin);
  TString getHistoName(short et, short eta, short ntt);
  Double_t findEfficiencyProgression(Double_t IEt, Double_t MinPt,
				     Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at);
  Double_t findEfficiencyProgressionForV3(Double_t IEt, Double_t MinPt,
				     Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at);
  void fillLUTProgression(const std::string option);
  short getUpdatedBin(short value, std::string type);
    
  UInt_t nBinsIEta;
  UInt_t nBinsIEt; 
  UInt_t nBinsnTT;

  UInt_t tmpFitMin;
  UInt_t tmpFitMax;
  
  TProfile* hprof_IEt;
  TProfile* hprof_IEta;
  TProfile* hprof_nTT;
  
  TH1F* pt_all;
  TH1F* eta_all;
  TH1F* nTT_all;
  
  TH1F* pt_large_eta;
  TH1F* eta_large_eta;
  TH1F* nTT_large_eta;
  TH1F* iso_large_eta;
  
  //Efficiency as function of pT
  std::map<Int_t,TH1F*> pt_pass_efficiency ;
  std::map<Int_t,TH1F*> pt_pass_efficiency_TH3 ;
  std::map<Int_t,TH1F*> eta_pass_efficiency ;
  std::map<Int_t,TH1F*> nTT_pass_efficiency ;
  
  std::map<short, short> lutMapIEt_;
  std::map<short, short> lutMapIEta_;
  std::map<short, short> lutMapnTT_;
   
  std::map<TString,TH1F*> Histos_PerBin ;
  std::map<Int_t,TH3F*> IsoCut_PerBin ;
  std::map<Int_t,std::map<TString,Int_t>> IsoCut_PerEfficiency_PerBin;
  
  std::vector<std::string> lutProgOptVec_;
  std::map<std::string, TH3F*> lutProgHistoMap_;
  std::map<std::string, TH3F*> lutProgHistoMap_v2_;
  std::map<std::string, TH3F*> lutProgHistoMap_v3_;
  std::vector<TH3F*> LUT_WP ;
  
  std::vector<short> lutIEtaVec_;
  std::vector<short> lutIEtVec_;
  std::vector<short> lutnTTVec_;
  std::vector<short>  updatedIEtaVec_;
  std::vector<short>  updatedIEtVec_;
  std::vector<short>  updatednTTVec_;  
    
  int doSuperCompression;
  bool doDynamicBinning;
  std::vector<Int_t> superCompressionToDefaultCompressionMap;
  std::vector<Float_t> puRewightMap;

  Int_t isoOffset;
  Int_t ietaMinForOffset;
  Int_t ietaMaxForOffset;

};

#endif
