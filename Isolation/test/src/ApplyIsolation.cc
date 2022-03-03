#include "ApplyIsolation.h"
#include<stdlib.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <TSelector.h>
#include <iostream>
#include <TLorentzVector.h>
#include <TH1.h>
#include <TH2.h>
#include <TH3.h>
#include <TMath.h>
#include <TCanvas.h>
#include <TLegend.h>
#include <TPaveText.h>
#include <TStyle.h>
#include <TROOT.h>
#include <sstream>
#include <TBranchElement.h>
#include <fstream>
#include <map>
#include <sstream>
#include <iostream>
#include <string>
#include <TGraphAsymmErrors.h>
#include <chrono>

ApplyIsolation::ApplyIsolation(std::string& inputFileName){
  readParameters(inputFileName);
  if (ntupleFileNameRate_.size() == 0) {
    std::cout << " Inputfile for rate list missing !!!" << ntupleFileNameRate_ << std::endl;
    return;
  }

  if (ntupleFileNameTurnOn_.size() == 0) {
    std::cout << " Inputfile for turn on list missing !!!" << ntupleFileNameTurnOn_ << std::endl;
    return;
  }

  accessTree(ntupleFileNameRate_ , ntupleFileNameTurnOn_);
  
  assert(fChain);

  bookedHistograms_ = false;

}

ApplyIsolation::~ApplyIsolation() {
}

void ApplyIsolation::accessTree(std::string & input_filelist_rate , std::string & input_filelist_turn_on) {
  std::ifstream myFileRate;
  myFileRate.open(input_filelist_rate.c_str(), std::ios::in);
  if (!myFileRate) {
    std::cout << "Input File: " << input_filelist_rate << " could not be opened!" << std::endl;
    fChain = 0;
  } else {
    fChain = new TChain("L1UpgradeTree");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileRate.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                     
      std::string line(buf);
      fChain->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
    }
    //    fChain->Print();                                                                                                             
  }

  std::ifstream myFileTurnOn;
  myFileTurnOn.open(input_filelist_turn_on.c_str(), std::ios::in);

  if (!myFileTurnOn) {
    std::cout << "Input File: " << input_filelist_turn_on << " could not be opened!" << std::endl;
    fChain1 = 0;
  } 
  else {
    fChain1 = new TChain("Ntuplizer/TagAndProbe");
    static constexpr int BUF_SIZE = 256;
    char buf[BUF_SIZE];
    while (myFileTurnOn.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character                                
      std::string line(buf);
      fChain1->AddFile(line.c_str(),-1);
      std::cout << "Adding file " << line << " Entries " << fChain1->GetEntries() <<  std::endl;
    }
    //    fChain1->Print();
  }
 }

void ApplyIsolation::loops() {
  short in_compressediEta;
  short in_compressediEt;
  short in_compressedNTT;
  Long64_t nbytes = 0, nb = 0;

  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
  optionsFile_->cd("Step2Histos");
  
  //Event loop for rate                                                                                                                                                                                  
  if (fChain == 0) return;
  nEntries_ = fChain->GetEntriesFast();
  std::cout << "Total available  Entries for Rate: " << nEntries_ << std::endl;
  if(maxEntriesForRate_ >0 ) nEntries_= maxEntriesForRate_ < nEntries_ ? maxEntriesForRate_ : nEntries_;
  std::cout<<"Processing a total of: "<<nEntries_<<" Entries \n";

  auto t_start = std::chrono::high_resolution_clock::now();
  auto t_end = std::chrono::high_resolution_clock::now();

  int den=0;
  for (Long64_t jentry=0; jentry < nEntries_; jentry++) {
    Long64_t ientry = fChain->LoadTree(jentry);

    if (ientry < 0) break;
    nb = fChain->GetEntry(jentry);   nbytes += nb;

    if(jentry%reportEvery_ == 0 )
      {
	t_end = std::chrono::high_resolution_clock::now();
	std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nEntries_<<"  [ "<<100.0*jentry/nEntries_<<"  % ]  "
		 << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
		 <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries_ - jentry)/(1e-9 + jentry)* 0.001
		 <<std::endl;
      }


    
    if (nEGs < 1) continue;
    for (auto it :lutProgOptVec_) {
      bool Filled_Progression= kFALSE;
      Int_t IsoCut_Progression;
      TString ResultProgressionName_= "LUT_Progression";
      TString pt_Progression_= "pt_Progression";
      ResultProgressionName_ += it;
      pt_Progression_ +=it;
      
      TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
      
      for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
	
	if (egBx[iEG]!=0)   continue;
	
	float EG_Et  = egEt[iEG];
	short EG_NTT = egNTT[iEG];
	short EG_TowerIEta = egTowerIEta[iEG];
	short EG_Iso_Et = egIsoEt[iEG];
	short EG_Raw_Et = egRawEt[iEG];

	std::map<short, short>::iterator EtaPos = lutMapEta.find(int(abs(EG_TowerIEta)));
	if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
	else in_compressediEta = nBinsIEta-1;
	
	std::map<short, short>::iterator EtPos = lutMapEt.find(int(EG_Raw_Et));
	if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
	else {
	  in_compressediEt = nBinsIEt-1;
	  EG_Et = EG_Et > ET_MAX ? ET_MAX : EG_Et;
	}
	
	std::map<short, short>::iterator NTTPos = lutMapNTT.find(int(EG_NTT));
	if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
	else in_compressedNTT = nBinsNTT-1;
        
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	
	if(!Filled_Progression && EG_Iso_Et<=IsoCut_Progression) {
	  ptMap_[pt_Progression_]->Fill(EG_Et);
	  Filled_Progression = kTRUE;
	}
      }
    }
    den++;
  }
  std::cout<<den<<","<<nEntries_<<std::endl;
  double scale = (11.2456 * bunches)/nEntries_;

  //Filling Rate Histos  
  for(UInt_t i=0;i< ET_MAX;i++) {
    for (auto it :lutProgOptVec_) {
      TString CurrentNameHisto = "pt_Progression";
      TString CurrentNameHisto1= "rate_Progression";
      CurrentNameHisto += it;
      CurrentNameHisto1 += it;
      rateMap_[CurrentNameHisto1]->SetBinContent(i+1, ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale);
    }
  }
  
  //Finding Et for fixed rate (frate)
  int xbin;
  for (auto it :lutProgOptVec_) {
    TString CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto1 += it;
    float vl,bl; //last value, last bin
    for(UInt_t i=0;i< ET_MAX;i++) {
      float rate = rateMap_[CurrentNameHisto1]->GetBinContent(i+1);
      if(rate > frate_) {
        bl =i;
	vl=rate;
      }
      
      if(rate <=frate_){
	xbin = (rate - frate_) < abs(rate-vl) ? i : bl;
	break;    
      }
    }     

    std:: cout<<"option: "<<it<<" , "<<"et: "<<xbin<<" +/- " <<vl/sqrt(vl/scale)<<std::endl;

    et_option.push_back(xbin);
  }


  //booking histograms for fixed Et for all options
  int n=-1;
  for (auto it :lutProgOptVec_) {    
    n++;
    TString PtPassName_fr_ = "pt_pass_option";
    TString turnOn_Option_fr_ = "turnOn_Option";
    
    PtPassName_fr_ += it;
    PtPassName_fr_ +="_Et_";
    PtPassName_fr_=PtPassName_fr_ + std::to_string(et_option[n])  + "_fr";
    td2->cd();
    TH1F*PtPass_fr= new TH1F(PtPassName_fr_, PtPassName_fr_, nBins_fine,xEdges_fine);
    pt_pass_Map_.insert(std::make_pair(PtPassName_fr_, PtPass_fr));


    turnOn_Option_fr_ += it;
    turnOn_Option_fr_ +="_Et_";
    turnOn_Option_fr_ = turnOn_Option_fr_ + std::to_string(et_option[n]) + "_fr";
    
    if(!check_turn_on_dir) {                                                                                                     
      td3 = outputFile_->mkdir("turn_on_progression");                                                                           
      check_turn_on_dir = true;                                                                                                  
    }                                                                                                                            
    td3->cd();
    TGraphAsymmErrors* turnOn_Option_fr;
    turnOn_Map_.insert(std::make_pair(turnOn_Option_fr_,turnOn_Option_fr));
  }

  
  optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
  optionsFile_->cd("Step2Histos");     

//Event loop for turnon
  if (fChain1 == 0) return;

  nEntries1_ = fChain1->GetEntries();
  Long64_t nbytes1 = 0, nb1 = 0;
  double sum=0;
  std::cout << "Total available  Entries For Efficiency " << nEntries1_ << std::endl;
  if(maxEntriesForEfficiency >0 ) nEntries1_= maxEntriesForEfficiency < nEntries1_ ? maxEntriesForEfficiency : nEntries1_;


  t_start = std::chrono::high_resolution_clock::now();
  for (Long64_t jentry=0; jentry< nEntries1_; jentry++) {
    Long64_t ientry = fChain1->LoadTree(jentry);
    if (ientry < 0) break;
    nb1 = fChain1->GetEntry(jentry);   nbytes1 += nb1;

 //   if(l1tEmuNTT<60 && nTTRange_) continue;  //define range
    if(l1tEmuRawEt < 0.) continue;
    
    if(jentry%reportEvery_ == 0 )
      {
	t_end = std::chrono::high_resolution_clock::now();
	std::cout<<"Processing Entry in event loop (TurnOn) : "<<jentry<<" / "<<nEntries1_<<"  [ "<<100.0*jentry/nEntries1_<<"  % ]  "
		 << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
		 <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries1_ - jentry)/(1e-9 + jentry)* 0.001
		 <<std::endl;
      }


    // if(l1tEmuNTT<60 && nTTRange_) continue;  //define range
    //    sum++;
    if(!( isProbeLoose==1 && fabs(eleProbeEta) < 2.5  && sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2))>0.6)) continue;
    // sum1++;
    //    pT_all->Fill(eleProbeSclEt);
    
    if(l1tEmuRawEt < 0.) continue;
    pT_all->Fill(eleProbeSclEt);
    sum++;

    
    std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(l1tEmuTowerIEta));
    if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
    else in_compressediEta = nBinsIEta-1;
    
    std::map<short, short>::iterator EtPos = lutMapEt.find(l1tEmuRawEt);
    if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
    else in_compressediEt = nBinsIEt-1;

    std::map<short, short>::iterator NTTPos = lutMapNTT.find(l1tEmuNTT);
    if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
    else in_compressedNTT = nBinsNTT-1;
  

    Int_t IsoCut_Progression;
    for(UInt_t e = 15 ; e <= 50 ; e += 5) {
      for (auto it :lutProgOptVec_) {
	TString ResultProgressionName_= "LUT_Progression";
        TString PtPassName_= "pt_pass_option";
	
	ResultProgressionName_ += it;
	PtPassName_ += it;
	PtPassName_ +="_Et_";
	PtPassName_ += std::to_string(e);
	
	TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
	IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
	
	if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression){

	  pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
	}
      }
    }
    int n=-1;
    for (auto it :lutProgOptVec_) {
      n++;
      TString PtPassName_fr= "pt_pass_option";
      TString ResultProgressionName_= "LUT_Progression";      
      
      PtPassName_fr += it;
      PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";
      ResultProgressionName_ += it;
      
      TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
      IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
      if(l1tEmuPt >= et_option[n] && l1tEmuIsoEt <= IsoCut_Progression){
	pt_pass_Map_[PtPassName_fr]->Fill(eleProbeSclEt);
      }
    }


    // using Run2 LUT                                                                                                              
    if(hasL1Emu_26==1) pt_pass_Map_Run2_["pt_pass_Et_26"]->Fill(eleProbeSclEt);
    if(hasL1Emu_24==1) pt_pass_Map_Run2_["pt_pass_Et_24"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_looseiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_loose"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso26==1) pt_pass_Map_Run2_["pt_pass_Et_26_tight"]->Fill(eleProbeSclEt);
    if(hasL1Emu_tightiso24==1) pt_pass_Map_Run2_["pt_pass_Et_24_tight"]->Fill(eleProbeSclEt);

  }

  t_end = std::chrono::high_resolution_clock::now();
  std::cout << " Turn Ons Done: Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0<<std::endl;

  //Turnons for various Ets
  for(UInt_t e = 15 ; e <= 50 ; e += 5) {
    for (auto it :lutProgOptVec_) {
      TString PtPassName_= "pt_pass_option";
      TString turnOn_Option_="turnOn_Option";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);
      
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);
      outputFile_->cd();
      if(!check_turn_on_dir) {
	td3 = outputFile_->mkdir("turn_on_progression");
        check_turn_on_dir = true;
      }
      td3->cd();
      turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp"); 
      turnOn_Map_[turnOn_Option_]->Write();

      double acceptance = pt_pass_Map_[PtPassName_]->GetEntries() ;
      th1fStore["FixedEtTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),e);
      th1fStore["FixedEtTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),acceptance);

    }
  }

  //Turnons for fixed Rate 
  n=-1;
  for (auto it :lutProgOptVec_) {
    n++;
    TString PtPassName_fr= "pt_pass_option";
    TString turnOn_Option_fr="turnOn_Option";
    PtPassName_fr += it;
    PtPassName_fr = PtPassName_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";
    
    td3->cd();
    turnOn_Option_fr += it;
    turnOn_Option_fr = turnOn_Option_fr + "_Et_" + std::to_string(et_option[n]) +"_fr";;
    
    turnOn_Map_[turnOn_Option_fr]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_fr],pT_all,"cp");
    turnOn_Map_[turnOn_Option_fr]->Write();

    double acceptance = pt_pass_Map_[PtPassName_fr]->GetEntries() ;
    th1fStore["FixedRateTurnons"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),et_option[n]);
    th1fStore["FixedRateTurnons_Acceptance"]->Fill(turnOn_Map_[turnOn_Option_fr]->GetName(),acceptance);

  }

  // using Run2 LUT                                                                                                                                      
  td4->cd();
  std::vector<int> v = {24, 26};
  for(int e : v) {
    TString PtPassName_= "pt_pass";
    TString turnOnName_="turnOn";
    PtPassName_ +="_Et_";
    PtPassName_ += std::to_string(e);
    TString PtPassLooseIsoName_=PtPassName_+"_loose";
    TString PtPassTightIsoName_=PtPassName_+"_tight";

    turnOnName_ +="_Et_";
    turnOnName_ += std::to_string(e);

    TString turnOnLooseIsoName_=turnOnName_+"_loose";
    TString turnOnTightIsoName_=turnOnName_+"_tight";
    turnOn_Map_Run2_[turnOnName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassName_],pT_all,"cp");
    turnOn_Map_Run2_[turnOnLooseIsoName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassLooseIsoName_],pT_all,"cp");
    turnOn_Map_Run2_[turnOnTightIsoName_]=new TGraphAsymmErrors( pt_pass_Map_Run2_[PtPassTightIsoName_],pT_all,"cp");

    turnOn_Map_Run2_[turnOnName_]->Write();
    turnOn_Map_Run2_[turnOnLooseIsoName_]->Write();
    turnOn_Map_Run2_[turnOnTightIsoName_]->Write();

    double acceptance = pt_pass_Map_Run2_[PtPassName_]->GetEntries() ;
    double acceptance_loose = pt_pass_Map_Run2_[PtPassLooseIsoName_]->GetEntries() ;
    double acceptance_tight = pt_pass_Map_Run2_[PtPassTightIsoName_]->GetEntries() ;

    th1fStore["Run2Turnons_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnName_]->GetName(),acceptance);
    th1fStore["Run2Turnons_LooseIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnLooseIsoName_]->GetName(),acceptance_loose);
    th1fStore["Run2Turnons_TightIso_Acceptance"]->Fill(turnOn_Map_Run2_[turnOnTightIsoName_]->GetName(),acceptance_tight);
  }

}

void ApplyIsolation::bookHistograms() {
  outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
  outputFile_->cd();

  for (auto it : lutProgOptVec_) {    
    TString  CurrentNameHisto = "pt_Progression";
    TString CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto += it;
    CurrentNameHisto1 += it;
    
    if(!check_pt_rate_dir) { 
      td = outputFile_->mkdir("pt_progression_for_rate");
      check_pt_rate_dir = true;
    }
    td->cd();
    TH1F* pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto ,ET_MAX ,  0.0 , ET_MAX);
    ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));
    
    if(!check_rate_dir) {
      td1 = outputFile_->mkdir("rate_progression");
      check_rate_dir = true;
    }
    td1->cd();   
    TH1F* rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX ,  0.0 , ET_MAX );
    rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
    
    for(UInt_t e = 15 ; e <= 50 ; e += 5) {
      TString PtPassName_= "pt_pass_option";
      TString turnOn_Option_="turnOn_Option";
      
      PtPassName_ += it;
      PtPassName_ +="_Et_";
      PtPassName_ += std::to_string(e);  
      
      if(!check_pt_turn_on_dir) {
	td2 = outputFile_->mkdir("pt_progression_for_turnon");
	check_pt_turn_on_dir = true;
      }
      td2->cd();
      TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, 38,binning);
      pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));
      
      turnOn_Option_ += it;
      turnOn_Option_ +="_Et_";
      turnOn_Option_ += std::to_string(e);
      
      TGraphAsymmErrors* turnOn_Option;
      turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
    }
  }
  td2->cd();
  pT_all = new TH1F("pT_all","pT_all",38,binning);

  th1fStore["FixedRateTurnons"] = new TH1F("FixedRateTurnons","",3,0.0,3.0);
  th1fStore["FixedRateTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedRateTurnons_Acceptance"] = new TH1F("FixedRateTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedRateTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);  

  th1fStore["FixedEtTurnons"]   = new TH1F("FixedEtTurnons","",3,0.0,3.0);
  th1fStore["FixedEtTurnons"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["FixedEtTurnons_Acceptance"]   = new TH1F("FixedEtTurnons_Acceptance","",3,0.0,3.0);
  th1fStore["FixedEtTurnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);


  //using Run2 LUT
  std::vector<int> v = {24, 26};
  for(int e : v) {
    TString PtPassName_= "pt_pass";
    TString turnOnName_="turnOn";
    PtPassName_ +="_Et_";
    PtPassName_ += std::to_string(e);
    TString PtPassLooseIsoName_=PtPassName_+"_loose";
    TString PtPassTightIsoName_=PtPassName_+"_tight";
    
    turnOnName_ +="_Et_";
    turnOnName_ += std::to_string(e);
    
    TString turnOnLooseIsoName_=turnOnName_+"_loose";
    TString turnOnTightIsoName_=turnOnName_+"_tight";
    
    if(!check_turn_on_dir_Run2) {
      td4 = outputFile_->mkdir("turnon_progression_Run2");
      check_turn_on_dir_Run2 = true;
    }
    td4->cd();
    //TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, 38,binning);
    //TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, 38,binning);
    //TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, 38,binning);
    TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
    TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, nBins_fine,xEdges_fine);
    TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, nBins_fine,xEdges_fine);
    
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassName_, PtPass));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassLooseIsoName_, PtPassLoose));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassTightIsoName_, PtPassTight));
    
    
    TGraphAsymmErrors* turnOn;
    TGraphAsymmErrors* turnOnLoose;
    TGraphAsymmErrors* turnOnTight;

    turnOn_Map_Run2_.insert(std::make_pair(turnOnName_,turnOn));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnLooseIsoName_,turnOnLoose));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnTightIsoName_,turnOnTight));
  }


  th1fStore["Run2Turnons_Acceptance"]   = new TH1F("Run2Turnons_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]   = new TH1F("Run2Turnons_LooseIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_TightIso_Acceptance"]   = new TH1F("Run2Turnons_TightIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_TightIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  

  //using Run2 LUT
  std::vector<int> v = {24, 26};
  for(int e : v) {
    TString PtPassName_= "pt_pass";
    TString turnOnName_="turnOn";
    PtPassName_ +="_Et_";
    PtPassName_ += std::to_string(e);
    TString PtPassLooseIsoName_=PtPassName_+"_loose";
    TString PtPassTightIsoName_=PtPassName_+"_tight";
    
    turnOnName_ +="_Et_";
    turnOnName_ += std::to_string(e);
    
    TString turnOnLooseIsoName_=turnOnName_+"_loose";
    TString turnOnTightIsoName_=turnOnName_+"_tight";
    
    if(!check_turn_on_dir_Run2) {
      td4 = outputFile_->mkdir("turnon_progression_Run2");
      check_turn_on_dir_Run2 = true;
    }
    td4->cd();
    TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, 38,binning);
    TH1F*PtPassLoose= new TH1F(PtPassLooseIsoName_, PtPassLooseIsoName_, 38,binning);
    TH1F*PtPassTight= new TH1F(PtPassTightIsoName_, PtPassTightIsoName_, 38,binning);
    
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassName_, PtPass));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassLooseIsoName_, PtPassLoose));
    pt_pass_Map_Run2_.insert(std::make_pair(PtPassTightIsoName_, PtPassTight));
    
    
    TGraphAsymmErrors* turnOn;
    TGraphAsymmErrors* turnOnLoose;
    TGraphAsymmErrors* turnOnTight;

    turnOn_Map_Run2_.insert(std::make_pair(turnOnName_,turnOn));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnLooseIsoName_,turnOnLoose));
    turnOn_Map_Run2_.insert(std::make_pair(turnOnTightIsoName_,turnOnTight));
  }


  th1fStore["Run2Turnons_Acceptance"]   = new TH1F("Run2Turnons_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]   = new TH1F("Run2Turnons_LooseIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_LooseIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  th1fStore["Run2Turnons_TightIso_Acceptance"]   = new TH1F("Run2Turnons_TightIso_Acceptance","",3,0.0,3.0);
  th1fStore["Run2Turnons_TightIso_Acceptance"]->SetCanExtend(TH1::kAllAxes);
  
  bookedHistograms_ = true;
}

void ApplyIsolation::readLUTTable(std::string& file_name, unsigned int& nbin, std::map<short, short>& lut_map) {
  std::cout << "Opening LUT file " << file_name << std::endl;
  std::ifstream lutFile(file_name.c_str());
  
  if (!lutFile) {
    std::cerr << "Input File: " << file_name << " could not be opened!" << std::endl;
    return;
  }
  std::string line;
  nbin = 0;
  if(lutFile.is_open()) {
    nbin = 0;
    while(std::getline(lutFile,line)) {
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      if (line.length() == 0) continue;
      std::vector<std::string> tokens;
      tokenize(line,tokens," ");
      unsigned int key   = std::atoi(tokens.at(0).c_str());
      unsigned int value = std::atoi(tokens.at(1).c_str());
      lut_map.insert({ key, value });
      if (nbin < value) nbin = value;
    }
    nbin += 1;
  }
  std::cout << " nbin " << nbin << std::endl;
}

void ApplyIsolation::saveHistograms() {
  if (outputFile_ && bookedHistograms_) {
    outputFile_->cd();
    for (std::map<TString,TH1F *>::iterator it=th1fStore.begin() ; it!=th1fStore.end(); ++it)
      {
        auto &ahist = *(it->second); 
        ahist.Write();
      }
    outputFile_->Write();
    outputFile_->Close();
  }
}

void ApplyIsolation::readParameters(const std::string jfile) {
  std::cout << jfile << std::endl;
  std::ifstream jobcardFile(jfile.c_str());
   if (!jobcardFile) {
    std::cerr << "Input File: " << jfile << " could not be opened!" << std::endl;
    return;
  }
  std::string line;
  if(jobcardFile.is_open()) {
    while(std::getline(jobcardFile,line)) {
      // enable '#' and '//' style comments                                                                                                                                                   
      if (line.substr(0,1) == "#" || line.substr(0,2) == "//") continue;
      std::vector<std::string> tokens;
      tokenize(line,tokens,"=");
      std::cout << tokens[0] << ":" << tokens[1] << std::endl;
      std::string key   = tokens.at(0);
      std::string value = tokens.at(1);

      if(key=="NtupleFileNameRate")        ntupleFileNameRate_= value;
      else if(key=="NtupleFileNameTurnOn")        ntupleFileNameTurnOn_= value;
      else if (key=="OptionsFileName") optionsFileName_ = value.c_str();
      else if (key=="OutputFileName")  outputFileName_ = value.c_str();
      else if (key=="EtLUTFileName")    readLUTTable(value,nBinsIEt, lutMapEt);
      else if (key=="EtaLUTFileName")   readLUTTable(value,nBinsIEta, lutMapEta);
      else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
      else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());

      else if (key=="MaxEntriesForEfficency") maxEntriesForEfficiency = atoi(value.c_str());
      else if (key=="MaxEntriesForRate")   maxEntriesForRate = atoi(value.c_str());
      else if (key=="FixedRate") frate = std::stof(value.c_str());
      else if (key=="ReportEvery")    reportEvery = atoi(value.c_str());
      else if (key=="LUTProgressionOptions")
        {
	  std::string tmp_string = value;
	  std::vector<std::string> tmp_vec;
          tokenize(tmp_string,tmp_vec,",");
          for (auto it : tmp_vec) {
            lutProgOptVec_.push_back(it);
          }
        }
      else 
	std::cout << " unknown option " << " key " << key << std::endl;
    }
  }
  jobcardFile.close();
}


void ApplyIsolation::readTree() {
  fChain->SetMakeClass(1);
  //for rate
  fChain->SetBranchAddress("nEGs", &nEGs);
  fChain->SetBranchAddress("egEt", &egEt);
  fChain->SetBranchAddress("egTowerIEta", &egTowerIEta); 
  fChain->SetBranchAddress("egNTT", &egNTT);
  fChain->SetBranchAddress("egBx", &egBx);
  fChain->SetBranchAddress("egIsoEt", &egIsoEt);
  fChain->SetBranchAddress("egRawEt", &egRawEt);

  //for turnon
  fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
  fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
  fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
  fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
  fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
  fChain1->SetBranchAddress("l1tEmuIsoEt",&l1tEmuIsoEt);


  fChain1->SetBranchAddress("eleProbeEta"  ,&eleProbeEta			);
  fChain1->SetBranchAddress("eleProbePhi"  ,&eleProbePhi			);
  fChain1->SetBranchAddress("eleTagEta"    ,&eleTagEta		    	);
  fChain1->SetBranchAddress("eleTagPhi"    ,&eleTagPhi			    );
  fChain1->SetBranchAddress("isProbeLoose" ,&isProbeLoose			);
  fChain1->SetBranchAddress("hasL1Emu_26",&hasL1Emu_26);
  fChain1->SetBranchAddress("hasL1Emu_24",&hasL1Emu_24);
  fChain1->SetBranchAddress("hasL1Emu_looseiso26",&hasL1Emu_looseiso26);
  fChain1->SetBranchAddress("hasL1Emu_looseiso24",&hasL1Emu_looseiso24);
  fChain1->SetBranchAddress("hasL1Emu_tightiso26",&hasL1Emu_tightiso26);

}


void ApplyIsolation::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
  // Skip delimiters at beginning.                                                                                          
  std::string::size_type lastPos = str.find_first_not_of(delimiters, 0);
  
  // Find first "non-delimiter".                                                                                            
  std::string::size_type pos = str.find_first_of(delimiters, lastPos);
  
  while (std::string::npos != pos || std::string::npos != lastPos)  {
    // Found a token, add it to the vector.                                                                      
    tokens.push_back(str.substr(lastPos, pos - lastPos));
    
    // Skip delimiters.  Note the "not_of"                                                                       
    lastPos = str.find_first_not_of(delimiters, pos);
    
    // Find next "non-delimiter"                                                                                    
    pos = str.find_first_of(delimiters, lastPos);
  }
}

int main(int argc,char *argv[]) {
  if (argc == 1) {
    std::cout << " No option provided!!!" << std::endl;
    return 1;
  }
  std::string data_file = argv[1];
  ApplyIsolation treeReader(data_file);
  treeReader.readTree();
  treeReader.bookHistograms();
  std::cout << " Calling Loop" << std::endl;
  treeReader.loops();
  treeReader.saveHistograms();
  return 0;
}


