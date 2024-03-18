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

ApplyIsolation::ApplyIsolation(std::string& inputFileName) {

    for(Int_t i=0 ; i<= nBins_fine; i++)
    {
        xEdges_fine[i]=0.0+i*1.0;
    }
    readParameters(inputFileName);
    if (ntupleFileNameRate_.size() == 0) {
        std::cout << " Inputfile for rate list missing !!!" << ntupleFileNameRate_ << std::endl;
        return;
    }

    if (ntupleFileNameTurnOn_.size() == 0) {
        std::cout << " Inputfile for turn on list missing !!!" << ntupleFileNameTurnOn_ << std::endl;
        return;
    }

    accessTree(ntupleFileNameRate_, ntupleFileNameTurnOn_);

    assert(fChain);

    bookedHistograms_ = false;

}

ApplyIsolation::~ApplyIsolation() {
}

void ApplyIsolation::accessTree(std::string & input_filelist_rate, std::string & input_filelist_turn_on) {
    std::ifstream myFileRate;
    myFileRate.open(input_filelist_rate.c_str(), std::ios::in);
    if (!myFileRate) {
        std::cout << "Input File: " << input_filelist_rate << " could not be opened!" << std::endl;
        fChain = 0;
    }
    else {

        fChain = new TChain("l1UpgradeEmuTree/L1UpgradeTree");
        fChain_1= new TChain("l1EventTree/L1EventTree");
        static constexpr int BUF_SIZE = 2569;
        char buf[BUF_SIZE];
        while (myFileRate.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character

            std::string line(buf);

            fChain->AddFile(line.c_str(),-1);
            fChain_1->AddFile(line.c_str(),-1);
            std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
            //fChain->Print();
            //fChain_1->Print();

        }
    }

    std::ifstream myFileTurnOn;
    myFileTurnOn.open(input_filelist_turn_on.c_str(), std::ios::in);

    if (!myFileTurnOn) {
        std::cout << "Input File: " << input_filelist_turn_on << " could not be opened!" << std::endl;
        fChain1 = 0;
    }
    else {
        //    fChain1 = new TChain("test");
        fChain1 =  new TChain("Ntuplizer/TagAndProbe");
        static constexpr int BUF_SIZE = 2569;
        char buf[BUF_SIZE];
        while (myFileTurnOn.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
            std::string line(buf);
            fChain1->AddFile(line.c_str(),-1);
            std::cout << "Adding file " << line << " Entries " << fChain1->GetEntries() <<  std::endl;
        }
        //    fChain1->Print();
    }
}

void ApplyIsolation::createEfficiencyHistograms(TH1* h_P, TH1* h_F, TH1* h_E)
{
    Double_t num, den, ratio, eL, eU;
    double a  = 0.3173;
    double aeff = (1-a)/2;

    for (int k=1; k<=h_P->GetNbinsX(); ++k) {
        num = h_P->GetBinContent(k);
        den = h_P->GetBinContent(k) + h_F->GetBinContent(k);
        eU = 0.0;
        eL = 0.0;
        h_E->SetBinErrorOption(TH1::kPoisson);
        if(den!=0 && num != 0) {
            ratio = num/den;
            if (num > 100 || den*(1-ratio) > 100) {
                eU = sqrt(ratio*(1-ratio)/den);
                eL = sqrt(ratio*(1-ratio)/den);
            } else {
                eU = (1-BetaInverse(aeff,den-num,num+1))-ratio;
                eL = ratio-(1-BetaInverse(1-aeff,den-num+1,num));
            }

        } else {
            ratio = 0.0;
            eL = 0.0;
            eU = 0;
        }
        h_E->SetBinContent(k,ratio);
        h_E->SetBinError(k,eL);
        h_E->SetStats(0);
    }
}


void ApplyIsolation::loops() {
    short in_compressediEta;
    short in_compressediEt;
    short in_compressedNTT;
    Long64_t nbytes = 0, nbytes_1=0, nb = 0, nb_1=0;
    //  weightFile_ = new TFile(weightFileName_.c_str(), "READ");
    //TH1D*w= (TH1D*)weightFile_->Get("ratio_hist");
    optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
    optionsFile_->cd("Step2Histos");

    //Event loop for rate  ( SIngle EG + Double EG)
    if (fChain == 0) return;
    nEntries_ = fChain->GetEntriesFast();
    std::cout << "Total available  Entries for Rate: " << nEntries_ << std::endl;
    if(maxEntriesForRate_ >0 ) nEntries_= maxEntriesForRate_ < nEntries_ ? maxEntriesForRate_ : nEntries_;
    std::cout<<"Processing a total of: "<<nEntries_<<" Entries \n";

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();
    double sum_weight=0.;
    for (Long64_t jentry=0; jentry < nEntries_; jentry++) {

        Long64_t ientry = fChain->LoadTree(jentry);
        Long64_t i1entry = fChain_1->LoadTree(jentry);

        if (ientry < 0) break;
        if (i1entry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        nb_1 = fChain_1->GetEntry(jentry);
        nbytes_1 += nb_1;

        if(jentry%reportEvery_ == 0 )
        {
            t_end = std::chrono::high_resolution_clock::now();
            std::cout<<"Processing Entry in event loop (Rate) : "<<jentry<<" / "<<nEntries_<<"  [ "<<100.0*jentry/nEntries_<<"  % ]  "
                     << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                     <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries_ - jentry)/(1e-9 + jentry)* 0.001
                     <<std::endl;
        }
        //    weight=w->GetBinContent(nPV_True+1);
        //    std::cout<<"weight ="<<weight<<","<<nPV_True<<std::endl;
        //    sum_weight=sum_weight+weight;

        ///    if (!(nPV_True <=55 && nPV_True >=49)) continue;
        if (!(run = 362616)) continue;
        sum_weight++;

        if (nEGs < 1) continue;

        for (auto it :lutProgOptVec_) {
            bool Filled_Progression= kFALSE;
            bool Filled_ProgressionD= kFALSE;
            Int_t IsoCut_Progression;
            TString ResultProgressionName_= "LUT_Progression_" +it;
            TString pt_Progression_= "pt_Progression" + it;
            TString pt_ProgressionD_ = "pt_Progression_double" +it;
            TString pt_ProgressionDI_= "pt_Progression_doubleIso" +it;
            TString pt_ProgressionDIER_= "pt_Progression_doubleIsoER" +it;
            TString pt_ProgressionDIXXp10ER_= "pt_Progression_doubleIsoXXp10ER" +it;
            TString pt_ProgressionDIXXp5ER_ = "pt_Progression_doubleIsoXXp5ER" +it;
             TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
            if(!ResultProgressionName) {
    
                        std::cout<<"LOLLLLLLL | "<<ResultProgressionName_<<std::endl;
                    }

            Int_t nIsoEGs=0;
            Int_t nIsoEREGs=0;
            Int_t nIsoERXXp10EGs=0;
            Int_t nIsoERXXp5EGs=0;

            float doubleIso_threshold(1e9);
            float doubleIso_thresholdER(1e9);
            float doubleIsoXXp10_thresholdER(1e9);
            float doubleIsoXXp5_thresholdER(1e9);

            for (UShort_t iEG=0; iEG < nEGs; ++iEG) {
                if (egBx[iEG]!=0)   continue;
                float EG_Et  = egEt[iEG];
                float EG_Eta  = egEta[iEG];
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

                
                if(nEGs >=2 && !Filled_ProgressionD) {
                    if(EG_Et > 10. && EG_Iso_Et<=IsoCut_Progression) {
                        if(iEG==0) {
                            Et_Double = std::min(egEt[iEG+1], (EG_Et-10));
                        }
                        else  Et_Double = EG_Et-10;
                        ptMap_[pt_ProgressionD_]->Fill(Et_Double);
                        Filled_ProgressionD = kTRUE;
                    }
                }

                if(nEGs >=2) {
                    if( (EG_Iso_Et<=IsoCut_Progression )  and (abs(EG_Eta) < 1.5) and (nIsoERXXp10EGs  < 2) ) {

                        if(nIsoERXXp10EGs==0) {
                            doubleIsoXXp10_thresholdER   =  EG_Et ;
                        }
                        else{
                            doubleIsoXXp10_thresholdER = std::min( EG_Et , doubleIsoXXp10_thresholdER-10 ) ;
                        }
                        nIsoERXXp10EGs++;
                    }

                    if( (EG_Iso_Et<=IsoCut_Progression) and  (abs(EG_Eta) < 1.5) and ( nIsoERXXp5EGs < 2 ) ) {
                        if(nIsoERXXp5EGs==0) {
                            doubleIsoXXp5_thresholdER = EG_Et;
                        }
                        else  {
                            doubleIsoXXp5_thresholdER = std::min( EG_Et , doubleIsoXXp5_thresholdER -1 );
                        }

                        nIsoERXXp5EGs+=1;
                    }
                }

                    if( EG_Iso_Et<=IsoCut_Progression  and nIsoEREGs < 2) {
                        
                        if( doubleIso_threshold > EG_Et ) 
                            doubleIso_threshold = EG_Et ;
                        nIsoEGs+=1;

                    }

                    if( EG_Iso_Et<=IsoCut_Progression and abs(EG_Eta) < 1.5  and nIsoEREGs < 2) {
                        if( doubleIso_thresholdER > EG_Et ) 
                            doubleIso_thresholdER = EG_Et ;
                        nIsoEREGs+=1;
                    }

            }

            if ( nIsoEGs >= 2)
            {
                        ptMap_[pt_ProgressionDI_]->Fill(doubleIso_threshold);
            }
            if ( nIsoEREGs >= 2)
            {
                        ptMap_[pt_ProgressionDIER_]->Fill(doubleIso_thresholdER);
            }

            if(nIsoERXXp5EGs > 1 and (doubleIsoXXp5_thresholdER >=0 ) )
            {
                        ptMap_[pt_ProgressionDIXXp5ER_]->Fill(doubleIsoXXp5_thresholdER);
            }

            if( nIsoERXXp10EGs > 1 and ( doubleIsoXXp10_thresholdER >= 0 ) )
            {
                        ptMap_[pt_ProgressionDIXXp10ER_]->Fill(doubleIsoXXp10_thresholdER );
            }
        }
    }

    std::cout<<"total entries: "<<nEntries_<<", sumweight:"<<sum_weight<<std::endl;
    double scale = (11.2456 * bunches)/sum_weight;

    //Filling Rate Histos ( Single EG + Double EG)
    for(UInt_t i=0; i< ET_MAX; i++) {
        for (auto it :lutProgOptVec_) {
            TString CurrentNameHisto   = "pt_Progression" + it;
            TString CurrentNameHisto1  = "rate_Progression" + it;
            TString CurrentNameHistoD  = "pt_Progression_double" + it;
            TString CurrentNameHisto1D = "rate_Progression_double" + it;
            TString CurrentNameHistoDI  = "pt_Progression_doubleIso" + it;
            TString CurrentNameHisto1DI = "rate_Progression_doubleIso" + it;
            TString CurrentNameHistoDIER  = "pt_Progression_doubleIsoER" + it;
            TString CurrentNameHisto1DIER = "rate_Progression_doubleIsoER" + it;
            TString CurrentNameHistoDXXp10IER  = "pt_Progression_doubleIsoXXp10ER" + it;
            TString CurrentNameHisto1DXXp10IER = "rate_Progression_doubleIsoXXp10ER" + it;
            TString CurrentNameHistoDXXp5IER  = "pt_Progression_doubleIsoXXp5ER" + it;
            TString CurrentNameHisto1DXXp5IER = "rate_Progression_doubleIsoXXp5ER" + it;

            rateMap_[CurrentNameHisto1]->SetBinContent(i+1, ptMap_[CurrentNameHisto]->Integral(i+1,ET_MAX)*scale); //Single EG Iso
            rateMap_[CurrentNameHisto1D]->SetBinContent(i+1, ptMap_[CurrentNameHistoD]->Integral(i+1,ET_MAX)*scale); // Double EG Iso
            rateMap_[CurrentNameHisto1DI]->SetBinContent(i+1, ptMap_[CurrentNameHistoDI]->Integral(i+1,ET_MAX)*scale); // Double EG DoubleIso
            rateMap_[CurrentNameHisto1DIER]->SetBinContent(i+1, ptMap_[CurrentNameHistoDIER]->Integral(i+1,ET_MAX)*scale); // Double EG DoubleIso ER
            rateMap_[CurrentNameHisto1DXXp10IER]->SetBinContent(i+1, ptMap_[CurrentNameHistoDXXp10IER]->Integral(i+1,ET_MAX)*scale); // Double EG DoubleIso ER
            rateMap_[CurrentNameHisto1DXXp5IER]->SetBinContent(i+1, ptMap_[CurrentNameHistoDXXp5IER]->Integral(i+1,ET_MAX)*scale); // Double EG DoubleIso ER

        }
    }

    /*  //Finding Et for fixed rate (frate)
    int xbin;
    double frate_=3.9;
    TString CurrentNameHisto1DI = "rate_Progression_double_inclusive";
    float vl,bl; //last value, last bin
    for(UInt_t i=0;i< ET_MAX;i++) {
      float rate = rateMap_[CurrentNameHisto1DI]->GetBinContent(i+1);
      if(rate > frate_) {
        bl =i;
        vl=rate;
      }

        if(rate <=frate_){
          xbin = (rate - frate_) < abs(rate-vl) ? i : bl;
          break;
        }
    }
    std:: cout<<"rate "<<frate_<<" , "<<"et: "<< xbin<<std::endl;
    */

    optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
    optionsFile_->cd("Step2Histos");

    //Event loop5B for turnon
    if (fChain1 == 0) return;
    nEntries1_ = fChain1->GetEntriesFast();
    Long64_t nbytes1 = 0, nb1 = 0;
    double sum=0;
    std::cout << "Total available  Entries For Efficiency " << nEntries1_ << std::endl;
    if(maxEntriesForEfficiency_ >0 ) nEntries1_= maxEntriesForEfficiency_ < nEntries1_ ? maxEntriesForEfficiency_ : nEntries1_;
    std::cout<<"Processing a total of "<<nEntries1_<<" Entries \n";
    t_start = std::chrono::high_resolution_clock::now();

    //  for (Long64_t jentry=0; jentry< 100; jentry++) {
    for (Long64_t jentry=0; jentry< nEntries1_; jentry++) {
        Long64_t ientry = fChain1->LoadTree(jentry);
        if (ientry < 0) break;
        nb1 = fChain1->GetEntry(jentry);
        nbytes1 += nb1;

        if(jentry%reportEvery_ == 0 )
        {
            t_end = std::chrono::high_resolution_clock::now();
            std::cout<<"Processing Entry in event loop (TurnOn) : "<<jentry<<" / "<<nEntries1_<<"  [ "<<100.0*jentry/nEntries1_<<"  % ]  "
                     << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                     <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nEntries1_ - jentry)/(1e-9 + jentry)* 0.001
                     <<std::endl;
        }
        if(!( isProbeLoose==1 && fabs(eleProbeEta) < 2.5  && sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2))>0.6)) continue;
        if(l1tEmuRawEt < 0.) continue;
        pT_all->Fill(eleProbeSclEt);
        sum++;

        std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(short(l1tEmuTowerIEta)));
        if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
        else in_compressediEta = nBinsIEta-1;

        std::map<short, short>::iterator EtPos = lutMapEt.find(short(l1tEmuRawEt));
        if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
        else in_compressediEt = nBinsIEt-1;

        std::map<short, short>::iterator NTTPos = lutMapNTT.find(short(l1tEmuNTT));
        if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
        else in_compressedNTT = nBinsNTT-1;


        Int_t IsoCut_Progression;
        for(UInt_t e = etMin_; e <= etMax_; e += 1) {

                //Filling pt Progression for Turnon
                TString PtPassName_= "pT_pass_option_Et" + std::to_string(e);
                if(l1tEmuPt >= e )	  pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
                for (auto it :lutProgOptVec_) {
                TString ResultProgressionName_= "LUT_Progression_" + it ;
                TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
                IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
                
                PtPassName_= "pT_pass_option" + it + "_Et" + std::to_string(e);
                //std::cout<<__LINE__<<" | "<<"Et : "<<l1tEmuPt<<" l1tEmuIsoEt : "<<l1tEmuIsoEt<<" IsoCut_Progression : "<<IsoCut_Progression<<" [  "<<ResultProgressionName_ <<" ] \n"; 
                if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression)	  pt_pass_Map_[PtPassName_]->Fill(eleProbeSclEt);
                

                //Fillling Nvtx/Eta with Iso histos
                if(eleProbeSclEt>32) {
                    if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression) {
                        TString NvtxPassIso_ ="Nvtx_Pass_Iso_" + it + "_Et_" + std::to_string(e);
                        TString EtaPassIso_  ="Eta_Pass_Iso_" + it + "_Et_" + std::to_string(e);
                        //TString nTTPassIso_  ="nTT_Pass_Iso_" + it + "_Et_" + std::to_string(e);
                        Nvtx_pass_Map_[NvtxPassIso_]->Fill(Nvtx,1);
                        eta_pass_Map_[EtaPassIso_]->Fill(eleProbeEta,1);
                       // nTT_pass_Map_[nTTPassIso_]->Fill(l1tEmuNTT,1);

                    }
                    else {
                        TString NvtxFailIso_ ="Nvtx_Fail_Iso_" + it + "_Et_" + std::to_string(e);
                        TString EtaFailIso_ ="Eta_Fail_Iso_" + it + "_Et_" + std::to_string(e);
                        //TString nTTFailIso_ ="nTT_Fail_Iso_" + it + "_Et_" + std::to_string(e);
                        Nvtx_fail_Map_[NvtxFailIso_]->Fill(Nvtx,1);
                        eta_fail_Map_[EtaFailIso_]->Fill(eleProbeEta,1);
                      //  nTT_fail_Map_[nTTFailIso_]->Fill(l1tEmuNTT,1);
                    }
                }

                if( ( eleProbeSclEt > (e-2) ) and ( eleProbeSclEt < (e+15)  ) )
                {
                    if(l1tEmuPt >= e && l1tEmuIsoEt <= IsoCut_Progression) {
                        TString nTTPassIso_  ="nTT_Pass_Iso_" + it + "_Et_" + std::to_string(e);
                        nTT_pass_Map_[nTTPassIso_]->Fill(l1tEmuNTT,1);

                    }
                    else {
                        TString nTTFailIso_ ="nTT_Fail_Iso_" + it + "_Et_" + std::to_string(e);
                        nTT_fail_Map_[nTTFailIso_]->Fill(l1tEmuNTT,1);
                    }

                }
                
            //Fillling Nvtx/Eta no Iso histos
            if(eleProbeSclEt>32) {
                if(l1tEmuPt >= e) {
                    TString NvtxPass_ ="Nvtx_Pass_Et_" + std::to_string(e);
                    TString EtaPass_ ="Eta_Pass_Et_" + std::to_string(e);
                    TString nTTPass_ ="nTT_Pass_Et_" + std::to_string(e);
                    Nvtx_pass_Map_[NvtxPass_]->Fill(Nvtx,1);
                    eta_pass_Map_[EtaPass_]->Fill(eleProbeEta,1);
                    nTT_pass_Map_[nTTPass_]->Fill(l1tEmuNTT,1);
                }
                else {
                    TString NvtxFail_ ="Nvtx_Fail_Et_" + std::to_string(e);
                    TString EtaFail_ ="Eta_Fail_Et_" + std::to_string(e);
                    TString nTTFail_ ="nTT_Fail_Et_" + std::to_string(e);
                    Nvtx_fail_Map_[NvtxFail_]->Fill(Nvtx,1);
                    eta_fail_Map_[EtaFail_]->Fill(eleProbeEta,1);
                    nTT_fail_Map_[nTTFail_]->Fill(l1tEmuNTT,1);
                }
            }
        }

      }
    } //End of Event loop for turnon



    //Turnons/Nvtx_Eff/Eta_Eff  for various Ets
    for(UInt_t e = etMin_ ; e <= etMax_ ; e += 1) {

            TString PtPassName_= "pT_pass_option_Et" + std::to_string(e);
            TString turnOn_Option_="turnOn_Option_Et_" + std::to_string(e);
            td3->cd();
            turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp");
            turnOn_Map_[turnOn_Option_]->Write();

        for (auto it :lutProgOptVec_) {
            outputFile_->cd();

            TString PtPassName_= "pT_pass_option" + it + "_Et" + std::to_string(e);
            TString turnOn_Option_="turnOn_Option" +it + "_Et_" + std::to_string(e);
            td3->cd();
            turnOn_Map_[turnOn_Option_]=new TGraphAsymmErrors( pt_pass_Map_[PtPassName_],pT_all,"cp");
            turnOn_Map_[turnOn_Option_]->Write();

            double acceptance = (pt_pass_Map_[PtPassName_]->GetEntries())/(pT_all->GetEntries()) ;
            double rate = FindRate(it, e);
            double rateD = FindRateD(it, e);

            th1fStore["EtForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),e);
            th1fStore["AcceptanceForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),acceptance);
            th1fStore["RateForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),rate);
            th1fStore["RateDForTurnons"]->Fill(turnOn_Map_[turnOn_Option_]->GetName(),rateD);

            //Nvtx/Eta Eff with  Iso   + nTT
            TString NvtxPassIso_ ="Nvtx_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString NvtxFailIso_ ="Nvtx_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString NvtxIsoEff_ ="Nvtx_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            createEfficiencyHistograms(Nvtx_pass_Map_[NvtxPassIso_],Nvtx_fail_Map_[NvtxFailIso_],Nvtx_Eff_Map_[NvtxIsoEff_]);

            TString EtaPassIso_ ="Eta_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString EtaFailIso_ ="Eta_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString EtaIsoEff_ ="Eta_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            createEfficiencyHistograms(eta_pass_Map_[EtaPassIso_],eta_fail_Map_[EtaFailIso_],eta_Eff_Map_[EtaIsoEff_]);

            TString nTTPassIso_ ="nTT_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString nTTFailIso_ ="nTT_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString nTTIsoEff_ ="nTT_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            createEfficiencyHistograms(nTT_pass_Map_[nTTPassIso_],nTT_fail_Map_[nTTFailIso_],nTT_Eff_Map_[nTTIsoEff_]);

        }

        //Nvtx/Eta without Iso + nTT
        TString NvtxPass_ ="Nvtx_Pass_Et_" + std::to_string(e);
        TString NvtxFail_ ="Nvtx_Fail_Et_" + std::to_string(e);
        TString NvtxEff_ ="Nvtx_Eff_Et_" + std::to_string(e);
        createEfficiencyHistograms(Nvtx_pass_Map_[NvtxPass_],Nvtx_fail_Map_[NvtxFail_],Nvtx_Eff_Map_[NvtxEff_]);

        TString EtaPass_ ="Eta_Pass_Et_" + std::to_string(e);
        TString EtaFail_ ="Eta_Fail_Et_" + std::to_string(e);
        TString EtaEff_ ="Eta_Eff_Et_" + std::to_string(e);
        createEfficiencyHistograms(eta_pass_Map_[EtaPass_],eta_fail_Map_[EtaFail_],eta_Eff_Map_[EtaEff_]);

        TString nTTPass_ ="nTT_Pass_Et_" + std::to_string(e);
        TString nTTFail_ ="nTT_Fail_Et_" + std::to_string(e);
        TString nTTEff_ ="nTT_Eff_Et_" + std::to_string(e);
        createEfficiencyHistograms(nTT_pass_Map_[nTTPass_],nTT_fail_Map_[nTTFail_],nTT_Eff_Map_[nTTEff_]);

    }
}



void ApplyIsolation::bookHistograms() {
    outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
    outputFile_->cd();

    for (auto it : lutProgOptVec_) {
        //For Single EG pt Progression with Iso
        TString  CurrentNameHisto = "pt_Progression" + it;
        if(!check_pt_rate_dir) {
            td = outputFile_->mkdir("pt_progression_for_rate");
            check_pt_rate_dir = true;
        }
        td->cd();
        TH1F* pt_Progression=  new TH1F(CurrentNameHisto, CurrentNameHisto,ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHisto,pt_Progression));

        //For Double EG pt Progression with iso
        if(!check_pt_rate_double_dir) {
            tdD = outputFile_->mkdir("pt_progression_for_rate_double");
            check_pt_rate_double_dir = true;
        }
        tdD->cd();
        TString  CurrentNameHistoD = "pt_Progression_double" + it;
        TH1F* pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD , ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));

        CurrentNameHistoD = "pt_Progression_doubleIso" + it;
        pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD,ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));

        CurrentNameHistoD = "pt_Progression_doubleIsoER" + it;
        pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD,ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));

        CurrentNameHistoD = "pt_Progression_doubleIsoXXp10ER" + it;
        pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD,ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));

        CurrentNameHistoD = "pt_Progression_doubleIsoXXp5ER" + it;
        pt_ProgressionD=  new TH1F(CurrentNameHistoD, CurrentNameHistoD,ET_MAX,  0.0, ET_MAX);
        ptMap_.insert(std::make_pair(CurrentNameHistoD,pt_ProgressionD));


        //For Single EG Rate with Iso
        if(!check_rate_dir) {
            td1 = outputFile_->mkdir("rate_progression");
            check_rate_dir = true;
        }
        td1->cd();
        TString CurrentNameHisto1= "rate_Progression" + it;
        TH1F* rate_Progression = new TH1F(CurrentNameHisto1, CurrentNameHisto1, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1,rate_Progression));
        
        //For Double EG Rate with Iso
        if(!check_rate_double_dir) {
            td1D = outputFile_->mkdir("rate_progression_double");
            check_rate_double_dir = true;
        }
        td1D->cd();
        
        TString CurrentNameHisto1D= "rate_Progression_double" + it;
        TH1F* rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));

        CurrentNameHisto1D= "rate_Progression_doubleIso" + it;
        rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));

        CurrentNameHisto1D= "rate_Progression_doubleIsoER" + it;
        rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));

        CurrentNameHisto1D= "rate_Progression_doubleIsoXXp10ER" + it;
        rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));

        CurrentNameHisto1D= "rate_Progression_doubleIsoXXp5ER" + it;
        rate_ProgressionD = new TH1F(CurrentNameHisto1D, CurrentNameHisto1D, ET_MAX,  0.0-0.5, ET_MAX-0.5 );
        rateMap_.insert(std::make_pair(CurrentNameHisto1D,rate_ProgressionD));


        for(UInt_t e = etMin_ ; e <= etMax_ ; e += 1) {
            //For pt Progression for turnons
            TString  PtPassName_ = "pT_pass_option" + it + "_Et" + std::to_string(e);
            if(!check_pt_turn_on_dir) {
                td2 = outputFile_->mkdir("pt_progression_for_turnon");
                check_pt_turn_on_dir = true;
            }
            td2->cd();
            TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
            pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));

            //For TurnOns
            TString turnOn_Option_="turnOn_Option" + it + "_Et_" + std::to_string(e);
            TGraphAsymmErrors* turnOn_Option;
            if(!check_turn_on_dir) {
                td3 = outputFile_->mkdir("turn_on_progression");
                check_turn_on_dir = true;
            }
            turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
        }
    }

       for(UInt_t e = etMin_ ; e <= etMax_ ; e += 1) {

            //For pt Progression for turnons
            TString  PtPassName_ = "pT_pass_option_Et" + std::to_string(e);
            if(!check_pt_turn_on_dir) {
                td2 = outputFile_->mkdir("pt_progression_for_turnon");
                check_pt_turn_on_dir = true;
            }
            td2->cd();
            TH1F*PtPass= new TH1F(PtPassName_, PtPassName_, nBins_fine,xEdges_fine);
            pt_pass_Map_.insert(std::make_pair(PtPassName_, PtPass));

            //For TurnOns
            TString turnOn_Option_="turnOn_Option_Et_" + std::to_string(e);
            TGraphAsymmErrors* turnOn_Option;
            if(!check_turn_on_dir) {
                td3 = outputFile_->mkdir("turn_on_progression");
                check_turn_on_dir = true;
            }
            turnOn_Map_.insert(std::make_pair(turnOn_Option_,turnOn_Option));
        }

    td2->cd();
    pT_all = new TH1F("pT_all","pT_all",nBins_fine,xEdges_fine); //pt all for turnon

    outputFile_->cd();
    th1fStore["EtForTurnons"]   = new TH1F("EtForTurnons","",3,0.0,3.0);
    th1fStore["EtForTurnons"]->SetCanExtend(TH1::kAllAxes);
    th1fStore["AcceptanceForTurnons"]   = new TH1F("AcceptanceForTurnons","",3,0.0,3.0);
    th1fStore["AcceptanceForTurnons"]->SetCanExtend(TH1::kAllAxes);
    th1fStore["RateForTurnons"]   = new TH1F("RateForTurnons","",3,0.0,3.0);
    th1fStore["RateForTurnons"]->SetCanExtend(TH1::kAllAxes);
    th1fStore["RateDForTurnons"]   = new TH1F("RateDForTurnons","",3,0.0,3.0);
    th1fStore["RateDForTurnons"]->SetCanExtend(TH1::kAllAxes);

    //For Vertex/Eta + nTT
    for(UInt_t e = etMin_; e <= etMax_ ; e += 1) {

        if(!check_Nvtx) {
            td5 = outputFile_->mkdir("Nvtx");
            check_Nvtx = true;
        }
        if(!check_Nvtx_Iso) {
            td6 = outputFile_->mkdir("Nvtx_Iso");
            check_Nvtx_Iso = true;
        }

        if(!check_eta) {
            td7 = outputFile_->mkdir("Eta");
            check_eta = true;
        }
        td7->cd();
        if(!check_eta_Iso) {
            td8 = outputFile_->mkdir("Eta_Iso");
            check_eta_Iso = true;
        }
        td8->cd();
        if(!check_nTT) {
            td9 = outputFile_->mkdir("nTT");
            check_nTT = true;
        }
        td9->cd();
        if(!check_nTT_Iso) {
            td10 = outputFile_->mkdir("nTT_Iso");
            check_nTT_Iso = true;
        }
        td10->cd();

        TString NvtxPass_ ="Nvtx_Pass_Et_" + std::to_string(e);
        TString NvtxFail_ ="Nvtx_Fail_Et_" + std::to_string(e);
        TString NvtxEff_  ="Nvtx_Eff_Et_" + std::to_string(e);
        td5->cd();
        TH1D*thNvtxP = new TH1D(NvtxPass_, NvtxPass_, 100, -0.5, 99.5);
        TH1D*thNvtxF = new TH1D(NvtxFail_, NvtxFail_, 100, -0.5, 99.5);
        TH1D*thNvtxE = new TH1D(NvtxEff_, NvtxEff_, 100, -0.5, 99.5);
        Nvtx_pass_Map_.insert(std::make_pair(NvtxPass_, thNvtxP));
        Nvtx_fail_Map_.insert(std::make_pair(NvtxFail_,thNvtxF));
        Nvtx_Eff_Map_.insert(std::make_pair(NvtxEff_,thNvtxE));

        TString EtaPass_ ="Eta_Pass_Et_" + std::to_string(e);
        TString EtaFail_ ="Eta_Fail_Et_" + std::to_string(e);
        TString EtaEff_ ="Eta_Eff_Et_" + std::to_string(e);
        td7->cd();
        TH1D*thetaP = new TH1D(EtaPass_, EtaPass_, 100, -2.5, 2.5);
        TH1D*thetaF = new TH1D(EtaFail_, EtaFail_, 100, -2.5, 2.5);
        TH1D*thetaE = new TH1D(EtaEff_, EtaEff_, 100, -2.5, 2.5);
        eta_pass_Map_.insert(std::make_pair(EtaPass_, thetaP));
        eta_fail_Map_.insert(std::make_pair(EtaFail_,thetaF));
        eta_Eff_Map_.insert(std::make_pair(EtaEff_,thetaE));

        TString nTTPass_ ="nTT_Pass_Et_" + std::to_string(e);
        TString nTTFail_ ="nTT_Fail_Et_" + std::to_string(e);
        TString nTTEff_ ="nTT_Eff_Et_" + std::to_string(e);
        td9->cd();
        TH1D*thnTTP = new TH1D(nTTPass_, nTTPass_, 260, -1, 259);
        TH1D*thnTTF = new TH1D(nTTFail_, nTTFail_, 260, -1, 259);
        TH1D*thnTTE = new TH1D(nTTEff_, nTTEff_, 260, -1, 259);
        nTT_pass_Map_.insert(std::make_pair(nTTPass_, thnTTP));
        nTT_fail_Map_.insert(std::make_pair(nTTFail_,thnTTF));
        nTT_Eff_Map_.insert(std::make_pair(nTTEff_,thnTTE));

        for (auto it :lutProgOptVec_) {    // with Iso
            TString NvtxPassIso_ ="Nvtx_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString NvtxFailIso_ ="Nvtx_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString NvtxIsoEff_ ="Nvtx_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            td6->cd();
            TH1D*thNvtxPI = new TH1D(NvtxPassIso_, NvtxPassIso_, 100, -0.5, 99.5);
            TH1D*thNvtxFI = new TH1D(NvtxFailIso_, NvtxFailIso_, 100, -0.5, 99.5);
            TH1D*thNvtxEI = new TH1D(NvtxIsoEff_, NvtxIsoEff_, 100, -0.5, 99.5);
            Nvtx_pass_Map_.insert(std::make_pair(NvtxPassIso_,thNvtxPI));
            Nvtx_fail_Map_.insert(std::make_pair(NvtxFailIso_,thNvtxFI));
            Nvtx_Eff_Map_.insert(std::make_pair(NvtxIsoEff_,thNvtxEI));

            TString EtaPassIso_ ="Eta_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString EtaFailIso_ ="Eta_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString EtaIsoEff_ ="Eta_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            td8->cd();
            TH1D*thetaPI = new TH1D(EtaPassIso_, EtaPassIso_, 100, -2.5, 2.5);
            TH1D*thetaFI = new TH1D(EtaFailIso_, EtaFailIso_, 100, -2.5, 2.5);
            TH1D*thetaEI = new TH1D(EtaIsoEff_, EtaIsoEff_, 100, -2.5, 2.5);
            eta_pass_Map_.insert(std::make_pair(EtaPassIso_,thetaPI));
            eta_fail_Map_.insert(std::make_pair(EtaFailIso_,thetaFI));
            eta_Eff_Map_.insert(std::make_pair(EtaIsoEff_,thetaEI));

            TString nTTPassIso_ ="nTT_Pass_Iso_" + it + "_Et_" + std::to_string(e);
            TString nTTFailIso_ ="nTT_Fail_Iso_" + it + "_Et_" + std::to_string(e);
            TString nTTIsoEff_ ="nTT_Eff_Iso_" + it + "_Et_" + std::to_string(e);
            td10->cd();
            TH1D*thnTTPI = new TH1D(nTTPassIso_, nTTPassIso_, 260, -1, 259);
            TH1D*thnTTFI = new TH1D(nTTFailIso_, nTTFailIso_, 260, -1, 259);
            TH1D*thnTTEI = new TH1D(nTTIsoEff_, nTTIsoEff_, 260, -1, 259);
            nTT_pass_Map_.insert(std::make_pair(nTTPassIso_,thnTTPI));
            nTT_fail_Map_.insert(std::make_pair(nTTFailIso_,thnTTFI));
            nTT_Eff_Map_.insert(std::make_pair(nTTIsoEff_,thnTTEI));

        }
    }

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
            else if(key=="WeightFileName") weightFileName_ = value.c_str();
            else if (key=="OptionsFileName") optionsFileName_ = value.c_str();
            else if (key=="OutputFileName")  outputFileName_ = value.c_str();
            else if (key=="EtLUTFileName")    readLUTTable(value,nBinsIEt, lutMapEt);
            else if (key=="EtaLUTFileName")   readLUTTable(value,nBinsIEta, lutMapEta);
            else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
            else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());
            else if (key=="ReportEvery")    reportEvery_ = atoi(value.c_str());
            else if (key=="MaxEntriesForEfficency") maxEntriesForEfficiency_ = atoi(value.c_str());
            else if (key=="MaxEntriesForRate")   maxEntriesForRate_ = atoi(value.c_str());
            else if (key=="EtMin") etMin_ = std::stoi(value.c_str());
            else if (key=="EtMax") etMax_ = std::stoi(value.c_str());
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
    fChain->SetBranchAddress("egEta", &egEta);
    fChain->SetBranchAddress("egTowerIEta", &egTowerIEta);
    fChain->SetBranchAddress("egNTT", &egNTT);
    fChain->SetBranchAddress("egBx", &egBx);
    fChain->SetBranchAddress("egIsoEt", &egIsoEt);
    fChain->SetBranchAddress("egRawEt", &egRawEt);
    fChain_1->SetMakeClass(1);
    fChain_1->SetBranchAddress("run",&run);
    fChain_1->SetBranchAddress("nPV_True", &nPV_True);

    //for turnon
    fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
    fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
    fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
    fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
    fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
    fChain1->SetBranchAddress("l1tEmuIsoEt",&l1tEmuIsoEt);

    fChain1->SetBranchAddress("eleProbeEta",&eleProbeEta);
    fChain1->SetBranchAddress("eleProbePhi",&eleProbePhi);
    fChain1->SetBranchAddress("eleTagEta",&eleTagEta);
    fChain1->SetBranchAddress("eleTagPhi",&eleTagPhi);
    fChain1->SetBranchAddress("isProbeLoose",&isProbeLoose);
    fChain1->SetBranchAddress("Nvtx",&Nvtx);

}

double ApplyIsolation::BetaInverse(double x,double p, double q)
{
    double result(0.0);
    double dy = 0.001;
    double eMin = 100;
    for(int i=0; i<1000; i++) {
        double y = i*dy;
        double e = fabs(TMath::BetaIncomplete(y,p,q)-x);
        if (e<eMin)
        {
            eMin = e;
            result = y;
        }
    }
    return result;
}

double ApplyIsolation::FindRate(std::string it, UInt_t e) {
    TString CurrentNameHisto1= "rate_Progression";
    CurrentNameHisto1 += it;
    double rate = rateMap_[CurrentNameHisto1]->GetBinContent(e);
    return rate;
}

double ApplyIsolation::FindRateD(std::string it, UInt_t e) {
    TString CurrentNameHisto1D= "rate_Progression_double";
    CurrentNameHisto1D += it;
    double rateD =  rateMap_[CurrentNameHisto1D]->GetBinContent(e);
    return rateD;
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


