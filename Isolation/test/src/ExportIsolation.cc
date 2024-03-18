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

    if (ntupleFileNameTurnOn_.size() == 0) {
        std::cout << " Inputfile for turn on list missing !!!" << ntupleFileNameTurnOn_ << std::endl;
        return;
    }

    accessTree(ntupleFileNameRate_, ntupleFileNameTurnOn_);

    assert(fChain1);

    bookedHistograms_ = false;

}

ApplyIsolation::~ApplyIsolation() {
}

void ApplyIsolation::accessTree(std::string & input_filelist_rate, std::string & input_filelist_turn_on) {
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
    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

    optionsFile_ = new TFile(optionsFileName_.c_str(), "READ");
    optionsFile_->cd("Step2Histos");

    if (fChain1 == 0) return;
    nEntries1_ = fChain1->GetEntriesFast();
    
    TFile outputTreeFile(outputFileName_.data(),"RECREATE");
    outputTreeFile.cd();
    auto outputTree = fChain1->CloneTree(0);

    std::map<TString,float> isolationValues;
    std::map<TString,float> isolationThresolds;
    
    for (auto it :lutProgOptVec_) {
        
        TString ResultProgressionName_= "ISO_LUT_Progression_" + it ;
        isolationValues[ResultProgressionName_] = 0.0 ;
        isolationThresolds[ResultProgressionName_] = 0 ;
        outputTree->Branch(ResultProgressionName_ , &isolationValues[ResultProgressionName_] ,ResultProgressionName_+"/F") ;
        outputTree->Branch(ResultProgressionName_+"_threashold" , &isolationThresolds[ResultProgressionName_] ,ResultProgressionName_+"_threashold/F") ;

    }
        outputTree->Branch("in_compressediEta" , &in_compressediEta ) ;
        outputTree->Branch("in_compressediEt" , &in_compressediEt ) ;
        outputTree->Branch("in_compressedNTT" , &in_compressedNTT) ;

    Int_t IsoCut_Progression;
    Long64_t nbytes1 = 0, nb1 = 0;
    std::cout << "Total available  Entries For Efficiency " << nEntries1_ << std::endl;
    if(maxEntriesForEfficiency_ >0 ) nEntries1_= maxEntriesForEfficiency_ < nEntries1_ ? maxEntriesForEfficiency_ : nEntries1_;
    std::cout<<"Processing a total of "<<nEntries1_<<" Entries \n";
    t_start = std::chrono::high_resolution_clock::now();

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

        std::map<short, short>::iterator EtaPos = lutMapEta.find(abs(l1tEmuTowerIEta));
        if (EtaPos != lutMapEta.end()) in_compressediEta = EtaPos->second;
        else in_compressediEta = nBinsIEta-1;

        std::map<short, short>::iterator EtPos = lutMapEt.find(l1tEmuRawEt);
        if (EtPos != lutMapEt.end()) in_compressediEt = EtPos->second;
        else in_compressediEt = nBinsIEt-1;

        std::map<short, short>::iterator NTTPos = lutMapNTT.find(l1tEmuNTT);
        if (NTTPos != lutMapNTT.end()) in_compressedNTT = NTTPos->second;
        else in_compressedNTT = nBinsNTT-1;


            for (auto it :lutProgOptVec_) {
                TString ProgressionName_="ISO_LUT_Progression_" + it ;
                TString ResultProgressionName_="Step2Histos/LUT_Progression_" + it ;
                //std::cout<<ResultProgressionName_<<"\n";
                //TH3F* ResultProgressionName = (TH3F*)gDirectory->Get(ResultProgressionName_.Data());
                TH3F* ResultProgressionName = (TH3F*)optionsFile_->Get(ResultProgressionName_.Data());
                IsoCut_Progression = ResultProgressionName->GetBinContent(in_compressediEta+1,in_compressediEt+1,in_compressedNTT+1);
                isolationThresolds[ProgressionName_] = IsoCut_Progression ;
                if(l1tEmuIsoEt <= IsoCut_Progression)	  isolationValues[ProgressionName_] = 1.0  ;
                else                                      isolationValues[ProgressionName_] = 0.0  ; 
                
          }

          outputTree->Fill();

    } //End of Event loop for turnon
    outputTreeFile.cd(); 
    outputTree->Write() ;
    outputTreeFile.Purge();
    outputTreeFile.Write();
    outputTreeFile.Close();
}



void ApplyIsolation::bookHistograms() 
{
        // pass
        return;
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

            if(key=="NtupleFileNameTurnOn")        ntupleFileNameTurnOn_= value;
            else if(key=="WeightFileName") weightFileName_ = value.c_str();
            else if (key=="OptionsFileName") optionsFileName_ = value.c_str();
            else if (key=="OutputFileName")  outputFileName_ = value.c_str();
            else if (key=="EtLUTFileName")    readLUTTable(value,nBinsIEt, lutMapEt);
            else if (key=="EtaLUTFileName")   readLUTTable(value,nBinsIEta, lutMapEta);
            else if (key=="NTTLUTFileName")   readLUTTable(value,nBinsNTT, lutMapNTT);
            else if (key=="nTTRange") nTTRange_ = std::stoi(value.c_str());
            else if (key=="ReportEvery")    reportEvery_ = atoi(value.c_str());
            else if (key=="MaxEntriesForEfficency") maxEntriesForEfficiency_ = atoi(value.c_str());
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
    //for turnon
    fChain1->SetBranchAddress("l1tEmuPt", &l1tEmuPt);
    fChain1->SetBranchAddress("l1tEmuNTT",&l1tEmuNTT);
    fChain1->SetBranchAddress("l1tEmuRawEt",&l1tEmuRawEt);
    fChain1->SetBranchAddress("l1tEmuTowerIEta",&l1tEmuTowerIEta);
    fChain1->SetBranchAddress("eleProbeSclEt",&eleProbeSclEt);
    fChain1->SetBranchAddress("l1tEmuIso",&l1tEmuIso);
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
    std::cout << " Calling Loop" << std::endl;
    treeReader.loops();
 //   treeReader.saveHistograms();
    return 0;
}


