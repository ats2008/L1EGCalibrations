#include "Util.h"
#include <iostream>
#include <fstream>
#include <string>

using namespace std;

void processOptionFile(TString fileName,std::fstream &result_file,Bool_t *filled_Baseline,Double_t baselineEt=24,TString prefix="",Bool_t chkQuality=false,Bool_t saveGraph=false) {
  std::cout<<"fileName : "<<fileName<<"\n";
  std::cout<<"prefix : "<<prefix<<"\n";
  Bool_t temp = *filled_Baseline;
  std::cout<<"bool baseline : "<<temp<<"\n";
  

  if(!temp) result_file <<"fileName , option , eTMin, eff , eTMax, Et_threshold,area, isBetterThanLoose, isBetterThanTight, Acceptance\n";
  
  TString folderName;
  Double_t eT_threshold,left_DX,right_DX,area,acceptance;
  string cmd;
  Bool_t isBetterThanLoose(false);
  Bool_t isBetterThanTight(false);
  int thresoldPassEvents(0),betterThanPassEvents(0),optBins(0);
  eT_threshold=30;
  left_DX=7.0;
  right_DX=7.0;

  Int_t isGoodChkBeg(29),isGoodChkEnd(54);
  
  TFile * f=TFile::Open(fileName,"READ");
  if (!f) std::cout<<"Input file ["<<fileName<<"] not found"<<std::endl;

  //Working with baseline
  TGraphAsymmErrors * baselineTurnOn=nullptr;
  TString baselinePlotName="divide_pt_pass_Et_"+to_string(int(baselineEt))+"_tight_by_pT_all";
  baselineTurnOn=(TGraphAsymmErrors *)f->Get("turnon_progression_Run2/"+baselinePlotName);
  if(!baselineTurnOn) {
    std::cout<<"Baseline plot not found !! ["<<baselinePlotName<<"] from "<<f<<"\n";
  }  
  else { 
    if(temp==false) {
      folderName=prefix+"/Baseline_Run2/";
      cmd="mkdir -p "+folderName;
      system(cmd.c_str());
      if(temp==false) {
	     
         folderName=prefix+"/Baseline_Run2/";
	     cmd="mkdir -p "+folderName;
	     system(cmd.c_str());
	     
	     std::cout<<" Obtained baseline as  : "<<baselineTurnOn->GetName()<<"\n";
	     eT_threshold=baselineEt;
	     area=getIntegral(baselineTurnOn, eT_threshold,eT_threshold-left_DX,eT_threshold+right_DX,folderName,nullptr,saveGraph);
	     std::cout<<"Baseline Area = "<<area<<"\n";

	     auto acceptanceHistogram_Run2 = (TH1F*) f->Get("Run2Turnons_TightIso_Acceptance");
	     acceptance =acceptanceHistogram_Run2->GetBinContent(1);
	     std::cout<<"Acceptance for Baseline = "<<acceptance<<std::endl;
      
      std::cout<<"XX"<<" , "<<"Baseline "
	       <<" ,  "
	       <<" XX "<<" , "
	       <<" XX "<<" , "
	       <<" XX "<<"  "
	       <<" , "<<eT_threshold 
	       <<" , "<<area
	       <<" , "<<" X "
	       <<" , "<<" X "
	       <<" , "<<acceptance
	       <<"\n";

      result_file <<fileName<<" "<<"Baseline "
	       <<"\t"
	       <<" XX "<<"\t"
	       <<" XX "<<"\t"
	       <<" XX "<<"\t"
	       <<"\t"<<eT_threshold 
	       <<"\t"<<area
	       <<"\t"<<" X "
	       <<"\t"<<" X "
	       <<"\t"<<acceptance
	       <<"\n";


    }
        
    }
    auto mindiff=0.2;
    for(Int_t ii=0;ii<=baselineTurnOn->GetN();ii++)
    {
           if(abs(baselineTurnOn->GetPointY(ii) - 0.8) < mindiff ) 
           { 
               isGoodChkBeg=ii  ;
               mindiff=abs(baselineTurnOn->GetPointY(ii) - 0.8);
           }
    }
    
    mindiff=0.2;
    for(Int_t ii=0;ii<=baselineTurnOn->GetN();ii++)
    {
        if(abs(baselineTurnOn->GetPointY(ii) - 0.95) < mindiff ) 
           { 
               isGoodChkEnd=ii  ;
               mindiff=abs(baselineTurnOn->GetPointY(ii) - 0.95);
           }
    }

    *filled_Baseline=true;
  }

  std::cout<<"Setting up isGoodChk Rage = [ "<<isGoodChkBeg<<" , "<<isGoodChkEnd<<" ]"<<"\n";
  
  
  //working with all the option in a file
  TGraphAsymmErrors * graphToIntegrate(nullptr);
  auto descriptionHistogram = (TH1F*) f->Get("FixedRateTurnons");
  auto acceptanceHistogram = (TH1F*) f->Get("FixedRateTurnons_Acceptance");
  // Getting the fixed Rate Metrics Printed
  auto xAxis = descriptionHistogram->GetXaxis();
  auto nBins = descriptionHistogram->GetNbinsX();
  std::vector<std::string> tokens;
  std::cout<<"idx , option , eTMin, eff , eTMax,eTThreshold ,area, isBetterThanLoose, isBetterThanTight, Acceptance\n";
  Bool_t passedThresoldCriteria;
  for(Int_t i=1;i<=nBins;i++)
    {
      isBetterThanTight=false;
      passedThresoldCriteria=false;
      TString histname(xAxis->GetBinLabel(i));
      std::string hname(xAxis->GetBinLabel(i));
      if(hname.size()==0) continue;
      optBins++;
      graphToIntegrate=(TGraphAsymmErrors *) f->Get("turn_on_progression/"+histname);
      if(! graphToIntegrate) {
	std::cout<<"\n\tturn_on_progression/"+histname<<" not available "<<"\n";
	continue;
      }
      tokens.clear();
      tokenize(hname,tokens,"_");
      std::replace( tokens[4].begin(), tokens[4].end(), 'p', '.');
      std::replace( tokens[5].begin(), tokens[5].end(), 'p', '.');
      std::replace( tokens[6].begin(), tokens[6].end(), 'p', '.');
      
      acceptance=acceptanceHistogram->GetBinContent(i);
      eT_threshold=descriptionHistogram->GetBinContent(i);
      if(eT_threshold <= baselineEt+2.0 ) passedThresoldCriteria=true;
      if(chkQuality and eT_threshold > baselineEt+2.0 ) continue;
      if(passedThresoldCriteria) thresoldPassEvents++;
      
      if(baselineTurnOn)  isBetterThanTight=isGoodTurnON(baselineTurnOn,graphToIntegrate,29,54);
      if(chkQuality and !isBetterThanTight) continue;
      if(isBetterThanTight) betterThanPassEvents++;
      
         folderName=prefix+"/"+tokens[3]+"/";
          if(saveGraph  or (isBetterThanTight and passedThresoldCriteria)) {
                cmd="mkdir -p "+folderName;
                 system(cmd.c_str());
            }
      area=getIntegral(graphToIntegrate, eT_threshold,eT_threshold-left_DX,eT_threshold+right_DX,folderName,baselineTurnOn,saveGraph or (isBetterThanTight and passedThresoldCriteria) );
      std::cout<<i<<" , "<<tokens[3]
	       <<" ,  "
	       <<tokens[4]<<" , "
	       <<tokens[5]<<" , "
	       <<tokens[6]<<"  "
	       <<" , "<<eT_threshold 
	       <<" , "<<area
	       <<" , "<<isBetterThanLoose
	       <<" , "<<isBetterThanTight
	       <<" , "<<acceptance
	       <<"\n";

      result_file <<fileName<<"\t"<<tokens[3]
		  <<"\t"
		  <<tokens[4]<<"\t"
		  <<tokens[5]<<"\t"
		  <<tokens[6]<<"\t"
		  <<"\t"<<eT_threshold
		  <<"\t"<<area
		  <<"\t"<<isBetterThanLoose
		  <<"\t"<<isBetterThanTight
		  <<"\t"<<acceptance
		  <<"\n";
    }
  
  f->Close();
  std::cout<<" Number of options passing the \" eT thresold \" : "<<thresoldPassEvents<<" / "<<optBins<<"\n";
  std::cout<<" Number of options passing the \" better than the baseline \" : "<<betterThanPassEvents<<" / "<<optBins<<"\n";
}

int  main(int argc,char *argv[])
{
  
  if(argc<2)
    {
      std::cout<<" Usage : \n"
	       <<" ./optEval.exe listof<fname>  prefix_to_save\n"
	       <<"\n";
      
      exit(1);
    }
  std::fstream result_file; 
  TString basePath(argv[2]);
  Bool_t filled_Baseline = false;
  std::ifstream file_names(argv[1]);
  std::string line;
  Bool_t saveGraph,chkQuality;
  TString savePrefix("");
  chkQuality=atoi(argv[2]);
  saveGraph=atoi(argv[3]);
  
 // if(saveGraph)
  {savePrefix=argv[4];}
  
  result_file.open(savePrefix+"/result.txt",ios::out);  
  std::cout<<" Processing Only Good Graphs  : "<<chkQuality<<"\n";
  std::cout<<" Saving Graphs : "<<saveGraph<<"\n";
  std::cout<<" Saving to : "<<savePrefix<<"\n";

  if(file_names.is_open()) {
    while (file_names>>line) {
      std::cout<<"Processing file : "<<line<<std::endl;
      processOptionFile(line,result_file,& filled_Baseline,27.0,savePrefix,chkQuality,saveGraph);
    }
  }
  std::cout<<"Results saved to file : "<<savePrefix+"/result.txt"<<"\n";
  result_file.close();
  return 0;
}



