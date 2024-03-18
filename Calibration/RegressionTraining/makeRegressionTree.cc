#include <TTree.h>
#include <TString.h>
#include <TFile.h>
#include <TChain.h>
#include <TH1F.h>

//#include "Helpers.C"

using namespace std;

#ifdef __MAKECINT__
#pragma link C++ class vector<float>+;
#endif


void makeRegressionTree(){ //vector<TString> input_names, TString output_name){

/*
  TChain* intree = new TChain("produceNtuple/eIDSimpleTree");
  for(unsigned int i=0; i<input_names.size();i++)
    intree->Add(input_names[i]);
*/

  TFile *f=TFile::Open("/grid_mnt/t3storage3/athachay/l1egamma/data/2024/TandP_EGamma0And1_2023D_reEmulated_v0.root");
  TTree *intree=(TTree*)f->Get("Ntuplizer/TagAndProbe");


  TString output_file("./regressionTrainerFile.root");
  TFile* f_new = TFile::Open(output_file,"RECREATE");

  TTree* tree = new TTree("eIDSimpleTree","eIDSimpleTree");

  //Branches new tree

  ULong64_t _event;
  Int_t _lumi;
  Int_t _run;  
  int _ieta;
  int _E;
  int _shape;
  float _target;
  float _target2;
  int _HoverERatio;
  int _nTT;
  int _Run2IdLevel;
  Float_t eleProbePt;
 
  //Branches in intree
  ULong64_t _nEvent;
  Int_t _nRun;
  Int_t _nLumi;
  float _ele_probePt;//[10];
  float _ele_sclEt;//[10];
  int _ele_L1Stage2_emul_ieta;//[10];
  int _ele_L1Stage2_emul_rawEt;//[10];
  int _ele_L1Stage2_emul_shape;//[10];
  int _ele_L1Stage2_emul_hOverERatio;//[10];
  int _ele_L1Stage2_emul_nTT;//[10];
  int _ele_L1Stage2_emul_hwIsoSum;//[10];
  int _ele_LooseIdDecisions;//[10];
 
  intree->SetBranchStatus("*",0);
  
  intree->SetBranchStatus("EventNumber",1);
  intree->SetBranchStatus("RunNumber",1);
  intree->SetBranchStatus("lumi",1);
  intree->SetBranchStatus("eleProbeSclEt",1);
  intree->SetBranchStatus("l1tEmuTowerIEta",1);
  intree->SetBranchStatus("l1tEmuRawEt",1);
  intree->SetBranchStatus("shape",1);
  intree->SetBranchStatus("TowerHoE",1);
  intree->SetBranchStatus("l1tEmuNTT",1);
  intree->SetBranchStatus("isProbeLoose",1);
  intree->SetBranchStatus("eleProbePt",1);

  intree->SetBranchAddress("EventNumber",&_nEvent);
  intree->SetBranchAddress("RunNumber",&_nRun);
  intree->SetBranchAddress("lumi",&_nLumi);
  intree->SetBranchAddress("eleProbePt",&_ele_probePt);
  intree->SetBranchAddress("eleProbeSclEt",&_ele_sclEt);
  intree->SetBranchAddress("l1tEmuTowerIEta",&_ele_L1Stage2_emul_ieta);
  intree->SetBranchAddress("l1tEmuRawEt",&_ele_L1Stage2_emul_rawEt);
  intree->SetBranchAddress("shape",&_ele_L1Stage2_emul_shape);
  intree->SetBranchAddress("TowerHoE",&_ele_L1Stage2_emul_hOverERatio);
  intree->SetBranchAddress("l1tEmuNTT",&_ele_L1Stage2_emul_nTT);
  intree->SetBranchAddress("isProbeLoose",&_ele_LooseIdDecisions);
  intree->SetBranchAddress("eleProbePt",&eleProbePt);

  tree->Branch("Run",&_run,"_run/I");
  tree->Branch("Event",&_event,"_event/l");
  tree->Branch("Lumi",&_lumi,"_lumi/I");
  tree->Branch("ieta",&_ieta,"_ieta/I");
  tree->Branch("E",&_E,"_E/I");
  tree->Branch("eleProbePt",&_ele_probePt,"_ele_probePt/F");
  tree->Branch("eleProbeSclEt",&_ele_sclEt,"_ele_sclEt/F");
  tree->Branch("shape",&_shape,"_shape/I");
  tree->Branch("target",&_target,"_target/F");
  tree->Branch("target2",&_target,"_target2/F");
  tree->Branch("HoverERatio",&_HoverERatio,"_HoverERatio/I");
  tree->Branch("nTT",&_nTT,"_nTT/I");
  tree->Branch("eleProbePt",&eleProbePt,"eleProbePt/F");
  tree->Branch("Run2IdLevel",&_Run2IdLevel,"_Run2IdLevel/I");


    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

  Long64_t nentries = intree->GetEntries();
  Long64_t nentries_beg = 0;
  std::cout<<" Available total "<<nentries<<" \n";
  std::cout<<" Processing total "<<nentries - nentries_beg<<" \n";

  
  for(int i=0;i<nentries;i++){
    if(i<nentries_beg) continue;
    if(i%50000==0){
         t_end = std::chrono::high_resolution_clock::now();
         cout<<"Processing Entry "<<i<<" / "<<nentries<<"  [ "<<100.0*i/nentries<<"  % ]  "
             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries  - i )/( i - nentries_beg) * 0.001
             <<endl;
    }

    intree->GetEntry(i);
    
      _run = 0;
      _event = 0;
      _lumi = 0;
      _ieta = 0;
      _E = 0;
      _shape = 0;
      _target2 = 0;      
      _HoverERatio = -1;      
      _nTT = -1;      
      _Run2IdLevel = -1;      
     
      if( _ele_L1Stage2_emul_rawEt>0 )
      {
	    _run = _nRun;
	    _event = _nEvent;
	    _lumi = _nLumi;

	    _ieta           = _ele_L1Stage2_emul_ieta;
	    _E              = _ele_L1Stage2_emul_rawEt;
	    _shape          = _ele_L1Stage2_emul_shape;
	    _target         = _ele_probePt/(0.5*_E);
	    _target2        = _ele_sclEt/(_E);
	    _HoverERatio    = _ele_L1Stage2_emul_hOverERatio;
	    _nTT            = _ele_L1Stage2_emul_nTT;
        
        if(_ele_LooseIdDecisions)
        {
            _Run2IdLevel = 1;
        }

	    tree->Fill();

      }
	
  }

  tree->Write();
  std::cout<<"\n";
  std::cout<<"Made File : "<<f_new->GetName()<<"\n";
  std::cout<<"\n";
  f_new->Close();

  return;
}


