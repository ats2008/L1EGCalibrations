#include <string>
#include <TH2.h>
#include <TF1.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TGraphAsymmErrors.h>
#include <TObjString.h>
#include <iostream>
#include <fstream>
#include "IsolationAnalysis.h"

using namespace std;

IsolationAnalysis::IsolationAnalysis(const std::string& inputFileName) {
    doSuperCompression = 0 ;    
    doDynamicBinning=false;
    reportEvery=5000;
    maxEntries=-1;
    isoOffset=0;       
    ietaMinForOffset = 1e3;
    ietaMaxForOffset =-1e3;
    readParameters(inputFileName);
    
    puRewightMap.resize(80);
                                             //0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15      
    puRewightMap={22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 22.47891203305216, 3.449855825118133, 2.878747795414462, 2.419966323904733, 2.077106966259453, 1.8384299149631131, 1.6506667244542683, 1.500594288265926, 1.379246865965319, 1.2916800913443331, 1.2119789760643236, 1.157706220299312, 1.1084728332840172, 1.0757190603963658, 1.0517704023160444, 1.027920722601076, 1.0088294373441054, 1.0043776563715558, 1.0, 1.0053763204322175, 1.0229739953353658, 1.0473979484264853, 1.0814206615872415, 1.1254789743791802, 1.1833533566708785, 1.2618375787561362, 1.360953617457179, 1.4660614614743055, 1.615769154622847, 1.7814322242664256, 2.008163946815709, 2.2666990695736704, 2.5804867031788152, 3.002312351372302, 3.462296025817791, 4.063789301465357, 4.7940880292032055, 5.666410434437611, 6.870565243535779, 8.288838913272153, 9.933275374918496, 12.285752688172042, 15.174646390862607, 18.58147666287201, 22.961716237942124, 28.131847839468175, 36.93470179408437, 45.11648568608095, 57.88120567375886, 70.52932098765432, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777, 219.01760041822777} ;
    
    superCompressionToDefaultCompressionMap.resize(16);
                                             //0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15      
    superCompressionToDefaultCompressionMap={ 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 } ;
    if(doSuperCompression==101){
                                             // 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 
      superCompressionToDefaultCompressionMap={ 0,1,2,3,4,5,7,7,8,8, 8, 8, 9, 10, 11, 12};  // custom 101 scheme
      std::cout<<"\t Loading the super compression map for custom 101 scheme ! \n";
    }
    if(doSuperCompression==3){
                                            // 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 
      superCompressionToDefaultCompressionMap={ 0,0,1,1,2,2,3,3,4,4, 4, 5, 5, 6, 6, 7};  // 3 bit scheme
      std::cout<<"\t Loading the super compression map for 3 bit compression ! \n";
    }
    if(doSuperCompression==2){
                                            // 0,1,2,3,4,5,6,7,8,9,10,11,12,13,14,15 
      superCompressionToDefaultCompressionMap={0,0,0,0,1,1,1,1,2,2, 2, 2, 2, 3, 3, 3};  // 2 dit scheme
      std::cout<<"\t Loading the super compression map for 2 bit compression ! \n";
    }

    tmpFitMin = 10 ;
    tmpFitMax = 27 ;

    if (ntupleFileName_.size() == 0) {
        std::cout << " Inputfile list missing !!!" << ntupleFileName_ << std::endl;
        return;
    }

    accessTree(ntupleFileName_);

    assert(fChain);

    bookedHistograms_ = false;

    short last = -1;
    for (auto ieta : lutMapIEta_)
    {
        short val = ieta.second;
        if (val != last) {
            lutIEtaVec_.push_back(ieta.first);
        }
        last = val;
    }
    lutIEtaVec_.push_back(31);
    std::cout << "LUTIEtaBins: ";
    for (auto ieta : lutIEtaVec_) std::cout << ieta << " , " ;
    std::cout << " (Size " << lutIEtaVec_.size() << ")"<< std::endl;
    last = -1;
    for (auto iet : lutMapIEt_)
    {
        short val = iet.second;
        if (val != last)
        {
            lutIEtVec_.push_back(iet.first);
        }
        last = val;
    }
    lutIEtVec_.push_back(255);
    std::cout << "LUTIEtBins: ";
    for (auto iet : lutIEtVec_) std::cout << iet << " , " ;
    std::cout << " (Size " << lutIEtVec_.size()<< ")"<< std::endl;
    last = -1;
    for (auto intt : lutMapnTT_)
    {
        short val = intt.second;
        if (val != last) {
            lutnTTVec_.push_back(intt.first);
        }
        last = val;
    }
    lutnTTVec_.push_back(255);
    std::cout << "LUTnTTBins: ";
    for (auto intt : lutnTTVec_) std::cout << intt << " , " ;
    std::cout << " (Size " << lutnTTVec_.size() <<  ")" << std::endl;
}


IsolationAnalysis::~IsolationAnalysis() {
    /*  for (auto it : Histos_PerBin) {
        if (it.second) it.second->Delete();
        }
        Histos_PerBin.clear();*/
}
void IsolationAnalysis::accessTree(std::string & input_filelist) {
    std::ifstream myFile;
    myFile.open(input_filelist.c_str(), std::ios::in);
    if (!myFile) {
        std::cout << "Input File: " << input_filelist << " could not be opened!" << std::endl;
        fChain = 0;
    } else {
        fChain = new TChain("Ntuplizer/TagAndProbe");
        static constexpr int BUF_SIZE = 256;
        char buf[BUF_SIZE];
        while (myFile.getline(buf, BUF_SIZE, '\n')) {  // Pops off the newline character
            std::string line(buf);
            fChain->AddFile(line.c_str(),-1);
            std::cout << "Adding file " << line << " Entries " << fChain->GetEntries() <<  std::endl;
        }
    }
}

void IsolationAnalysis::analyse() {

    if (fChain == 0) return;

    Long64_t nentries = fChain->GetEntriesFast();
    std::cout<<"Total Available Entries : "<<nentries<<std::endl;
    Long64_t nbytes = 0, nb = 0;
    if(maxEntries > -1) nentries = maxEntries < nentries ? maxEntries :nentries ;
    std::cout<<"Processing a total of "<<nentries<<" Entries \n";

    auto t_start = std::chrono::high_resolution_clock::now();
    auto t_end = std::chrono::high_resolution_clock::now();

    for (Long64_t jentry=0; jentry<nentries; jentry++) {
        if (jentry < 0) break;
        nb = fChain->GetEntry(jentry);
        nbytes += nb;
        if(jentry%reportEvery == 0 )
        {
            t_end = std::chrono::high_resolution_clock::now();
            std::cout<<"  Loop 1/3 : Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
                     <<"  Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                     <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
                     <<std::endl;

        }

        if( isProbeLoose==0 ) continue;
        if( fabs(eleProbeEta) >= 2.5) continue;
        if( sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2)) < 0.6 ) continue;
        hprof_IEt->Fill(et, iso);
        hprof_IEta->Fill(eta, iso);
        hprof_nTT->Fill(ntt, iso);

        if(et < 0 )  continue;
        //if(Nvtx < 20 )  continue;
        //if(Nvtx > 69 )  continue;
        //auto wei = puRewightMap[int(Nvtx)];
        //std::cout<<" PU : "<<Nvtx<<" , weight : "<<wei<<"\n";
        auto wei=1.0;
        TString Name_Histo = getHistoName(eta,et, ntt);
        //std::cout<<__LINE__<<" | Name Histo : "<<Name_Histo<<" | eta "<<eta<<" | et "<<et<<" | ntt "<<ntt<<"\n";
        std::map<TString, TH1F*>::iterator iPos = Histos_PerBin.find(Name_Histo);
        if (iPos != Histos_PerBin.end()) iPos->second->Fill(iso,wei);
        if (eta > 31.5) {
            pt_large_eta->Fill(et);
            eta_large_eta->Fill(eta);
            nTT_large_eta->Fill(ntt);
            iso_large_eta->Fill(iso);
        }
    }

    std::cout << "  Total Number of Histograms " << Histos_PerBin.size() << std::endl;
    Int_t NumberOfHistosWithLowStats = 0;

    Int_t c=0;
    t_start = std::chrono::high_resolution_clock::now();
    t_end = std::chrono::high_resolution_clock::now();

    for (auto it : Histos_PerBin) {
        TString Name_Histo = it.first;
        TH1F* th = it.second;

        short ibin_eta;
        short ibin_et;
        short ibin_ntt;
        bool name_flag = getHistoBin(Name_Histo, ibin_eta, ibin_et, ibin_ntt);
        if (!name_flag) {
            std::cout << " ===> " <<  Name_Histo.Data() << ibin_eta << " " << ibin_et<<  " " << ibin_ntt << std::endl;
            continue;
        }
        if(th->GetEntries()<40) {
            //std::cout<<"Low stat Bin : "<<ibin_et<<" , "<<ibin_eta<<" , "<< ibin_ntt <<" , "<< th->GetEntries()<<"\n";
            NumberOfHistosWithLowStats++;
        }
        for(UInt_t iEff = 1 ; iEff < 101 ; ++iEff) {
            Float_t Efficiency = 0.01*iEff;

            for(UInt_t iIso = 0 ; iIso < 100 ; ++iIso) {
                if(th->Integral(1,iIso+1)/th->Integral(1,100+1)>=Efficiency) {
                    if (IsoCut_PerEfficiency_PerBin[iEff][Name_Histo] == -1) {
                        IsoCut_PerEfficiency_PerBin[iEff][Name_Histo] =  iIso;
                        IsoCut_PerBin[iEff]->SetBinContent(ibin_eta+1,ibin_et+1,ibin_ntt+1,iIso);
                    }
                }
            }


        }
        c++;
        if(c%200 == 0)
        {
            t_end = std::chrono::high_resolution_clock::now();
            std::cout<<"Loop 2/3 : Processing Hist in event loop : "<<c<<" / "<<Histos_PerBin.size()<<"  [ "<<100.0*c/Histos_PerBin.size()<<"  % ]  "
                     << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                     <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( Histos_PerBin.size() - c)/(1e-9 + c )* 0.001
                     <<std::endl;
        }
    }

    // Obtaining the fits and projections from the histogram
    int iMax=16*16*100;
    t_start = std::chrono::high_resolution_clock::now();
    t_end = std::chrono::high_resolution_clock::now();
    c=0;
    for(Int_t iEff = 100 ; iEff >= 0 ; --iEff)
    {   
        std::cout<<"\t Fitting loop ! processing for efficiency : "<<iEff<<"\n";
        for ( short iet =0 ; iet < lutIEtVec_.size()-1; iet++)
        {
            for (short ieta=0; ieta < lutIEtaVec_.size()-1 ; ieta++ )
            {
                TString projName = "pz_"+to_string(iEff)+"_eta"+to_string(ieta)+"_e"+to_string(iet);
                TH1D* projection = IsoCut_PerBin[iEff]->ProjectionZ(projName, ieta+1, ieta+1, iet+1, iet+1, "e");

                TString fitName = "fit_pz_"+to_string(iEff)+"_eta"+to_string(ieta)+"_e"+to_string(iet);
                TF1* projection_fit = new TF1(fitName,"[0]+[1]*x+[2]*x*x", tmpFitMin, tmpFitMax);
                projection_fit->SetParameter(0,0.0);
                projection_fit->SetParameter(1,0.5);
                projection_fit->SetParameter(2,0);
                //projection_fit->SetParLimits(1,0.0,20.0);
                projection_fit->SetParLimits(2,0.0,10.0);

                projection->Fit(projection_fit,"QR");
      //          std::cout<<" 1/2/3 : "<<projection_fit->GetParameter(0)<<" / "<<projection_fit->GetParameter(1) <<" / "<<projection_fit->GetParameter(2)<<"\n";              

                if (doDynamicBinning){
                    if ( projection_fit->GetParameter(1) < 0 )
                     {
                        for( short dX=1; dX <5  ; dX++)
                        {
                             projection = IsoCut_PerBin[iEff]->ProjectionZ(projName, 
                                                                             max(ieta-dX+1 , 0), 
                                                                             min(ieta+dX+1 , int(lutIEtVec_.size()-1)), 
                                                                             max(iet -dX+1 , 0 ),
                                                                             min(iet +dX+1 , int(lutIEtVec_.size()-1)), 
                                                                             "e");
                             projection->Fit(projection_fit,"QR");
                             if ( projection_fit->GetParameter(1) > 0 )  {
                                std::cout<<"\t CASE 2  [ dX = "<<dX<<" ] : Fit for eff : "<<iEff<<" , et : "<<iet<<" , eta : "<<ieta<<" to be taken from exteded region "
                                        <<"et [ "<<max(iet -dX , 0 )<<" , "<<min(iet +dX , int(lutIEtVec_.size()-1))<<" ] " 
                                        <<"eta [ "<<max(ieta-dX , 0)<<" , "<<min(ieta+dX , int(lutIEtVec_.size()-1))<<" ] ,"
                                        <<"\n";
                                 break;
                             }
                        } 
                     }

                     if ( projection_fit->GetParameter(1) < 0  and iEff < 100)
                     {
                         std::cout<<"\t CASE 3 : Fit for eff "<<iEff<<" , "<<iet<<" , "<<ieta<<" to be taken from eff "<<iEff-1<<"\n";
                         TString fitName_pre = "fit_pz_"+to_string(iEff+1)+"_eta"+to_string(ieta)+"_e"+to_string(iet);
                         
                         projection_fit = (TF1*) gDirectory->Get(fitName_pre)->Clone();    
                         projection_fit->SetName(fitName);
                         projection_fit->SetParameter(0 , projection_fit->GetParameter(0) -1 );
                     }
                      
                }

                if( projection_fit->GetParameter(1) < 0.01 )
                {
                    std::cout<<"\t Fit for iet : "<<iet<<" , ieta "<<ieta<<" is  "<<projection_fit->GetParameter(0)<<" + "<<projection_fit->GetParameter(1)<<" * x "
                             <<" | integral of the projection TH1D "<<projection->Integral()<<" \n";
                }
                projection->Write();
                projection_fit->Write();
 
                c++;
                if(c%100 == 0)
                {
                    t_end = std::chrono::high_resolution_clock::now();
                    std::cout<<"Loop 3/3 : Processing Hist in event loop : "<<c<<" / "<<iMax <<"  [ "<<100.0*c/iMax<<"  % ]  "
                             << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                             <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( iMax - c)/(1e-9 + c )* 0.001
                             <<std::endl;
                }
            }
        }

    }
    /*
          // for each bin in eta/e project icoCut on nTT and make a fit of it
          // this approach makes the isoCut strictly increasing and less dependent on statistic fluctuations

          for(UInt_t i = 0 ; i < NbinsIEta-1 ; ++i)
          {
            for(UInt_t j = 0 ; j < NbinsIEt-1 ; ++j)
              {
                TString projName = "pz_"+to_string(iEff)+"_eta"+to_string(i)+"_e"+to_string(j);
                TH1D* projection = IsoCut_PerBin[iEff]->ProjectionZ(projName, i+1, i+1, j+1, j+1, "e");
                projection->Write();

                TString fitName = "fit_pz_"+to_string(iEff)+"_eta"+to_string(i)+"_e"+to_string(j);
                TF1* projection_fit = new TF1(fitName,"[0]+[1]*x", FitMin, FitMax);
                projection->Fit(projection_fit);
                projection_fit->Write();
              }
          }



    */
    t_start = std::chrono::high_resolution_clock::now();
    if(false)
        for (Long64_t jentry=0; jentry<nentries; jentry++) {

            if(jentry%reportEvery == 0 )
            {
                t_end = std::chrono::high_resolution_clock::now();
                std::cout<<"Loop 3/3 : Processing Entry in event loop : "<<jentry<<" / "<<nentries<<"  [ "<<100.0*jentry/nentries<<"  % ]  "
                         << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                         <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( nentries - jentry)/(1e-9 + jentry)* 0.001
                         <<std::endl;
            }
            if (jentry < 0) break;

            nb = fChain->GetEntry(jentry);
            nbytes += nb;

            if( isProbeLoose==0 ) continue;
            if( fabs(eleProbeEta) >= 2.5) continue;
            if( sqrt(pow(eleProbeEta-eleTagEta,2)+pow(eleProbePhi-eleTagPhi,2)) < 0.6 ) continue;

            pt_all->Fill(et);
            if ( et < 0) continue;

            eta_all->Fill(eta);
            nTT_all->Fill(ntt);

            TString Name_Histo = getHistoName(eta,et, ntt);

            short ibin_eta;
            short ibin_et;
            short ibin_ntt;

            getHistoBin(Name_Histo, ibin_eta,  ibin_et, ibin_ntt);

            std::map<TString,TH1F*>::iterator iPos = Histos_PerBin.find(Name_Histo);
            if (iPos == Histos_PerBin.end()) continue;
            TH1F* th = iPos->second;
            if (!th->GetEntries()) continue;
            float eff_at_isolation = th->Integral(1,iso)/th->Integral(1,100+1);

            for(UInt_t iEff = 1 ; iEff < 101 ; ++iEff) 	{
                std::map<Int_t,std::map<TString,Int_t>>::iterator jPos = IsoCut_PerEfficiency_PerBin.find(iEff);

                if (jPos == IsoCut_PerEfficiency_PerBin.end()) continue;
                std::map<TString,Int_t>::iterator kPos = jPos->second.find(Name_Histo);
                if (kPos == jPos->second.end())  continue;
                if(iso<=kPos->second) {
                    eta_pass_efficiency[iEff]->Fill(eta);
                    pt_pass_efficiency[iEff]->Fill(et);
                    nTT_pass_efficiency[iEff]->Fill(ntt);
                }
                if(iso<=IsoCut_PerBin[iEff]->GetBinContent(ibin_eta+1, ibin_et+1, ibin_ntt+1)) pt_pass_efficiency_TH3[iEff]->Fill(et);
                //      if (eff_at_isolation <= iEff*0.01  ){
                //	eta_pass_efficiency[iEff]->Fill(eta);
                //	pt_pass_efficiency[iEff]->Fill(et);
                //	nTT_pass_efficiency[iEff]->Fill(ntt);
                //
                //	IsoCut_PerBin[iEff]->SetBinContent(ibin_eta+1, ibin_et+1, ibin_ntt+1, isolation);
                //	}
            }
        }
    outputFile_->cd();
    outputFile_->cd("Step1Histos");
    for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff) {
        TGraphAsymmErrors* temp_histo1 = new TGraphAsymmErrors(pt_pass_efficiency[iEff],pt_all,"cp");
        TGraphAsymmErrors* temp_histo_TH3 = new TGraphAsymmErrors(pt_pass_efficiency_TH3[iEff],pt_all,"cp");

        TGraphAsymmErrors* temp_histo2 = new TGraphAsymmErrors(eta_pass_efficiency[iEff],eta_all,"cp");

        TGraphAsymmErrors* temp_histo3 = new TGraphAsymmErrors(nTT_pass_efficiency[iEff],nTT_all,"cp");

        temp_histo1->Write();
        temp_histo_TH3->Write();
        temp_histo2->Write();
        temp_histo3->Write();
    }

    std::cout<<"NumberOfHistosWithLowStats/Tot = "<<NumberOfHistosWithLowStats<<"/"
             << Histos_PerBin.size()<<std::endl;
}

void IsolationAnalysis::bookHistograms(std::string option)
{

    outputFile_ = new TFile(outputFileName_.c_str(), "RECREATE");
    outputFile_->cd();

    if(option == "do_1" || option == "do_all") {
        TDirectory* td_step1 =   outputFile_->mkdir("Step1Histos");
        if (td_step1) td_step1->cd();
        for(UInt_t i = 0 ; i < nBinsIEta ; ++i) {
            for(UInt_t j = 0 ; j < nBinsIEt ; ++j) {
                for(UInt_t k = 0 ; k < nBinsnTT ; ++k) {
                    TString Name_Histo = "Hist_";
                    Name_Histo +=  i;
                    Name_Histo += "_";
                    Name_Histo += j;
                    Name_Histo += "_";
                    Name_Histo += k;

                    TH1F* temp_histo = new TH1F(Name_Histo.Data(),Name_Histo.Data(),100,0.,100.);
                    Histos_PerBin.insert({Name_Histo,temp_histo});
                }
            }
        }
        hprof_IEt  = new TProfile("hprof_IEt","Profile L1_Iso vs. L1_IEt",100,0.,200.,0,20);
        hprof_IEta  = new TProfile("hprof_IEta","Profile L1_Iso vs. L1_IEta",32,0.,32.,0,20);
        hprof_nTT  = new TProfile("hprof_nTT","Profile L1_Iso vs. L1_IEta",180,20.,200.,0,20);
        std::map<TString, Int_t> tempMap;
        for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff) {

            for (auto it : Histos_PerBin) {
                tempMap.insert({it.first, -1});
            }
            IsoCut_PerEfficiency_PerBin.insert({iEff, tempMap});

            TString NameEff = "Eff_";
            NameEff += iEff;
            TH3F* temp = new TH3F(NameEff,NameEff,nBinsIEta,0,nBinsIEta,nBinsIEt,0,nBinsIEt,nBinsnTT,0,nBinsnTT);
            IsoCut_PerBin.insert({iEff,temp});

            TString nameHisto1 = "pt_pass_efficiency_";
            nameHisto1 += iEff;

            TString nameHisto_TH3 = "pt_pass_efficiency_TH3_";
            nameHisto_TH3 += iEff;

            TH1F* temp_histo = new TH1F(nameHisto1,nameHisto1,100,0,200);
            TH1F* temp_histo_TH3 = new TH1F(nameHisto_TH3,nameHisto_TH3,100,0,200);

            pt_pass_efficiency.insert({iEff,temp_histo});
            pt_pass_efficiency_TH3.insert({iEff,temp_histo_TH3});

            TString nameHisto2 = "eta_pass_efficiency_";
            nameHisto2 += iEff;

            TH1F* temp_histo2 = new TH1F(nameHisto2,nameHisto2,100,0,100);
            eta_pass_efficiency.insert({iEff,temp_histo2});

            TString nameHisto3 = "nTT_pass_efficiency_";
            nameHisto3 += iEff;

            TH1F* temp_histo3 = new TH1F(nameHisto3,nameHisto3,180,20,200);
            nTT_pass_efficiency.insert({iEff,temp_histo3});
        }

        pt_large_eta = new  TH1F("Pt_large_eta",  "Pt_large_eta",  100, 0., 200.);
        eta_large_eta = new TH1F("Eta_large_eta", "Eta_large_eta", 100, 0., 100.);
        nTT_large_eta = new TH1F("NTT_large_eta", "NTT_large_eta", 180, 20., 200.);
        iso_large_eta = new TH1F("Isolation_large_eta", "Isolation_large_eta", 100, 0., 100.);

        pt_all = new  TH1F("Pt_all",  "Pt_all",  100, 0., 200.);
        eta_all = new TH1F("Eta_all", "Eta_all", 100, 0., 100.);
        nTT_all = new TH1F("NTT_all", "NTT_all", 180, 20., 200.);

    }
    if(option == "do_2" || option == "do_all") {
        TDirectory* td_step2 =   outputFile_->mkdir("Step2Histos");
        if (td_step2) td_step2->cd();
        int i=0;
        for (auto it : lutProgOptVec_) {
            std::string tmp_string = it;
            std::vector<std::string> options;
            tokenize(tmp_string,options,":");

            TString hname = "LUT_Progression_";
            hname += options[0];
            TH3F* th3 = new TH3F(hname,hname,
                                 lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
                                 lutIEtVec_.size()-1, 0, lutIEtVec_.size()-1,
                                 lutnTTVec_.size()-1, 0, lutnTTVec_.size()-1);
            lutProgHistoMap_.insert({it, th3});

            // Booking the Extrapolated LUT option
            hname = "LUT_Progression_v2_";
            hname += options[0];
            th3 = new TH3F(hname,hname,
                           lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
                           lutIEtVec_.size()-1, 0, lutIEtVec_.size()-1,
                           lutnTTVec_.size()-1, 0, lutnTTVec_.size()-1);
            lutProgHistoMap_v2_.insert({it, th3});

            // Booking the Extrapolated LUT option with 0 supression below the min threashold
            hname = "LUT_Progression_v3_";
            hname += options[0];
            th3 = new TH3F(hname,hname,
                           lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
                           lutIEtVec_.size()-1, 0, lutIEtVec_.size()-1,
                           lutnTTVec_.size()-1, 0, lutnTTVec_.size()-1);
            lutProgHistoMap_v3_.insert({it, th3});


        }
        for(UInt_t iEff = 0 ; iEff <= 100 ; ++iEff)
        {
            TString NameHisto = "LUT_WP";
            NameHisto += iEff;
            TH3F* LUT_temp = new TH3F(NameHisto.Data(),NameHisto.Data(),
                                      lutIEtaVec_.size()-1, 0, lutIEtaVec_.size()-1,
                                      lutIEtVec_.size()-1 , 0, lutIEtVec_.size()-1 ,
                                      lutnTTVec_.size()-1 , 0, lutnTTVec_.size()-1  );
            LUT_WP.push_back(LUT_temp);
        }
    }
    bookedHistograms_ = true;
}
void IsolationAnalysis::fillLUTProgression(std::string option) {
    Double_t AvEt;
    TFile*WPFile;
    if(option == "do_2") {
        WPFile = TFile::Open(outputWPFileName_.c_str(), "READ");
        if (!WPFile) return;
        std::cout<<"Filling LUT only"<<std::endl;
    }
    std::cout<<"Filling LUT"<<std::endl;
    Int_t maxEntriesForLUT=(lutIEtaVec_.size()-1)*(lutIEtVec_.size()-1)*(lutProgHistoMap_.size());
    std::cout<<"Number of pT , eta , option triplets :  "<<(lutIEtaVec_.size()-1)<<"*"<<(lutIEtVec_.size()-1)<<"*"<< lutProgHistoMap_.size()  <<" => "<<maxEntriesForLUT<<"\n";
    Int_t count=0;


    auto  t_start = std::chrono::high_resolution_clock::now();
    auto  t_end = std::chrono::high_resolution_clock::now();
    for (auto it : lutProgHistoMap_)
    {

        std::string tmp_string = it.first;
        std::vector<std::string> options;
        tokenize(tmp_string,options,":");
        Double_t minPt =  std::stod(options[1]);
        Double_t effLowMinPt = std::stod(options[2]);
        Double_t reach100pc= std::stod(options[3]);
        
        std::cout<<"Is using super Compressed vars : "<<doSuperCompression<<"\n";

        for(Int_t j = 0 ; j < lutIEtVec_.size()-1 ; j++)
        {


            // Obtaining the Efficiency valut from the LUT option
            Double_t Efficiency_Progression       = findEfficiencyProgression((lutIEtVec_[j]+lutIEtVec_[j+1])/2.0, minPt,effLowMinPt,reach100pc);
            Double_t Efficiency_Progression_forV3 = findEfficiencyProgressionForV3((lutIEtVec_[j]+lutIEtVec_[j+1])/2.0, minPt,effLowMinPt,reach100pc);

            if(Efficiency_Progression >= 0.9999) Efficiency_Progression = 1.0001;
            Int_t Int_Efficiency_Progression = int(Efficiency_Progression*100);

            if(Efficiency_Progression_forV3 >= 0.9999) Efficiency_Progression = 1.0001;
            Int_t Int_Efficiency_Progression_forV3 = int(Efficiency_Progression_forV3*100);
            
            //std::cout<<"\t for v2 : "<<Int_Efficiency_Progression<<" |  for v3 : "<< Int_Efficiency_Progression_forV3<<" @ Et : "<< 0.5*( lutIEtVec_[j]+lutIEtVec_[j+1]) <<"\n";
            // Loading corresponding Efficiency WP
            TH3F* eff_histo;
            if(option == "do_2") {
                TString WPName = "Eff_";
                WPName +=std::to_string(Int_Efficiency_Progression);
                eff_histo = dynamic_cast<TH3F*>(WPFile->Get("Step1Histos/"+WPName));
                //std::cout<<WPName<<std::endl;
            }
            else
                eff_histo = IsoCut_PerBin[Int_Efficiency_Progression];

            // Loading corresponding Efficiency WP for the Zero supressed 
            TH3F* eff_histo_forV3;
            if(option == "do_2") {
                TString WPName = "Eff_";
                WPName +=std::to_string(Int_Efficiency_Progression_forV3);
                eff_histo_forV3 = dynamic_cast<TH3F*>(WPFile->Get("Step1Histos/"+WPName));
                //std::cout<<WPName<<std::endl;
            }
            else
                eff_histo_forV3 = IsoCut_PerBin[Int_Efficiency_Progression_forV3];

            Int_t i = 0 ;
            //for(Int_t ii = 0 ; ii < 16 ; ii++)
            for(Int_t ii = 0 ; ii < lutIEtaVec_.size()-1; ii++)
            {
                if (doSuperCompression > 0)     i = superCompressionToDefaultCompressionMap[ii];
                else                        i = ii ;

                count++;
                if(count%500==0){
                    t_end = std::chrono::high_resolution_clock::now();
                    std::cout<<"LUT Filling : Processing  loop id : "<<count<<" / "<<maxEntriesForLUT<<"  [ "<<100.0*count/maxEntriesForLUT<<"  % ]  "
                         << " Elapsed time : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()/1000.0
                         <<"  Estimated time left : "<< std::chrono::duration<double, std::milli>(t_end-t_start).count()*( maxEntriesForLUT - count)/(1e-9 + count)* 0.001
                         <<std::endl;
                }
                

                // Loading the fit for correponding [ Eff , Et , Eta ]  V2 Scheme
                TString fitName = "fit_pz_"+to_string(Int_Efficiency_Progression)+"_eta"+to_string(i)+"_e"+to_string(j);
                TF1* currentFit = (TF1*) (WPFile->Get("Step1Histos/"+fitName));

                if ( not currentFit)
                {
                    std::cout<<" Fit not found for "<< fitName<<" ! \n \t Exiting !!  \n";
                    exit(1);
                }
                for(Int_t k = 0 ; k < lutnTTVec_.size()-1; k++)
                {
                    Int_t IsoCut_Progression = eff_histo->GetBinContent(i+1,j+1, k+1) ;
                    if( ( i <= ietaMaxForOffset ) and ( i>= ietaMinForOffset) )
                        {
                         std::cout<<"DEF Offseting for ieta : "<<i<<" \n";
                         IsoCut_Progression+=isoOffset;
                        }
                    if(Int_Efficiency_Progression==100) IsoCut_Progression = 1000;
                    else if(Int_Efficiency_Progression== 0) IsoCut_Progression = 0;
                    it.second->SetBinContent(ii+1,j+1,k+1,IsoCut_Progression);

                    // Filling from the Fit extrapolation
                    if(Int_Efficiency_Progression==100) IsoCut_Progression = 1000 ;
                    else if(Int_Efficiency_Progression== 0 ) IsoCut_Progression = 0 ;
                    else{
                            //IsoCut_Progression = max(Int_t(currentFit->GetParameter(0) + currentFit->GetParameter(1) *k  ) , 0);
                            IsoCut_Progression = Int_t(currentFit->GetParameter(0) + currentFit->GetParameter(1) *k  + currentFit->GetParameter(2) *k*k);
                            if (IsoCut_Progression==0) IsoCut_Progression=1;
                            if( ( i <= ietaMaxForOffset ) and ( i>= ietaMinForOffset) )
                            {
                                std::cout<<"V2 Offseting for ieta : "<<i<<" \n";
                                IsoCut_Progression+=isoOffset;
                            }
                           // std::cout<<"For [ "<<i<<","<<j<<" ]  nTT  = "<<k<<" : "<<
                    }
                    lutProgHistoMap_v2_[it.first]->SetBinContent(ii+1,j+1,k+1,IsoCut_Progression);


               }

                // Loading the fit for correponding [ Eff , Et , Eta ]  V3 Scheme
                fitName = "fit_pz_"+to_string(Int_Efficiency_Progression_forV3)+"_eta"+to_string(i)+"_e"+to_string(j);
                currentFit = (TF1*) (WPFile->Get("Step1Histos/"+fitName));

                if ( not currentFit)
                {
                    std::cout<<" Fit not found for "<< fitName<<" for v3 ! \n \t Exiting !!  \n";
                    exit(1);
                }

                for(Int_t k = 0 ; k < lutnTTVec_.size()-1; k++)
                {
                    // Filling from the Fit extrapolation in V3 Scheme
                    Int_t IsoCut_Progression=1000;
                    if(Int_Efficiency_Progression==100) IsoCut_Progression      = 1000 ;
                    else if(Int_Efficiency_Progression== 0 ) IsoCut_Progression = 0 ;
                    else{
                            IsoCut_Progression = max(Int_t(currentFit->GetParameter(0) + currentFit->GetParameter(1) *k + currentFit->GetParameter(2) *k*k) , 0) ;
                            if( ( i <= ietaMaxForOffset ) and ( i>= ietaMinForOffset) )
                            {
                                std::cout<<"V3 Offseting for ieta : "<<i<<" \n";
                                IsoCut_Progression+=isoOffset;
                            }
                    }

                    lutProgHistoMap_v3_[it.first]->SetBinContent(ii+1,j+1,k+1,IsoCut_Progression);
                }

            }
        }
    }
}

Double_t IsolationAnalysis::findEfficiencyProgression(Double_t IEt, Double_t MinPt,
        Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at) {
    
    Double_t Efficiency = 0;
    Double_t Pt         = IEt/2.;
    
    if(Pt>=Reaching_100pc_at) Efficiency = 1.;
    else if(Pt<MinPt) Efficiency = Efficiency_low_MinPt;
    else
    {
        Double_t  Slope = (1.-Efficiency_low_MinPt)/(Reaching_100pc_at-MinPt);
        Efficiency = Slope*Pt + (1. - Slope*Reaching_100pc_at);
    }

    if(Efficiency<0) Efficiency = 0.;
    if(Efficiency>=1) Efficiency = 1.;

    return Efficiency ;
}

Double_t IsolationAnalysis::findEfficiencyProgressionForV3(Double_t IEt, Double_t MinPt,
        Double_t Efficiency_low_MinPt, Double_t Reaching_100pc_at) {
    
    Double_t Efficiency = 0;
    Double_t Pt         = IEt/2.;
    
    if(Pt>=Reaching_100pc_at) Efficiency = 1.;
    else if(Pt<MinPt) Efficiency = 0.0;
    else
    {
        Double_t  Slope = (1.-Efficiency_low_MinPt)/(Reaching_100pc_at-MinPt);
        Efficiency = Slope*Pt + (1. - Slope*Reaching_100pc_at);
    }

    if(Efficiency<0) Efficiency = 0.;
    if(Efficiency>=1) Efficiency = 1.;
    
  //   std::cout<<"\t v3 | "<<IEt<<" , "<<MinPt<<" , "<<Efficiency_low_MinPt<<" , "<<Reaching_100pc_at<<" Eff . : "<<Efficiency <<"\n";

    return Efficiency ;
}




void IsolationAnalysis::readLUTTable(std::string& file_name, unsigned int& nbin,
                                     std::map<short, short>& lut_map) {

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
            // enable '#' and '//' style comments
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
void IsolationAnalysis::saveHistograms() {
    if (outputFile_ && bookedHistograms_) {
        /*    for(UInt_t iEff = 0 ; iEff < 101 ; ++iEff) {
          pt_pass_efficiency[iEff]->Delete();
          pt_pass_efficiency_TH3[iEff]->Delete();
          eta_pass_efficiency[iEff]->Delete();
          nTT_pass_efficiency[iEff]->Delete();

          }*/
        /*    for (auto it : Histos_PerBin) {
          if (it.second) it.second->Delete();
        }*/

        outputFile_->Purge();
        outputFile_->Write();
        outputFile_->Close();
    }
}
void IsolationAnalysis::readParameters(const std::string jfile) {
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
            if(key=="NtupleFileName")        ntupleFileName_= value;
            else if (key=="OutputWPStepFileName") outputWPFileName_ = value.c_str();
            else if (key=="DoSuperCompression")  doSuperCompression = atoi(value.c_str());
            else if (key=="ietaMinForOffset")  ietaMinForOffset = atoi(value.c_str());
            else if (key=="ietaMaxForOffset")  ietaMaxForOffset = atoi(value.c_str());
            else if (key=="IsoOffset")  isoOffset = atoi(value.c_str());
            else if (key=="DoDynamicBinning")  doDynamicBinning = atoi(value.c_str()) > 0;
            else if (key=="OutputFileName")  outputFileName_ = value.c_str();
            else if (key=="EtLUTFileName")	readLUTTable(value,nBinsIEt, lutMapIEt_);
            else if (key=="EtaLUTFileName")	readLUTTable(value,nBinsIEta, lutMapIEta_);
            else if (key=="NTTLUTFileName")	readLUTTable(value,nBinsnTT, lutMapnTT_);
            else if (key=="LUTProgressionOptions")
            {
                std::string tmp_string = value;
                std::vector<std::string> tmp_vec;
                tokenize(tmp_string,tmp_vec,",");
                for (auto it : tmp_vec) {
                    lutProgOptVec_.push_back(it);
                }
            }
            else if (key=="ReportEvery")	    {
                reportEvery= atoi(value.c_str());
            }
            else if (key=="MaxEntries")	    {
                maxEntries= atoi(value.c_str());
            }
            else
                std::cout << " unknown option " << " key " << key << std::endl;
        }
    }
    jobcardFile.close();

}
void IsolationAnalysis::readTree()
{
    fChain->SetMakeClass(1);

    fChain->SetBranchAddress("l1tEmuNTT",&ntt);
    fChain->SetBranchAddress("l1tEmuRawEt",&et);
    fChain->SetBranchAddress("l1tEmuTowerIEta",&eta);
    fChain->SetBranchAddress("l1tEmuIsoEt",&iso);

    fChain->SetBranchAddress("eleProbeEta",&eleProbeEta			);
    fChain->SetBranchAddress("eleProbePhi",&eleProbePhi			);
    fChain->SetBranchAddress("eleTagEta",&eleTagEta		    	);
    fChain->SetBranchAddress("eleTagPhi",&eleTagPhi			    );
    fChain->SetBranchAddress("Nvtx",&Nvtx			    );
    fChain->SetBranchAddress("isProbeLoose",&isProbeLoose			);
}

bool IsolationAnalysis::getHistoBin(TString str, short& eta_bin, short& et_bin, short& ntt_bin) {
    bool rval = false;
    if (str.Length() == 0) rval;

    str.Remove(0,5);
    TObjArray* toa = str.Tokenize("_");
    if (toa->GetEntries() >= 3) {
        eta_bin =  (dynamic_cast<TObjString*>(toa->At(0)))->GetString().Atoi();
        et_bin  =  (dynamic_cast<TObjString*>(toa->At(1)))->GetString().Atoi();
        ntt_bin =  (dynamic_cast<TObjString*>(toa->At(2)))->GetString().Atoi();
        rval = true;
    } else rval = false;
    return rval;
}
TString IsolationAnalysis::getHistoName(short eta, short et, short ntt) {

    TString Name_Histo = "Hist_";

    unsigned int jEta, jEt, jnTT;
    std::map<short, short >::iterator jEtaPos = lutMapIEta_.find(abs(eta));
    if (jEtaPos != lutMapIEta_.end()) jEta = jEtaPos->second;
    else                          jEta = nBinsIEta-1;

    std::map<short, short >::iterator jEtPos = lutMapIEt_.find(et);
    if (jEtPos != lutMapIEt_.end()) jEt = jEtPos->second;
    else                          jEt = nBinsIEt-1;

    std::map<short, short >::iterator jnTTPos = lutMapnTT_.find(ntt);
    if (jnTTPos != lutMapnTT_.end()) jnTT = jnTTPos->second;
    else                          jnTT = nBinsnTT-1;

    Name_Histo += jEta;
    Name_Histo += "_";
    Name_Histo += jEt;
    Name_Histo += "_";
    Name_Histo += jnTT;
    return Name_Histo;
}
void IsolationAnalysis::tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {

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
    std::string to_do = argv[2];

    IsolationAnalysis treeReader(data_file);
    treeReader.readTree();

    treeReader.bookHistograms(to_do);

    if(to_do == "do_1" || to_do == "do_all") {
        std::cout << " Calling Analyse" << std::endl;
        treeReader.analyse();
    }
    if(to_do =="do_2" || to_do == "do_all") {
        std::cout << " Calling LUT Options" << std::endl;
        treeReader.fillLUTProgression(to_do);
    }
    treeReader.saveHistograms();
    return 0;
}
