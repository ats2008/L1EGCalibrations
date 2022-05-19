#include "TROOT.h"
#include "TGraphAsymmErrors.h"
#include "TLegend.h"
#include "TFile.h"
#include "TH1F.h"
#include "TMath.h"
#include "TCanvas.h"
#include "TMultiGraph.h"
#include "TPaveText.h"
#include "TAxis.h"
#include "iostream"
#include "string"
#include "vector"

using namespace std;
//Double_t getIntegral(TGraphAsymmErrors* aGraph, Double_t midX,Double_t lX,Double_t rX,TString);
//void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) ;


void tokenize(const std::string& str, std::vector<std::string>& tokens, const std::string& delimiters) {
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

// Assumes x vals in graph are sorted in acending order
Double_t getIntegral(TGraphAsymmErrors* aGraph, Double_t midX,Double_t lX,Double_t rX,TString prefix,TGraphAsymmErrors* baseline=nullptr,Bool_t saveGraph=true)
{
    Int_t iX_beg=0;
    Int_t iX_end=0;
    Int_t iX_mid=0;
    Int_t nX= aGraph->GetN();
    Double_t a1=0.0;

    Double_t xArea[10000];
    Double_t yArea[10000];
    
    for(Int_t i=1; i< nX; i++)
    {
        if( (aGraph->GetPointX(i) >= lX ) and ( (aGraph->GetPointX(i-1)) <= lX) ){
            iX_beg=i;
         }
        if( (aGraph->GetPointX(i) >= midX ) and ( (aGraph->GetPointX(i-1)) <= midX) ){
            iX_mid=i;
         }
        if( (aGraph->GetPointX(i) > rX ) and ( (aGraph->GetPointX(i-1)) <= rX) ){
            iX_end=i;
         }
    }

    double x=0;
    Double_t y_mid = aGraph->GetPointY(iX_mid-1) +  
                        ( aGraph->GetPointY(iX_mid) - aGraph->GetPointY(iX_mid-1) ) * ( midX - aGraph->GetPointX(iX_mid-1))/(aGraph->GetPointX(iX_mid) - aGraph->GetPointX(iX_mid-1)) ;
    Double_t y_left = aGraph->GetPointY(iX_beg-1) +  
                        ( aGraph->GetPointY(iX_beg) - aGraph->GetPointY(iX_beg-1) ) * ( lX - aGraph->GetPointX(iX_beg-1))/(aGraph->GetPointX(iX_beg) - aGraph->GetPointX(iX_beg-1));
    Double_t y_right  = aGraph->GetPointY(iX_end-1) +  
                        ( aGraph->GetPointY(iX_end) - aGraph->GetPointY(iX_end-1) ) * ( rX - aGraph->GetPointX(iX_end-1))/(aGraph->GetPointX(iX_end) - aGraph->GetPointX(iX_end-1)) ;

    // LEFT AREA
    a1=0.5*(aGraph->GetPointX(iX_beg) - lX)*(aGraph->GetPointY(iX_beg) + y_left);
    Int_t idx=0;
    xArea[idx]=lX;  yArea[idx]=y_left;    idx++;
    for(Int_t i=iX_beg; i<iX_mid ; i++)
    {
        a1+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        x+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        xArea[idx]=aGraph->GetPointX(i);
        yArea[idx]=aGraph->GetPointY(i);
        idx++;
    }
   
    a1 -= 0.5*(aGraph->GetPointX(iX_mid) - midX)*(aGraph->GetPointY(iX_mid) + y_mid);
    a1=(midX-lX)*(y_mid)-a1;
    xArea[idx]=midX;  yArea[idx]=y_mid;    idx++;

    // RIGHT AREA
    Double_t a2=0.5*(aGraph->GetPointX(iX_mid) - midX)*( y_mid + aGraph->GetPointY(iX_mid) );
    x=0;
    for(Int_t i=iX_mid; i<iX_end ; i++)
    {
        a2+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        x+=0.5*(aGraph->GetPointX(i+1) - aGraph->GetPointX(i))*(aGraph->GetPointY(i+1) + aGraph->GetPointY(i));
        xArea[idx]=aGraph->GetPointX(i);
        yArea[idx]=aGraph->GetPointY(i);
        idx++;
    }
   
    a2 -= 0.5*(aGraph->GetPointX(iX_end) - rX)*(aGraph->GetPointY(iX_end) + y_right);
    a2 -= y_mid*( rX - midX) ;

    xArea[idx]=rX;    yArea[idx]=y_right;  idx++;
    xArea[idx]=rX;    yArea[idx]=y_mid;  idx++;
    xArea[idx]=lX;    yArea[idx]=y_mid;   idx++;
    
    Double_t A= a1+a2;
    
   if(saveGraph)
   {

    TCanvas * acanvas=new TCanvas("aCanvas","aCanvas");
    TMultiGraph *multigraph = new TMultiGraph();
    multigraph->SetTitle("Efficiency vs pT");
    multigraph->GetXaxis()->SetTitle("Oflline p_{T} [ GeV ]");
    multigraph->GetYaxis()->SetTitle("#epsilon");
    
    auto excl1 = new TGraph(idx,xArea,yArea);
    excl1->SetLineColor(41);
    excl1->SetLineWidth(3);
    excl1->SetFillColor(kSpring);
    excl1->SetFillStyle(1001);
    
    aGraph->SetLineWidth(2);
    aGraph->SetLineColor(kBlue);

    multigraph->Add(excl1,"F");
    multigraph->Add(aGraph,"L");
    multigraph->SetMinimum(0.);
    multigraph->SetMaximum(1.0);
    aGraph->Draw();
    if(baseline)
    {
     baseline->SetLineWidth(2);
     baseline->SetLineColor(kRed);
     multigraph->Add(baseline,"L");
    }
    multigraph->Draw("A* same");
    multigraph->GetXaxis()->SetLimits(0.0,80.0);
    
   
      TPaveText *pt = new TPaveText(45.0,0.3,75.0,0.6);
      
      std::string hname(aGraph->GetName());
      std::vector<string> tockens;
      tockens.clear();
      tokenize(hname,tockens,"_");
      std::string tag;
      if(tockens.size()>6)
      {
        std::replace( tockens[4].begin(), tockens[4].end(), 'p', '.');
        std::replace( tockens[5].begin(), tockens[5].end(), 'p', '.');
        std::replace( tockens[6].begin(), tockens[6].end(), 'p', '.');
        tag="[ "+tockens[4]+" | "+tockens[5]+" | "+tockens[6]+ " ]";
      }
      else
      {
        tag="Baseline";
      }
      pt->AddText(tockens[3].c_str());
      pt->AddText(tag.c_str());
      TString midXstr("");
      midXstr+=midX;
      pt->AddText("E_{T}^{L1EG} > "+midXstr);
      
      std::string st=std::to_string(A);
      pt->AddText(TString("A = ")+ st.c_str());
      pt->SetTextColor(1);
      pt->SetTextSize(0.04);
      pt->SetShadowColor(0);
      pt->Draw();

      //TLegend *leg = new TLegend(40.0, 0.15, 78.0,0.25);
      TLegend *leg = new TLegend(0.55, 0.23, 0.85,0.33);
      leg->SetFillColor(0);
      leg->AddEntry(aGraph,tockens[3].c_str(), "l");
      if (baseline)
       leg->AddEntry(baseline,"Baseline", "l");
      leg->Draw();

       acanvas->SaveAs(prefix+TString(aGraph->GetName())+".png","q");
    }
    return A;
}

bool isGoodTurnON(TGraphAsymmErrors* baseline, TGraphAsymmErrors* newOpt,Int_t point1=8, Int_t point2=14,bool verbose=false) {
// From L1 Tau  isolation LUT optimization workflow : Jona Motta , Olivier
// https://github.com/jonamotta/TauObjectsOptimization/blob/master/PlotGoodGridSearch/GoodTurnOns_gridSearch.C

    Double_t y;
    Double_t th;
    Double_t _;
    bool good = true;

    for (int i = point1; i <= point2; ++i)
    {
        th=baseline->GetPointY(i);
        y=newOpt->GetPointY(i);
        if(verbose)
            std::cout<<"th = "<<th<<" , y = "<<y<<"\n";
        if ( (th - y) > 0.04 ) 
        {
            good = false;
            break;
        }
    }
    if(verbose)
    std::cout<<"\t Returning : "<<good<<"\n";
    
    return good;
}


