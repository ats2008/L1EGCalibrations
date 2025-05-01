#include <fstream>
#include <iostream>
#include <sstream>
#include <TString.h>

using namespace std;

TString header_=R"""(# Calibration vs |ieta|,shape,E. Derived from Run 283478 data, with semi-parametric regression 
# The LUT output is (ET_off/ET_L1) between 0 and 2, encoded on 9 bits 
# Index is compressedShape+compressedE<<4+compressedIeta<<8. 
# Compression version is v4 
#anything after # is ignored with the exception of the header 
#the header is first valid line starting with #<header> versionStr(unused but may be in future) nrBitsAddress nrBitsData </header> 
#<header> V8 12 10 </header> ;

)""";
const int HeaderLines=8;

vector<int> convert_LUT_calibr(TString LUTfile){
  std::ifstream infile(LUTfile);
  vector<int> converted_LUT;
  int i=0,x=0;
  string line;
  while(getline(infile,line)){
    i++;
    
    if (i<=HeaderLines) std::cout<<"skipping line ["<<i<<"] : "<<line<<"\n";
    if(i>HeaderLines){
      std::istringstream iss(line);
      int a,b;
      iss>>a>>b;
      b &= (0x1ff);
      converted_LUT.push_back(b);
       x++;
    }
  }
  std::cout<<" Read "<<x<<" items from LUT : "<<LUTfile<<"\n";

  return converted_LUT;

}





vector<int> convert_shapeLUT(TString LUTfile){

  std::ifstream infile(LUTfile);
  vector<int> converted_LUT;

  int i=0;
  string line;
  int x=0;
  while(getline(infile,line)){
    i++;
    if(i>HeaderLines){
      std::istringstream iss(line);
      int a,b;
      iss>>a>>b;
      converted_LUT.push_back(b);
      x++;
    }
  }

  std::cout<<" Read "<<x<<" items from LUT : "<<LUTfile<<"\n";
  return converted_LUT;

}


void write_calibrLUT(TString LUTfile_calib, TString LUTfile_shapeID, TString LUTfile_out){

  std::ofstream out(LUTfile_out);
  out<<header_;
  vector<int> converted_LUT=convert_LUT_calibr(LUTfile_calib);
  vector<int> converted_shapeLUT=convert_shapeLUT(LUTfile_shapeID);


  for(unsigned int ieta=0; ieta<16; ieta++){
    for(unsigned int iET=0; iET<16; iET++){
      for(unsigned int shape=0; shape<16; shape++){  

	int i=(ieta*16*16)+iET*16+shape;
	
	int LUT_entry=converted_LUT[i];

	//cout<<"converted_shapeLUT[i]="<<converted_shapeLUT[i]<<endl;
	if(converted_shapeLUT[i]==1)
	  LUT_entry+=(1<<9);
	out<<i<<" "<<LUT_entry<<" # compressedIeta="<<ieta<<",compressedE="<<iET<<",compressedShape="<<shape<<" shapeID="<<converted_shapeLUT[i]<<",calibr="<<converted_LUT[i]<<'\n';


      }
    }
  }


}

void makeCorrectedLUT( TString LUTfile_calib  ="calibration_april2024_newcaloL1_oldHcal.txt" ,
                       TString LUTfile_out    ="correctedLUT.txt"
                       ){

  TString LUTfile_shapeID="data/shapeIdentification_adapt0.99_compressedieta_compressedE_compressedshape_v17.05.19.txt" ;
  
  write_calibrLUT(LUTfile_calib,LUTfile_shapeID,LUTfile_out);

}
