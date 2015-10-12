// compile with 
// g++ -o inputFileCheck inputFileCheck.cpp `root-config --cflags --glibs` 

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 

int main (int argc, char** argv)
{
  
  //----------------------------------------------------------//
  //  Check input args                                        //
  //----------------------------------------------------------//
  if(argc == 0)
  {
    std::cout << "ERROR: YOU NEED TO PROVIDE AT LEAST AN INPUT FILE!" << std::endl;
  }
  else//if there is at least one file, go!
  {
    //play with strings for the output files...
    std::string firstFileName = argv[1];
    std::size_t found = firstFileName.find_first_of("_",firstFileName.find_first_of("_")+1);
    std::string fileRoot;
    fileRoot = firstFileName.substr(found+1,(firstFileName.length()-(found+1)) -5 );
    //std::cout << "fileRoot " << fileRoot << std::endl;
    //     if (Params.deepAnalysis)
    //       fileRoot += "_deep_";
    //     //std::string fileRoot;
    //     if (correctingForSaturation)
    //       fileRoot += "_saturation_corrected.root";
    //     else
    //       fileRoot += ".root";
    
    //a function for fitting...
    TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",0,20000); //quite incredibly, this means that [2] can be negative... so the fix will be to take always the module of [2]...
    
    
    //----------------------//
    //  TChain              //
    //----------------------//
    //tchain to merge the input root ttrees
    TChain *chain =  new TChain("adc");
    for (int i = 1 ; i < argc ; i++)
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      chain->Add(argv[i]);
    }
    
    
    
    
    
    
    
    
    std::stringstream snames[32],stypes[32];
    std::string names[32],types[32];
    //variables
    long long int t1_ExtendedTimeTag;
    long long int t1_DeltaTimeTag;
    int TriggerChannel[2];
    //Float_t t1_charge[32];
    int t1_charge[32];
    long long int counter = 0;
    float floodx[2],floody[2],firstonsecond[2];
    
    //branches
    TBranch *b_ExtendedTimeTag;
    TBranch *b_DeltaTimeTag;
    TBranch *b_charge[32];
    TBranch *b_TriggerChannel;
    TBranch *b_floodx;
    TBranch *b_floody;
    TBranch *b_firstonsecond;
    
    chain->SetBranchAddress("ExtendedTimeTag", &t1_ExtendedTimeTag, &b_ExtendedTimeTag);
    chain->SetBranchAddress("DeltaTimeTag", &t1_DeltaTimeTag, &b_DeltaTimeTag);  
    
    for(int i=0; i<32; i++)
    {
      snames[i].str(std::string());
      stypes[i].str(std::string());
      snames[i] << "ch" << i;
      names[i] = snames[i].str();
      chain->SetBranchAddress(names[i].c_str(), &t1_charge[i], &b_charge[i]);
    }
    chain->SetBranchAddress("TriggerChannel_0", &TriggerChannel[0], &b_TriggerChannel);
    chain->SetBranchAddress("FloodX_0", &floodx[0], &b_floodx);  
    chain->SetBranchAddress("FloodY_0", &floody[0], &b_floody);  
    chain->SetBranchAddress("FirstOnSecond_0", &firstonsecond[0], &b_firstonsecond);
    
//     Int_t nevent = chain->GetEntries();
    Int_t nevent = 100;
    
    for (Int_t i=0;i<nevent;i++) //loop on all the entries of tchain
    {
      chain->GetEvent(i);
      for(int i=0; i<32; i++)
      {
	if(!(i%2))
	  std::cout << t1_charge[i] << "\t";
      }
      std::cout << TriggerChannel[0] << "\t";
      std::cout << floodx[0] << "\t";
      std::cout << floody[0] << std::endl;
    }
    
    
  }
  
}