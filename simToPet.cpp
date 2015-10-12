//simple program to translate sim output to "adc" format

// compile with 
// g++ -o simToPet simToPet.cpp `root-config --cflags --glibs`
// syntax
// simToPet `ls out*`

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  
  //TFile* f1 = new TFile(argv[1]);
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  
  
  
  long int Seed;
  int Run;
  int Event;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  std::vector<float> CryEnergyDeposited[64]; 
  std::vector<float> *pCryEnergyDeposited[64];
  std::vector<float> PosXEnDep[64]; 
  std::vector<float> *pPosXEnDep[64];
  std::vector<float> PosYEnDep[64]; 
  std::vector<float> *pPosYEnDep[64];
  std::vector<float> PosZEnDep[64]; 
  std::vector<float> *pPosZEnDep[64];
  short RunDetectorHit[16];
  
  
  std::vector<float> *pEdep[64];
  std::vector<float> *px[64];
  std::vector<float> *py[64];
  std::vector<float> *pz[64];
  
  for (int i = 0 ; i < 64 ; i++)
  {
    pEdep[i] = 0; 
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
  }
  Short_t  detector[16];
  
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  
  for (int i = 0 ; i < 64 ; i++)
  {
    std::stringstream snames;
    snames << "cry" << i;
    tree->SetBranchAddress(snames.str().c_str(),&pEdep[i]);
    snames.str("");
    snames<< "cry" << i << "PosXEnDep";    
    tree->SetBranchAddress(snames.str().c_str(),&px[i]);
    snames.str("");
    snames<< "cry" << i << "PosYEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&py[i]);
    snames.str("");
    snames<< "cry" << i << "PosZEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&pz[i]);
  }
  for (int i = 0 ; i < 16 ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  
  
  
  
  
  //output ttree
  
  long long int DeltaTimeTag,ExtendedTimeTag;
  Short_t charge[32];
  Float_t RealX,RealY,RealZ;
  
  TTree* t1 = new TTree("adc","adc");
  
  t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
  //branches of the 32 channels data
  for (int i = 0 ; i < 32 ; i++)
  {
    //empty the stringstreams
    std::stringstream snames,stypes;
    charge[i] = 0;
    snames << "ch" << i;
    stypes << "ch" << i << "/S";  
    t1->Branch(snames.str().c_str(),&charge[i],stypes.str().c_str());
  }
  t1->Branch("RealX",&RealX,"RealX/F"); 
  t1->Branch("RealY",&RealY,"RealY/F"); 
  t1->Branch("RealZ",&RealZ,"RealZ/F"); 
  
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  for(int i = 0; i < nEntries ; i++)
  {  
    
    tree->GetEvent(i);
    
    ExtendedTimeTag = 1e-9;
    DeltaTimeTag = 1e-9;
    
    for(int i = 0; i < 16 ; i++)
    {
      //convert to ADC channels, as if it was data from a digitizer
      //mppc gain = 1.25e6
      //adc channel binning 156e-15 C
      double adcCh = detector[i]*1.25e6*1.6e-19/156e-15;
      charge[i*2] = (Short_t) adcCh;
    }
    
    RealX = RealY = RealZ = 0;
    
    // calculate a weigthed energy deposition in x,y,z
    for(int i = 0; i < 64 ; i++) //first total energy deposited
    {
      for(int j = 0; j < px[i]->size(); j++)
      {
	RealX += (px[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
      }
      for(int j = 0; j < px[i]->size(); j++)
      {
	RealY += (py[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
      }
      for(int j = 0; j < px[i]->size(); j++)
      {
	RealZ += (pz[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
      }
    }
    
    t1->Fill();
    
    counter++;
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
    
  }
  
  std::string outFile = "Tree_OUT.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  
  t1->Write();
  
//   f1->Close();
  
  fOut->Close();
  
  
  
  return 0;
}
    