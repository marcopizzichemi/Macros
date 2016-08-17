//simple program to translate sim output to "adc" format

// compile with 
// g++ -o ../build/simToPet simToPet.cpp `root-config --cflags --glibs`
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
#include "TObjArray.h"
#include "TObject.h"
#include <algorithm>

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
  
  //play with input names
//   std::string inputFileName = 
  
  // find the number of channels directly from the tchain file
  // before creating the variables
  // first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
//   std::cout << nLeaves << std::endl;
  for(int i = 0 ; i < nLeaves ; i++)
  {
//     std::cout << i << std::endl;
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  int numOfCry = 0;
  std::string det_prefix("detector");
  std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
//     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
      numOfCh++;
    if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
      numOfCry++;
    
  }
  
  //the string "cry" appears 4 times per crystal..
  numOfCry = numOfCry / 4;
  
  std::cout << numOfCh << std::endl;
  std::cout << numOfCry << std::endl;
  
  
  Long64_t Seed;
  int Run;
  int Event;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  std::vector<float> *CryEnergyDeposited;
  std::vector<float> **pCryEnergyDeposited;
  std::vector<float> *PosXEnDep; 
  std::vector<float> **pPosXEnDep;
  std::vector<float> *PosYEnDep; 
  std::vector<float> **pPosYEnDep;
  std::vector<float> *PosZEnDep; 
  std::vector<float> **pPosZEnDep;
  
//   DetectorHit         = new Short_t             [numOfCh];
  CryEnergyDeposited  = new std::vector<float>  [numOfCry];
  pCryEnergyDeposited = new std::vector<float>* [numOfCry];
  PosXEnDep           = new std::vector<float>  [numOfCry];
  pPosXEnDep          = new std::vector<float>* [numOfCry];
  PosYEnDep           = new std::vector<float>  [numOfCry];
  pPosYEnDep          = new std::vector<float>* [numOfCry];
  PosZEnDep           = new std::vector<float>  [numOfCry];
  pPosZEnDep          = new std::vector<float>* [numOfCry];
  
//   short RunDetectorHit[16];
  
  
  std::vector<float> **pEdep;
  std::vector<float> **px;
  std::vector<float> **py;
  std::vector<float> **pz;
  
  pEdep = new std::vector<float>* [numOfCry];
  px    = new std::vector<float>* [numOfCry];
  py    = new std::vector<float>* [numOfCry];
  pz    = new std::vector<float>* [numOfCry];
  
  for (int i = 0 ; i < numOfCry ; i++)
  {
    pEdep[i] = 0; 
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
  }

  Short_t  *detector;
  detector = new Short_t [numOfCh];
  
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  
  for (int i = 0 ; i < numOfCry ; i++)
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
  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  
  
  
  //output ttree
  long long int DeltaTimeTag,ExtendedTimeTag;
  Short_t charge[32]; //adc type is always 32 channels
  Float_t RealX,RealY,RealZ;
  Short_t CrystalsHit;
  Short_t NumbOfInteractions;
  std::vector <float> TotalCryEnergy;	
  std::vector <float>* pTotalCryEnergy; 
  pTotalCryEnergy = &TotalCryEnergy;
  
  TTree* t1 = new TTree("adc","adc");
  
  t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
  t1->Branch("TotalCryEnergy","std::vector<float>",&pTotalCryEnergy);	
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
  t1->Branch("CrystalsHit",&CrystalsHit,"CrystalsHit/S"); 
  t1->Branch("NumbOfInteractions",&NumbOfInteractions,"NumbOfInteractions/S"); 

  
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  for(int i = 0; i < nEntries ; i++)
  {  
    
    tree->GetEvent(i);
    
    ExtendedTimeTag = 1e-9;
    DeltaTimeTag = 1e-9;
    
    NumbOfInteractions = 0;
    CrystalsHit = 0;

    
    
    for(int i = 0; i < numOfCh ; i++)
    {
      //convert to ADC channels, as if it was data from a digitizer
      //mppc gain = 1.25e6
      //adc channel binning 156e-15 C
      double adcCh = detector[i]*1.25e6*1.6e-19/156e-15;
      charge[i*2] = (Short_t) adcCh; 
    }
    
    RealX = RealY = RealZ = 0;
    
    // calculate a weigthed energy deposition in x,y,z
    for(int i = 0; i < numOfCry ; i++) //first total energy deposited
    {
      Float_t SumEnergy = 0;
      NumbOfInteractions += px[i]->size();
      if(px[i]->size()) CrystalsHit++;
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
      for(int j = 0; j < px[i]->size(); j++)
      {
		SumEnergy += pEdep[i]->at(j);	
      }
      TotalCryEnergy.push_back(SumEnergy);
    }
/*
    //find crystal with max energy deposition
    Float_t MaxEnergyCry = 0;
    Short_t MaxEnergyCryNum = -1;
    MaxEnergyCry = *std::max_element(TotalCryEnergy.begin(), TotalCryEnergy.end());
    MaxEnergyCryNum = std::distance(TotalCryEnergy.begin(), (std::max_element(TotalCryEnergy.begin(), TotalCryEnergy.end())));
*/
    
    if(NumbOfInteractions > 0) // discard events with no energy deposition (they would never trigger the detectors anyway..)
    {
      t1->Fill();
    }
    
    counter++;
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
  TotalCryEnergy.clear();
  
  }
  std::cout << std::endl;
  std::string outFile = "Tree_OUT.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  
  t1->Write();
  
//   f1->Close();
  
  fOut->Close();
  
  
  
  return 0;
}
    