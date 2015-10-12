//simple analysis

// compile with 
// g++ -o simAnalysis simAnalysis.cpp `root-config --cflags --glibs`


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TChain.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main (int argc, char** argv)
{
  
  int nCrystals = 1;
  int nDetectors = 16;
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
  
  std::vector<float> CryEnergyDeposited[nCrystals]; 
  std::vector<float> *pCryEnergyDeposited[nCrystals];
  std::vector<float> PosXEnDep[nCrystals]; 
  std::vector<float> *pPosXEnDep[nCrystals];
  std::vector<float> PosYEnDep[nCrystals]; 
  std::vector<float> *pPosYEnDep[nCrystals];
  std::vector<float> PosZEnDep[nCrystals]; 
  std::vector<float> *pPosZEnDep[nCrystals];
  short RunDetectorHit[nDetectors];
  
  
  std::vector<float> *pEdep[nCrystals];
  std::vector<float> *px[nCrystals];
  std::vector<float> *py[nCrystals];
  std::vector<float> *pz[nCrystals];
  
  for (int i = 0 ; i < nCrystals ; i++)
  {
    pEdep[i] = 0; 
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
  }
  Short_t  detector[nDetectors];
  
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  
  for (int i = 0 ; i < nCrystals ; i++)
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
  for (int i = 0 ; i < nDetectors ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  
  
  
  
  
  //output ttree
  
  long long int DeltaTimeTag,ExtendedTimeTag;
  Short_t charge[32];
  
//   TTree* t1 = new TTree("adc","adc");
//   
//   t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
//   t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
//   //branches of the 32 channels data
//   for (int i = 0 ; i < 32 ; i++)
//   {
//     //empty the stringstreams
//     std::stringstream snames,stypes;
//     charge[i] = 0;
//     snames << "ch" << i;
//     stypes << "ch" << i << "/S";  
//     t1->Branch(snames.str().c_str(),&charge[i],stypes.str().c_str());
//   }
  
  // correct one
  Double_t xmppc[32]={-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0};
  Double_t ymppc[32]={-4.8,0,-4.8,0,-4.8,0,-4.8,0,-1.6,0,-1.6,0,-1.6,0,-1.6,0,1.6,0,1.6,0,1.6,0,1.6,0,4.8,0,4.8,0,4.8,0,4.8,0};
  
  // wrong geometry
//   Double_t xmppc[32]={-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0};
//   Double_t ymppc[32]={-4.65,0,-4.65,0,-4.65,0,-4.65,0,-1.55,0,-1.55,0,-1.55,0,-1.55,0,1.55,0,1.55,0,1.55,0,1.55,0,4.65,0,4.65,0,4.65,0,4.65,0};
  
  
  TH2F *flood = new TH2F("FloodHisto","FloodHisto",1000,-7,7,1000,-7,7);
  flood->GetXaxis()->SetTitle("X [mm]");
  flood->GetYaxis()->SetTitle("Y [mm]");
  flood->GetZaxis()->SetTitle("N");
//   recClean->SetTitle("Reconstruction of entire dataset");
//   varStream << "FloodY_" << k << ":FloodX_" << k << ">>" << histoString ;
  
  TH2F *positions = new TH2F("Positions","Positions",1000,-7,7,1000,-7,7);
  positions->GetXaxis()->SetTitle("X [mm]");
  positions->GetYaxis()->SetTitle("Y [mm]");
  positions->GetZaxis()->SetTitle("N");
  
  
  TH2F *HitPositions = new TH2F("HitPositions","HitPositions",1000,-7,7,1000,-7,7);
  HitPositions->GetXaxis()->SetTitle("X [mm]");
  HitPositions->GetYaxis()->SetTitle("Y [mm]");
  HitPositions->GetZaxis()->SetTitle("N");
//
  
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  for(int iEntry = 0; iEntry < nEntries ; iEntry++)
  {  
    tree->GetEvent(iEntry);
    double columsum = 0; 
    double rowsum = 0;
    double total = 0;
    double floodx = 0;
    double floody = 0;
    
    //first clean the array
    for(int i = 0; i < 32 ; i++)
    {
      charge[i] = 0;
    }
    
    //then fill it with the detector data
    for(int i = 0; i < nDetectors ; i++)
    {
      charge[i*2] = detector[i];
    }
    
    //calculate weighted energy and fill 2d histo
    for(int i = 0; i < 32 ; i++)
    {
      columsum += charge[i]*ymppc[i];
      rowsum += charge[i]*xmppc[i];
      total += charge[i];
    }
    floodx = rowsum/total;
    floody = columsum/total;
    flood->Fill(floodx,floody);
    
    
    //calculate the average x and y deposition position 
    double averageX = 0;
    double averageY = 0;
    double energyColumnSum = 0;
    double energyRowSum = 0;
    
    for(int i = 0; i < nCrystals ; i++)
    {
      for (int j = 0 ; j < pEdep[i]->size() ; j++)
      {
        energyRowSum    += px[i]->at(j) * pEdep[i]->at(j);
        energyColumnSum += py[i]->at(j) * pEdep[i]->at(j);
      }
    }
    averageX = energyRowSum / totalEnergyDeposited;
    averageY = energyColumnSum / totalEnergyDeposited ;
    
    positions->Fill(averageX,averageY);
    //t1->Fill();
    
    counter++;
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
    
  }
  
  std::string outFile = "analysis_OUT.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  
  flood->Write();
  positions->Write();
//   f1->Close();
  
  fOut->Close();
  
  
  
  return 0;
}
    