//simple analysis

// compile with 
// g++ -o doiSimAnalysis doiSimAnalysis.cpp `root-config --cflags --glibs`


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TChain.h"
#include "TCanvas.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>

int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  std::string baseFileName = argv[1];
  std::string baseOutFile = "spot_" + baseFileName.substr(0,baseFileName.length()-5);
  std::string baseOutFileRoot = baseOutFile + ".root";
  
  int nCrystalsX = 1;
  int nCrystalsY = 2;
  int nDetectorsX = 1;
  int nDetectorsY = 1;
  
  long int            Seed;
  Int_t               Run;
  Int_t               Event;
  Int_t               NumOptPhotons;
  Int_t               NumCherenkovPhotons;
  Float_t             totalEnergyDeposited;
  
  
  Short_t*            DetectorHit;          
  //   int            DetectorHit[16];
  
  Int_t               TagNumInteractions;
  Int_t               CryNumInteractions;
  Float_t             CryTotalEnergyDeposited;
  Float_t             TagTotalEnergyDeposited;
  Float_t             GeneratedSourceX;
  Float_t             GeneratedSourceY;
  Float_t             GeneratedSourceZ;
  Float_t             GeneratedSourcePhi;
  Float_t	      GeneratedSourceTheta;
  
  std::vector<float>*  CryEnergyDeposited;   
  std::vector<float>** pCryEnergyDeposited;
  
  std::vector<float>*  PosXEnDep; 
  std::vector<float>** pPosXEnDep;
  std::vector<float>*  PosYEnDep; 
  std::vector<float>** pPosYEnDep;
  std::vector<float>*  PosZEnDep; 
  std::vector<float>** pPosZEnDep;
  
  std::vector<float> PositionX;
  std::vector<float> *pPositionX;
  std::vector<float> PositionY;
  std::vector<float> *pPositionY;
  std::vector<float> PositionZ;
  std::vector<float> *pPositionZ;
  
  std::vector<float> PreMomentumX;
  std::vector<float> *pPreMomentumX;
  std::vector<float> PreMomentumY;
  std::vector<float> *pPreMomentumY;
  std::vector<float> PreMomentumZ;
  std::vector<float> *pPreMomentumZ;
  
  std::vector<float> PostMomentumX;
  std::vector<float> *pPostMomentumX;
  std::vector<float> PostMomentumY;
  std::vector<float> *pPostMomentumY;
  std::vector<float> PostMomentumZ;
  std::vector<float> *pPostMomentumZ;
  
  std::vector<float> GlobalTime;
  std::vector<float> *pGlobalTime;
  
  std::vector<int> PhotonType;
  std::vector<int> *pPhotonType;
  
  std::vector<float> PhotonEnergy;
  std::vector<float> *pPhotonEnergy;
  
  
  DetectorHit         = new Short_t [nDetectorsX*nDetectorsY];
  CryEnergyDeposited  = new std::vector<float> [nCrystalsX*nCrystalsY];
  pCryEnergyDeposited = new std::vector<float>* [nCrystalsX*nCrystalsY];
  PosXEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  pPosXEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  PosYEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  pPosYEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  PosZEnDep           = new std::vector<float> [nCrystalsX*nCrystalsY];
  pPosZEnDep          = new std::vector<float>* [nCrystalsX*nCrystalsY];
  
  
  for(int i = 0; i < nCrystalsX*nCrystalsY ; i++) 
  {
    pCryEnergyDeposited[i] = &CryEnergyDeposited[i];
    pPosXEnDep[i] = &PosXEnDep[i];
    pPosYEnDep[i] = &PosYEnDep[i];
    pPosZEnDep[i] = &PosZEnDep[i];
  }
  
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  
  tree->SetBranchAddress("TagNumInteractions"     ,&TagNumInteractions     );
  tree->SetBranchAddress("CryNumInteractions"     ,&CryNumInteractions     );
  tree->SetBranchAddress("CryTotalEnergyDeposited",&CryTotalEnergyDeposited);
  tree->SetBranchAddress("TagTotalEnergyDeposited",&TagTotalEnergyDeposited);
  tree->SetBranchAddress("GeneratedSourceX"       ,&GeneratedSourceX       );
  tree->SetBranchAddress("GeneratedSourceY"       ,&GeneratedSourceY       );
  tree->SetBranchAddress("GeneratedSourceZ"       ,&GeneratedSourceZ       );
  tree->SetBranchAddress("GeneratedSourcePhi"     ,&GeneratedSourcePhi     );
  tree->SetBranchAddress("GeneratedSourceTheta"   ,&GeneratedSourceTheta   );
  
  
  
  for (int i = 0 ; i < nCrystalsX*nCrystalsY ; i++)
  {
    std::stringstream snames;
    snames << "cry" << i;
    tree->SetBranchAddress(snames.str().c_str(),&pCryEnergyDeposited[i]);
    snames.str("");
    snames<< "cry" << i << "PosXEnDep";    
    tree->SetBranchAddress(snames.str().c_str(),&pPosXEnDep[i]);
    snames.str("");
    snames<< "cry" << i << "PosYEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&pPosYEnDep[i]);
    snames.str("");
    snames<< "cry" << i << "PosZEnDep";
    tree->SetBranchAddress(snames.str().c_str(),&pPosZEnDep[i]);
  }
  for (int i = 0 ; i < nDetectorsX*nDetectorsY ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&DetectorHit[i]);
  }
  
  
  
  
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  
  
  TH1F *CryX = new TH1F("CryX","CryX",100,-8,8);
  CryX->GetXaxis()->SetTitle("Crystal Length [mm]");
  CryX->SetTitle("DOI bench - beam profile");
  CryX->GetYaxis()->SetTitle("N");
  
  TH1F *CryY = new TH1F("CryY","CryY",100,-8,8);
  CryY->GetXaxis()->SetTitle("Crystal Length [mm]");
  CryY->SetTitle("DOI bench - beam profile");
  CryY->GetYaxis()->SetTitle("N");
  
  TH1F *CryZ = new TH1F("CryZ","CryZ",100,-8,8);
  CryZ->GetXaxis()->SetTitle("Crystal Length [mm]");
  CryZ->SetTitle("DOI bench - beam profile");
  CryZ->GetYaxis()->SetTitle("N");
  
  
  for(int i = 0; i < nEntries ; i++)
  {  
    tree->GetEvent(i);
    
    Double_t RealX = 0;
    Double_t RealY = 0;
    Double_t RealZ = 0;
    Double_t endep = 0;
    
    if(TagTotalEnergyDeposited > 0.5)
    {
      if(CryTotalEnergyDeposited > 0.5)
      {
	for(int j = 0; j < pPosXEnDep[1]->size(); j++)
	{
	  
	  RealX += (pPosXEnDep[1]->at(j) * pCryEnergyDeposited[1]->at(j))/CryTotalEnergyDeposited;
	}
	for(int j = 0; j < pPosXEnDep[1]->size(); j++)
	{
	  RealY += (pPosYEnDep[1]->at(j) * pCryEnergyDeposited[1]->at(j))/CryTotalEnergyDeposited;
	}
	for(int j = 0; j < pPosXEnDep[1]->size(); j++)
	{
	  RealZ += (pPosZEnDep[1]->at(j) * pCryEnergyDeposited[1]->at(j))/CryTotalEnergyDeposited;
	}
	CryX->Fill(RealX);
	CryY->Fill(RealY);
	CryZ->Fill(RealZ);
      }
    } 
    
    
    
    
    
    
    
    
    
    
    
    
  }
  
//   TCanvas *C_CryX = new TCanvas("C_CryX","",800,600);
//   C_CryX->cd();
//   CryX->Draw();
//   C_CryX->Print("spotX.gif");
//   TCanvas *C_CryY = new TCanvas("C_CryY","",800,600);
//   C_CryY->cd();
//   CryY->Draw();
//   C_CryY->Print("spotY.gif");
//   TCanvas *C_CryZ = new TCanvas("C_CryZ","",800,600);
//   C_CryZ->cd();
//   CryZ->Draw();
//   C_CryZ->Print("spotZ.gif");
  
  TFile* fOut = new TFile(baseOutFileRoot.c_str(),"recreate");
  
  CryX->Write();
  CryY->Write();
//   CryZ->Write();
  
  return 0;
}