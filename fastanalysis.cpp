// g++ -o fastanalysis fastanalysis.cpp `root-config --cflags --glibs`


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
#include "TEllipse.h"
#include "TGraph2D.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TChain.h"
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 


int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>");
  
  //TFile f1(argv[1]);
  
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  
  
  //TTree *tree = (TTree*)f1.Get("tree");
  
  long int Seed;
  int Run;
  int Event;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  Short_t charge[32];
  Short_t hit[16];
  
  std::vector<float> CryEnergyDeposited[1]; 
  std::vector<float> *pCryEnergyDeposited[1];
  std::vector<float> PosXEnDep[1]; 
  std::vector<float> *pPosXEnDep[1];
  std::vector<float> PosYEnDep[1]; 
  std::vector<float> *pPosYEnDep[1];
  std::vector<float> PosZEnDep[1]; 
  std::vector<float> *pPosZEnDep[1];
  short RunDetectorHit[16];
  
  
  std::vector<float> *pEdep[1];
  std::vector<float> *px[1];
  std::vector<float> *py[1];
  std::vector<float> *pz[1];
  
    
  std::vector<float> TransmissionX;
  std::vector<float> TransmissionY;
  std::vector<float> TransmissionZ;
  std::vector<float> *pTransmissionX;
  std::vector<float> *pTransmissionY;
  std::vector<float> *pTransmissionZ;
  
  std::vector<float> PositionX;
  std::vector<float> *pPositionX;
  std::vector<float> PositionY;
  std::vector<float> *pPositionY;
  std::vector<float> *pGlobalTime;
  
  
  for (int i = 0 ; i < 1 ; i++)
  {
    pEdep[i] = 0; 
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
    pTransmissionX = 0;
    pTransmissionY = 0;
    pTransmissionZ = 0;
    pPositionX = 0;
    pPositionY = 0;
    pGlobalTime = 0;
  }
  Short_t detector[16];
  
  tree->Branch("Seed",&Seed,"Seed/L");
  tree->Branch("Run",&Run,"Run/I");
  tree->Branch("Event",&Event,"Event/I");
  tree->Branch("totalEnergyDeposited",&totalEnergyDeposited,"totalEnergyDeposited/F");
  tree->Branch("NumOptPhotons",&NumOptPhotons,"NumOptPhotons/I");
  tree->Branch("NumCherenkovPhotons",&NumCherenkovPhotons,"NumCherenkovPhotons/I");
  
//   for (int i = 0 ; i < 1 ; i++)
//   {
//     std::stringstream snames;
//     snames << "cry" << i;
//     tree->SetBranchAddress(snames.str().c_str(),&pEdep[i]);
//     snames.str("");
//     snames<< "cry" << i << "PosXEnDep";    
//     tree->SetBranchAddress(snames.str().c_str(),&px[i]);
//     snames.str("");
//     snames<< "cry" << i << "PosYEnDep";
//     tree->SetBranchAddress(snames.str().c_str(),&py[i]);
//     snames.str("");
//     snames<< "cry" << i << "PosZEnDep";
//     tree->SetBranchAddress(snames.str().c_str(),&pz[i]);
//   }
  for (int i = 0 ; i < 16 ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  
//   tree->SetBranchAddress("TransmissionX",&pTransmissionX);
//   tree->SetBranchAddress("TransmissionY",&pTransmissionY);
//   tree->SetBranchAddress("TransmissionZ",&pTransmissionZ);
  
  tree->SetBranchAddress("PositionX",&pPositionX);
  tree->SetBranchAddress("PositionY",&pPositionY);
  
  tree->SetBranchAddress("GlobalTime",&pGlobalTime);
  
  int nEntries = tree->GetEntries();
  
  TH1F *histoLR = new TH1F("histoLR","histoLR",20,-10,10);
  TH1F *histoTB = new TH1F("histoTB","histoTB",20,-10,10);
  TH2F *flood = new TH2F("flood","flood",1000,-7,7,1000,-7,7);
  flood->GetXaxis()->SetTitle("X [mm]");
  flood->GetYaxis()->SetTitle("Y [mm]");
  flood->GetZaxis()->SetTitle("N");
  
  
  TH2F *HitPositions = new TH2F("HitPositions","HitPositions",1000,-7,7,1000,-7,7);
  HitPositions->GetXaxis()->SetTitle("X [mm]");
  HitPositions->GetYaxis()->SetTitle("Y [mm]");
  HitPositions->GetZaxis()->SetTitle("N");
  
  TH2F *AvgPositions = new TH2F("AvgPositions","AvgPositions",1000,-7,7,1000,-7,7);
  AvgPositions->GetXaxis()->SetTitle("X [mm]");
  AvgPositions->GetYaxis()->SetTitle("Y [mm]");
  AvgPositions->GetZaxis()->SetTitle("N");
  
  TH2F *floodHit = new TH2F("floodHit","floodHit",1000,-7,7,1000,-7,7);
  floodHit->GetXaxis()->SetTitle("X [mm]");
  floodHit->GetYaxis()->SetTitle("Y [mm]");
  floodHit->GetZaxis()->SetTitle("N");
  
  //TGraph2D graph = new TGraph2D(
  
  // correct one
  Double_t xmppc[32]={-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0,-4.8,0,-1.6,0,1.6,0,4.8,0};
  Double_t ymppc[32]={-4.8,0,-4.8,0,-4.8,0,-4.8,0,-1.6,0,-1.6,0,-1.6,0,-1.6,0,1.6,0,1.6,0,1.6,0,1.6,0,4.8,0,4.8,0,4.8,0,4.8,0};
  
  Short_t fakeMppc[4][4];
  double xlimit[8] = {-6.3,-3.3,-3.1,-0.1,0.1,3.1,3.3,6.3};
  double ylimit[8] = {-6.3,-3.3,-3.1,-0.1,0.1,3.1,3.3,6.3};
  
  long long int counter = 0;
   
  TGraph2D *graph = new TGraph2D();
  
  //std::cout << "qui" << std::endl;
  for(int i = 0; i < nEntries ; i++)
  {    
    tree->GetEvent(i);
//     std::cout << "---------------------------------------------------------------------------------" <<std::endl;
    int hitCounter0 = 0;
    for(int ix = 0; ix < 4 ; ix++)
    {
      for(int iy = 0; iy < 4 ; iy++)
      {  
	fakeMppc[ix][iy] = 0;
	hit[hitCounter0] = 0;
	hitCounter0++;
      }
    }
//     int nPoints = 0;
    int left = 0;
    int right = 0;
    int bottom = 0;
    int top = 0;
//     for (int i = 0 ; i < pTransmissionX->size() ; i++)
//     {
//       if(pTransmissionX->at(i) < -6.35)
// 	left++;
//       if(pTransmissionX->at(i) > 6.35)
// 	right++;
//       if(pTransmissionY->at(i) < -6.35)
// 	bottom++;
//       if(pTransmissionY->at(i) > 6.35)
// 	top++;
//       
//       
//     }
    //std::cout  << left << " " << right << std::endl;
//     if(left != 0 && right != 0) histoLR->Fill(right - left);
//     if(top != 0 && bottom != 0) histoTB->Fill(top - bottom);
    
    double avgX = 0;
    double avgY = 0;
   
    
    for(int i = 0 ; i < pPositionX->size() ; i++)
    {
      HitPositions->Fill(pPositionX->at(i),pPositionY->at(i));
      
      for(int ix = 0; ix < 4 ; ix++)
      {
	for(int iy = 0; iy < 4 ; iy++)
	{
	  if( (pPositionX->at(i) > xlimit[ix*2]) && (pPositionX->at(i) < xlimit[ix*2+1]) && (pPositionY->at(i) > ylimit[iy*2]) && (pPositionY->at(i) < ylimit[iy*2+1]) )
	    fakeMppc[ix][iy]++;
	}
      }
      
//       std::cout <<  pPositionX->at(i) << "\t";
//       std::cout <<  pPositionY->at(i)<< "\t";
//       std::cout <<  pGlobalTime->at(i)<< std::endl;
      
      //graph->SetPoint(i,posx,posy,gtime);
      
      avgX += pPositionX->at(i);
      avgY += pPositionY->at(i);
    }
    
    avgX = avgX / pPositionX->size();
    avgY = avgY / pPositionY->size();
    
    AvgPositions->Fill(avgX,avgY);
    
    
    double columsum = 0; 
    double rowsum = 0;
    double total = 0;
    double floodx = 0;
    double floody = 0;
    
    for(int i = 0; i < 32 ; i++)
    {
      charge[i] = 0;
    }
    
    //then fill it with the detector data
    for(int i = 0; i < 16 ; i++)
    {
      charge[i*2] = detector[i];
//       std::cout << detector[i] << " ";
    }
//     std::cout << std::endl;
    
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
    
    int hitCounter = 0;
    
    for(int i = 0; i < 4 ; i++)
    {
      for(int j = 0; j < 4 ; j++)
      {
	hit[hitCounter] = fakeMppc[i][j];
// 	std::cout << fakeMppc[i][j] << " ";
	hitCounter++;
      }
    }
//     std::cout << std::endl;
    double Hitcolumsum = 0; 
    double Hitrowsum = 0;
    double Hittotal = 0;
    double Hitfloodx = 0;
    double Hitfloody = 0;
    
    
    for(int i = 0; i < 16 ; i++)
    {
      Hitcolumsum += hit[i]*ymppc[i*2];
      Hitrowsum += hit[i]*xmppc[i*2];
      Hittotal += hit[i];
    }
    Hitfloodx = Hitrowsum/Hittotal;
    Hitfloody = Hitcolumsum/Hittotal;
    floodHit->Fill(Hitfloodx,Hitfloody);
    
//     if(!pEdep[1]->empty()) //if there is energy deposited in this crystal
//     {
//       //calculate total energy dep in crystal
//       float totalEnergy = 0;
//       
//       for(int iCry = 0 ; iCry < pEdep[1]->size() ; iCry++)
//       {
// 	totalEnergy += pEdep[1]->at(iCry);
//       }
//       if( totalEnergy > 0.510 /*&& pz[1]->at(0) > -2.0 && pz[1]->at(0) < 2.0*/ ) //if all energy was deposited in this crystal and if it was deposited in the centre
//       {
// 	histo->Fill(detector[1]);
// 	int sum = 0;
// 	for(int j = 0 ; j < 1 ; j++)
// 	{
// 	  sum += detector[j];
// 	}
// 	histoSum->Fill(sum);
//       }
//       if( totalEnergy > 0.510)
//       {
// 	float weigtheZPos = 0;
// 	for(int iCry = 0 ; iCry < pEdep[1]->size() ; iCry++)
//         {
// 	  weigtheZPos += pEdep[1]->at(iCry)*pz[1]->at(iCry);
//         }
//         weigtheZPos = weigtheZPos / totalEnergy;
// 	scatter->Fill(weigtheZPos,detector[1]);
//       }
    
    
    
    
    counter++;
    
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  } 
  
  
  
  
  //std::cout << histo->GetMean() << std::endl;
  
  
  std::string outFile = "fast_out.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  fOut->cd();
//   histoLR->Write();
//   histoTB->Write();
  HitPositions->Write();
  flood->Write();
  AvgPositions->Write();
  
  floodHit->Write();
  //graph->Write();
//   histoSum->Write();
//   scatter->Write();
  fOut->Close();
  
  
  
  
  return 0;
   
}