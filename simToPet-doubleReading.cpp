//simple program to translate sim output to "adc" format

// compile with 
// g++ -o ../build/simToPet simToPet.cpp `root-config --cflags --glibs`
// syntax
// simToPet `ls out*`

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TGraph.h"
#include "TF1.h"
#include "TCanvas.h"
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
  
  TChain *tree =  new TChain("tree");
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  

  //play with input names
  //std::string inputFileName = 
  
  //find the number of channels directly from the tchain file
  //before creating the variables
  //first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  //fill a vector with the leaves names
  //std::cout << nLeaves << std::endl;
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //std::cout << i << std::endl;
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  //count the entries that start with "ch"
  int numOfCh = 0;
  int numOfCry = 0;
  std::string det_prefix("detector");
  std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
      numOfCh++;
    if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
      numOfCry++;   
  }
  
  //the string "cry" appears 5 times per crystal..
  numOfCry = numOfCry / 5;
  
  std::cout << "number of channels: \t" << numOfCh << std::endl;
  std::cout << "number of crystals: \t" <<numOfCry << std::endl;
  
  
  Long64_t Seed;
  int Run;
  int Event;
  float totalEnergyDeposited;
  int NumOptPhotons;
  int NumCherenkovPhotons;
  
  std::vector<float> *CryEnergyDeposited;
  std::vector<float> **pCryEnergyDeposited;
  std::vector<float> *CryGlobalTime;
  std::vector<float> **pCryGlobalTime;
  std::vector<float> *PosXEnDep; 
  std::vector<float> **pPosXEnDep;
  std::vector<float> *PosYEnDep; 
  std::vector<float> **pPosYEnDep;
  std::vector<float> *PosZEnDep; 
  std::vector<float> **pPosZEnDep;
  
  CryEnergyDeposited  = new std::vector<float>  [numOfCry];
  pCryEnergyDeposited = new std::vector<float>* [numOfCry];
  CryGlobalTime       = new std::vector<float>  [numOfCry];
  pCryGlobalTime      = new std::vector<float>* [numOfCry];
  PosXEnDep           = new std::vector<float>  [numOfCry];
  pPosXEnDep          = new std::vector<float>* [numOfCry];
  PosYEnDep           = new std::vector<float>  [numOfCry];
  pPosYEnDep          = new std::vector<float>* [numOfCry];
  PosZEnDep           = new std::vector<float>  [numOfCry];
  pPosZEnDep          = new std::vector<float>* [numOfCry];
    
  std::vector<float> **pEdep;
  std::vector<float> **pTime;
  std::vector<float> **px;
  std::vector<float> **py;
  std::vector<float> **pz;
  
  pEdep = new std::vector<float>* [numOfCry];
  pTime = new std::vector<float>* [numOfCry];
  px    = new std::vector<float>* [numOfCry];
  py    = new std::vector<float>* [numOfCry];
  pz    = new std::vector<float>* [numOfCry];
  
  for (int i = 0 ; i < numOfCry ; i++)
  {
    pEdep[i] = 0; 
    pTime[i] = 0;
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
    snames<< "cry" << i << "GlobalTime";
    tree->SetBranchAddress(snames.str().c_str(),&pTime[i]);
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
  Short_t charge[numOfCh]; //128 channels
  Float_t RealX,RealY,RealZ;
  Short_t CrystalsHit;
  Short_t NumbOfInteractions;
  
  std::vector<float> TotalCryEnergy;	
  std::vector<float>* pTotalCryEnergy; 
  pTotalCryEnergy = &TotalCryEnergy;
  std::vector<float> DoiDetected;  
  std::vector<float>* pDoiDetected; 
  pDoiDetected = &DoiDetected;
  std::vector<float> DoiSimulation;  
  std::vector<float>* pDoiSimulation; 
  pDoiSimulation = &DoiSimulation;

  
  TTree* t1 = new TTree("adc","adc");
  
  t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
  t1->Branch("TotalCryEnergy","std::vector<float>",&pTotalCryEnergy); 
  t1->Branch("DoiDetected","std::vector<float>",&pDoiDetected); 
  t1->Branch("DoiSimulation","std::vector<float>",&pDoiSimulation); 
  //branches of the channels data
  for (int i = 0 ; i < numOfCh ; i++)
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

  //create canvas
  TCanvas* Canvas = new TCanvas("Canvas", "Canvas", 1200, 800); 


  TH2F* DOIscatter = new TH2F ("DOIscatter", "DOIscatter", 50, 0, 15, 50, 8, -8);
  DOIscatter->GetXaxis()->SetTitle("DOI Detected");
  DOIscatter->GetYaxis()->SetTitle("DOI Simulation");


  TGraph* conversionNphotonsEnergy = new TGraph();
  conversionNphotonsEnergy->SetNameTitle ("Conversion Edep vs detector hits", "Conversion Energy vs N photons detected");


  
  long int counter = 0;
  Short_t pointN = 0;
  Short_t goodCounter=0;


  int nEntries = tree->GetEntries();
  std::cout << "number of entries: \t" << nEntries << std::endl;

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
      //charge[i*2] = (Short_t) adcCh; 
      charge[i] = (Short_t) adcCh; 
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



    //prepare DoiDetected and DoiSimulation
    //fill DOIscatter plot to compare the two DOIs
    //fill conversion n photons to energy graph to find conversion factor, needed to set a minimum for accepting events
    Float_t doiDet;
    Float_t doiSim;
    Float_t CrystalLengthY = 15; //in mm
    Float_t min = 0;

    Float_t wZ = 0;
    Float_t totEnergyCry = 0;
    

    for(int i=0; i<numOfCry; i++)
    {
      if((detector[i] + detector[i+64]) != min)
      {
        for(int j=0; j<px[i]->size(); j++)
        {
          wZ += (pz[i]->at(j) * pEdep[i]->at(j));
          totEnergyCry += pEdep[i]->at(j);
        }
        doiSim = wZ/totEnergyCry;
        if (totEnergyCry != 0)
        {
          doiDet = (Float_t) CrystalLengthY*detector[i]/((Float_t) (detector[i] + detector[i+64]));
          DoiDetected.push_back(doiDet);
          DoiSimulation.push_back(doiSim);
          DOIscatter->Fill(doiDet, doiSim);
          conversionNphotonsEnergy->SetPoint(pointN, totEnergyCry, detector[i]+detector[i+64]);
          pointN++;
        }
        wZ = 0;
        totEnergyCry = 0;
      }
    }

    /*
    //work on detected data
    //find events where the energy deposited detected is > threshold in 2 crystals, and the sum of those energies is about 511 keV (1 compton scattering)

    Float_t cryThreshold = 0.1; //MeV
    Float_t totalThreshold = 0.5; //MeV
    Short_t NumbOfGoodInteractions = 0; //number of interaction where the energy detected is > threshold
    Float_t SumGoodInteractionsEnergy = 0; //sum of the energies of the good interactions
    Float_t convFactor = 5411; //comes from conversionFactor in conversionNphotonsEnergy->Fit

    for(int i=0; i<numOfCry; i++)
    {
      if(((Float_t) ((detector[i] + detector[i+64])))/convFactor > 0)
      {
        NumbOfGoodInteractions++;
        SumGoodInteractionsEnergy+= ((Float_t) (detector[i] + detector[i+64]))/convFactor;
        //put somewhere number of crystal 1 and 2
      }
    }

    if(NumbOfGoodInteractions==2 && SumGoodInteractionsEnergy>totalThreshold)
    {
      goodCounter++;
      std::cout << "\n SumGoodInteractionsEnergy: " << SumGoodInteractionsEnergy << std::endl;
      //case p1: crystal1 first, then crystal2
      //case p2: crystal2 first, then crystal1
      //factors in common will cancel out
      //do p1/p2
    }
    */


    //work on simulation data
    //find events where the energy deposited is > threshold in 2 crystals, and the sum of those energies is about 511 keV (1 compton scattering)
    Float_t cryThreshold = 0.1; //MeV
    Float_t totalThreshold = 0.5; //MeV
    Short_t NumbOfGoodInteractions = 0; //number of interaction where the energy deposited is > threshold
    Float_t SumGoodInteractionsEnergy = 0; //sum of the energies of the good interactions
    Float_t  *goodCrystals; //will store the number of crystals where good interactions happened
    goodCrystals = new Float_t [NumbOfGoodInteractions];  
    Float_t  *goodInteractionsEnergy; //will store the energy deposited in crystals where good interactions happened
    goodInteractionsEnergy = new Float_t [NumbOfGoodInteractions];  
    Float_t  *goodInteractionsDOI; //will store the DOI of the energy deposited in crystals where good interactions happened
    goodInteractionsDOI = new Float_t [NumbOfGoodInteractions];  


    for(int i=0; i<numOfCry; i++)
    {
      for(int j=0; j<px[i]->size(); j++)
      {
        totEnergyCry += pEdep[i]->at(j);
        wZ += (pz[i]->at(j) * pEdep[i]->at(j));
      }
      if(totEnergyCry > cryThreshold)
      {
        NumbOfGoodInteractions++;
        SumGoodInteractionsEnergy+= totEnergyCry;
        goodInteractionsEnergy[NumbOfGoodInteractions]=totEnergyCry;
        goodCrystals[NumbOfGoodInteractions]=i;
        goodInteractionsDOI[NumbOfGoodInteractions]=wZ/totEnergyCry;
      }
      totEnergyCry = 0;
      wZ = 0;
    }


    //filter events with exactly 2 good interactions whose energy sums up to >totalThreshold
    //the numbers of the two crystals are stored in goodCrystals[1] and goodCrystals[2]
    //the energies deposited in the two crystals are stored in goodInteractionsEnergy[1] and goodInteractionsEnergy[2]
    //the DOIs are stored in goodInteractionsDOI[1] and goodInteractionsDOI[2]
    if(NumbOfGoodInteractions==2 && SumGoodInteractionsEnergy>totalThreshold)
    {
      goodCounter++;
      //case p1: crystal1 first, then crystal2
      //case p2: crystal2 first, then crystal1
      //factors in common will cancel out
      //do p1/p2
    }







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
    }
    
    TotalCryEnergy.clear();
    DoiDetected.clear();
    DoiSimulation.clear();
  }

  std::cout << std::endl;
  std::cout << "number of good compton events: " << goodCounter << std::endl;
  std::string outFile = "Tree_OUT.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  

  DOIscatter->Write();
  
  Canvas->cd();
  TF1* line = new TF1 ("line", "[0]*x", 0, 1);
  conversionNphotonsEnergy->Fit("line", "Q");
  //conversion factor number of photons on detector and energy 
  Float_t conversionFactor = line->GetParameter(0);
  conversionNphotonsEnergy->Draw("A*");
  conversionNphotonsEnergy->GetXaxis()->SetTitle("Energy Dedposited");
  conversionNphotonsEnergy->GetYaxis()->SetTitle("N of detector hits");
  conversionNphotonsEnergy->Write();
  

  t1->Write();
  fOut->Close();
  
  return 0;
}
    