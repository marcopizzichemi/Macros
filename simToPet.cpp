// compile with
// g++ -o ../build/simToPet simToPet.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/
// syntax
// simRead_base `ls out*`

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TObjArray.h"
#include "TObject.h"
#include <algorithm>    // std::sort
#include <numeric>      // std::accumulate
#include <vector>
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"

#include "../code/struct.hh"


int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors

  //HACK to use the dictionary easily
  std::string fullFileName = "";
  // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
  std::string path = "";
  pid_t pid = getpid();
  char buf[20] = {0};
  sprintf(buf,"%d",pid);
  std::string _link = "/proc/";
  _link.append( buf );
  _link.append( "/exe");
  char proc[512];
  int ch = readlink(_link.c_str(),proc,512);
  if (ch != -1) {
    proc[ch] = 0;
    path = proc;
    std::string::size_type t = path.find_last_of("/");
    path = path.substr(0,t);
  }
  fullFileName = path + std::string("/");
  //now even worse, assuming the executable is in build and the .C macro in code
  // std::string command = ".L " + fullFileName.substr(0,fullFileName.size()-7) + "/code/structDictionary.C+";
  std::string command = ".L " + fullFileName + "structDictionary.C+";
  // std::cout << fullFileName << std::endl;
  // std::cout << "command " << command << std::endl;
  gROOT->ProcessLine(command.c_str());


  //-------------------
  // Input Files
  //-------------------
  TChain *tree =  new TChain("tree"); // read input files
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  // find the number of channels directly from the tchain file
  // before creating the variables
  // first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++)
  {
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  // int numOfCry = 0;
  std::string det_prefix("detector");
  // std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix))
    numOfCh++;
    // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
    // numOfCry++;
  }
  //the string "cry" appears 4 times per crystal..
  // numOfCry = numOfCry / 4;
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  // std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;


  //------------------
  // Input TTree
  //------------------

  //create the branches in the input ttree and connect to the variables

  // global variables
  // these are 1 number per TTree entry - so 1 number per gamma shot
  Long64_t Seed;                      // seed of the simulation (read every time, but always the same)
  int Run;                            // run id (usually just 1)(read every time, but always the same)
  int Event;                          // event id
  float totalEnergyDeposited;         // total energy deposited in this event, in all the matrix
  int NumOptPhotons;                  // number of optical photons generated in this event, in the entire matrix
  int NumCherenkovPhotons;            // number of Cherenkov photons generated in this event, in the entire matrix

  // energy deposition, each gamma 511 event has a std::vector of struct (type enDep) with all the data of each energy deposition
  std::vector<enDep> *energyDeposition = 0;

  // Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];

  // optical photons. for each gamma 511 event, every optical photon detected is a struct of type optPhot. a std::vector<optPhot> is saved for each gamma 511
  std::vector<optPhot> *photons = 0;

  //------------------------
  // Set Branch Addresses
  //------------------------
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);

  tree->SetBranchAddress("optical",&photons);
  tree->SetBranchAddress("energyDeposition",&energyDeposition);

  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }

  //output ttree

  std::string outFileName = "treeout.root"; //+ std::string(argv[1]);
  // long long int DeltaTimeTag,ExtendedTimeTag;
  // Short_t charge[32]; //adc type is always 32 channels
  // Float_t RealX,RealY,RealZ;
  // Short_t CrystalsHit;
  // Short_t NumbOfInteractions;
  //
  // TTree* t1 = new TTree("adc","adc");
  //
  // t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  // t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
  // //branches of the 32 channels data
  // for (int i = 0 ; i < 32 ; i++)
  // {
  //   //empty the stringstreams
  //   std::stringstream snames,stypes;
  //   charge[i] = 0;
  //   snames << "ch" << i;
  //   stypes << "ch" << i << "/S";
  //   t1->Branch(snames.str().c_str(),&charge[i],stypes.str().c_str());
  // }
  // t1->Branch("RealX",&RealX,"RealX/F");
  // t1->Branch("RealY",&RealY,"RealY/F");
  // t1->Branch("RealZ",&RealZ,"RealZ/F");
  // t1->Branch("CrystalsHit",&CrystalsHit,"CrystalsHit/S");
  // t1->Branch("NumbOfInteractions",&NumbOfInteractions,"NumbOfInteractions/S");
  long long int DeltaTimeTag,ExtendedTimeTag;
  Short_t charge[16];

  // correct geometry for 4x4 mppc
  // Double_t xmppc[16]={-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};
  // Double_t ymppc[16]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};

  // Double_t xmppc[16]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  // Double_t ymppc[16]={-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};

  // TH2F *flood = new TH2F("FloodHisto","FloodHisto",1000,-7,7,1000,-7,7);
  // flood->GetXaxis()->SetTitle("X [mm]");
  // flood->GetYaxis()->SetTitle("Y [mm]");
  // flood->GetZaxis()->SetTitle("N");
  //   recClean->SetTitle("Reconstruction of entire dataset");
  //   varStream << "FloodY_" << k << ":FloodX_" << k << ">>" << histoString ;

  // TH2F *positions = new TH2F("Positions","Positions",1000,-7,7,1000,-7,7);
  // positions->GetXaxis()->SetTitle("X [mm]");
  // positions->GetYaxis()->SetTitle("Y [mm]");
  // positions->GetZaxis()->SetTitle("N");


  // TH2F *HitPositions = new TH2F("HitPositions","HitPositions",1000,-7,7,1000,-7,7);
  // HitPositions->GetXaxis()->SetTitle("X [mm]");
  // HitPositions->GetYaxis()->SetTitle("Y [mm]");
  // HitPositions->GetZaxis()->SetTitle("N");
  // double columsum = 0;
  // double rowsum = 0;
  // double total = 0;
  // double floodx = 0;
  // double floody = 0;
  // double floodz = 0;


  // TBranch *b_floodx = tree->Branch("FloodX",&floodx,"floodx/F");
  // TBranch *b_floody = tree->Branch("FloodY",&floody,"floody/F");
  // TBranch *b_floodz = tree->Branch("FloodZ",&floodz,"floodz/F");

  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;

  // TH2F *averageZvsRatio = new TH2F("AverageZ vs Ratio","AverageZ vs Ratio",100,0,1,100,-8,8);
  // averageZvsRatio->GetXaxis()->SetTitle("Ratio (Fi/(Fi+Bi))");
  // averageZvsRatio->GetYaxis()->SetTitle("AverageZ [mm]");
  //
  // TH2F *normalizedZvsRatio = new TH2F("NormalizedZ vs Ratio","NormalizedZ vs Ratio",100,0,1,100,0,1);
  // normalizedZvsRatio->GetXaxis()->SetTitle("Ratio (Fi/(Fi+Bi))");
  // normalizedZvsRatio->GetYaxis()->SetTitle("NormalizedZ [fraction of crystal z]");
  //
  // TH2F *wivsRatio = new TH2F("wi vs Ratio","wi vs Ratio",100,0,1,100,0,1);
  // wivsRatio->GetXaxis()->SetTitle("Ratio (Fi/(Fi+Bi))");
  // wivsRatio->GetYaxis()->SetTitle("Wi");
  //
  // TH2F *wMeasuredvsRatioW = new TH2F("wMeasuredvsRatioW","Measured w vs. Weighted average w_0 w_1",100,0,1,100,0,1);
  // wMeasuredvsRatioW->GetXaxis()->SetTitle("RatioW = Sum_i(Wi*Ei/Etot)");
  // wMeasuredvsRatioW->GetYaxis()->SetTitle("Measured W");
  //
  // TH2F *uMeasuredVsRatioU = new TH2F("uMeasuredVsRatioU","Measured u vs. Weighted average u_0 u_1",100,-1.3,-1,100,-1.3,-1);
  // uMeasuredVsRatioU->GetXaxis()->SetTitle("RatioU = Sum_i(Ui*Ei/Etot)");
  // uMeasuredVsRatioU->GetYaxis()->SetTitle("Measured U");
  //
  // TH2F *vMeasuredVsRatioV = new TH2F("vMeasuredVsRatioV","Measured v vs. Weighted average v_0 v_1",100,1,2,100,1,2);
  // vMeasuredVsRatioV->GetXaxis()->SetTitle("RatioV = Sum_i(Vi*Ei/Etot)");
  // vMeasuredVsRatioV->GetYaxis()->SetTitle("Measured V");
  //
  // TH2F *averageZvsWi = new TH2F("AverageZ vs Wi","AverageZ vs Wi",100,0,1,100,-8,8);
  // averageZvsWi->GetXaxis()->SetTitle("Wi");
  // averageZvsWi->GetYaxis()->SetTitle("AverageZ [mm]");
  //
  // TH2F *kMeasuredvsRatioF = new TH2F("kMeasuredveRatioF","kMeasuredveRatioF",100,0,1,100,0,1);
  // kMeasuredvsRatioF->GetXaxis()->SetTitle("RatioF");
  // kMeasuredvsRatioF->GetYaxis()->SetTitle("Measured K");
  //
  // TH2F *pmaxiVsRatioF = new TH2F("pmaxiVsRatioF","Pmax_i vs. F_i",100,0,2000,100,0,2000);
  // pmaxiVsRatioF->GetXaxis()->SetTitle("F");
  // pmaxiVsRatioF->GetYaxis()->SetTitle("Pmax_i");
  //
  // TH2F *kMeasuredvsPmaxi = new TH2F("kMeasuredvsPmaxi","k vs. Ratio(Pmax_i)",100,0,1,100,0,1);
  // kMeasuredvsPmaxi->GetXaxis()->SetTitle("Ratio(Pmax_i)");
  // kMeasuredvsPmaxi->GetYaxis()->SetTitle("Measured K");



  // for(int i =0 ; i < 2; i++)
  // {
  //   std::stringstream name;
  //   name << "AverageZ vs. (Fi/(Fi+Bi)) - Crystal " << i
  //   averageZvsRatio[i] = new TH2F()
  // }
  // long int foundCandidate = 0;
  // long int singleCounter = 0;
  // long int doubleCounter = 0;
  // long int tripleCounter = 0;
  // long int multipleCounter = 0;
  for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  {
    tree->GetEvent(iEvent);

    // columsum = 0;
    // rowsum = 0;
    // total = 0;
    // floodx = 0;
    // floody = 0;
    // floodz = 0;



    //first clean the array
    // for(int i = 0; i < 16 ; i++)
    // {
    //   charge[i] = 0;
    // }

    //then fill it with the detector data
    // for(int i = 0; i < numOfCh ; i++)
    // {
    //   charge[i] = detector[i];
    // }

    //calculate weighted energy and fill 2d histo
    // for(int i = 0; i < 16 ; i++)
    // {
    //   columsum += detector[i]*ymppc[i];
    //   rowsum += detector[i]*xmppc[i];
    //   total += detector[i];
    // }
    // floodx = rowsum/total;
    // floody = columsum/total;
    // // b_floodx->Fill();
    //
    // flood->Fill(floodx,floody);




    // std::vector<int> crystals;
    // // std::vector<enDep> EventDepositions;
    //
    //
    //
    //
    // for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
    // {
    //   //read the crystal where energy was deposited
    //   int cry = energyDeposition->at(eEvent).CrystalID;
    //   //loop in the crystals found
    //   //look for the same id
    //   bool sameID = false;
    //   for(int j = 0 ; j < crystals.size(); j++)
    //   {
    //     if(crystals[j] == cry) sameID = true;
    //   }
    //   if(!sameID) crystals.push_back(cry);
    // }
    //
    // //now take only events where energy was deposited in 2 crystals and total energy deposited is 511KeV
    // //very not general, but this is only for testing: choose crystals 28 and 29 (2 crystals on the same mppc)
    // //and run only on the 9 relevant mppcs to compute w_near
    // if((crystals[0] == 28 && crystals[1] == 29) /*| (crystals[0] == 29 && crystals[1] == 28)*/)
    // // if((crystals[0] == 27 && crystals[1] == 28) | (crystals[0] == 28 && crystals[1] == 27))
    //
    // // if(true)
    // {
    //   if(totalEnergyDeposited > 0.510)
    //   {
    //     if(crystals.size() == 2)
    //     {
    //       //consider only non-lateral crystals - not needed if it's restricted to 28 and 29...
    //       bool isCandidate = true;
    //       for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and check is not on borders
    //       {
    //         bool check_b = energyDeposition->at(eEvent).CrystalI < 2 | energyDeposition->at(eEvent).CrystalI > 5 | energyDeposition->at(eEvent).CrystalJ < 2 | energyDeposition->at(eEvent).CrystalJ > 5;
    //         if(check_b)
    //         {
    //           isCandidate = false;
    //         }
    //       }
    //       if(isCandidate)
    //       {
    //         //calculate average Z of production in the 2 crystals
    //         // foundCandidate++;
    //
    //         //calculate the measured W
    //         int Measured_pmaxID = 0;
    //         int Measured_pmax = 0;
    //         int Measured_pSecond = 0;
    //         int Measured_totalDetCounts = 0;
    //         std::vector<int> vToSort;
    //         for(int i = 0; i < numOfCh ; i++)
    //         {
    //           if(i == 1 | i == 2 | i == 3 | i == 5 | i == 6 | i == 7 | i == 9 | i == 10 | i == 11) //just channel 6 (where 28 and 29 are) and surrounding
    //           {
    //             Measured_totalDetCounts+=detector[i];
    //           }
    //           vToSort.push_back(detector[i]);
    //         }
    //         std::sort (vToSort.begin(),vToSort.end());
    //         Measured_pmax = vToSort[vToSort.size()-1];
    //         Measured_pSecond = vToSort[vToSort.size()-2];
    //
    //         // {
    //         //   Measured_totalDetCounts += detector[i];
    //         //   if(detector[i] > Measured_pmax1)
    //         //   {
    //         //     Measured_pmax1 = detector[i];
    //         //     Measured_pmaxID = i;
    //         //   }
    //         // }
    //
    //
    //
    //         double Measured_w = ((double) Measured_pmax) / ((double)Measured_totalDetCounts);
    //         // double Measured_w = ((double) Measured_pmax + (double)Measured_pSecond) / (2.0*(double)Measured_totalDetCounts);
    //         // double Measured_w = ((double) Measured_pmax + (double)Measured_pSecond) / ((double)Measured_totalDetCounts);
    //         double RatioW = 0.0;
    //         double RatioF = 0.0;
    //
    //         //compute k
    //
    //         // we have u v coordinates for this event (floodx and floody)
    //         // we roughly measured the center of the two spots,
    //         // P28 = (-1.087,1.17) - P29 = (-1.087,1.97)
    //         // now we should calculate the line connecting the two spots in u,v, then the normal, then intersection and bla bla.
    //         // but we took same x for the two points, so the only coordinate that matters is y
    //         double SpotsDistance = 1.97 - 1.17;
    //         // std::cout << floodx << " "<<  floody << std::endl;
    //         double kDistance = fabs(floody - 1.97);// 28 is 0, 29 is 1 - but maybe i've inverted. whatever
    //         double Measured_k = kDistance/SpotsDistance;
    //
    //
    //         std::vector<double> totalEnergyPerCrystal;
    //         std::vector<double> averageZ;
    //         std::vector<double> normalizedZ;
    //         std::vector<double> wi;
    //         std::vector<double> ratio;
    //         std::vector<long int> Fi;
    //         std::vector<long int> Bi;
    //         std::vector<double> pmax_i;
    //         std::vector<double> ui;
    //         std::vector<double> vi;
    //
    //         for(int j = 0 ; j < crystals.size() ; j++) // we are in a specific crystal
    //         {
    //           double temp_totalEnergyPerCrystal = 0.0;
    //           for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++)
    //           {
    //
    //             if(energyDeposition->at(eEvent).CrystalID == crystals[j])
    //             temp_totalEnergyPerCrystal += energyDeposition->at(eEvent).EnergyDeposited;
    //           }
    //
    //           double temp_averageZ = 0.0;
    //           for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++)
    //           {
    //             if(energyDeposition->at(eEvent).CrystalID == crystals[j])
    //             temp_averageZ += (energyDeposition->at(eEvent).DepositionZ *  energyDeposition->at(eEvent).EnergyDeposited)/temp_totalEnergyPerCrystal;
    //           }
    //
    //           // now on opticals
    //           // 1.is averageZ correlated to Fi/(Fi+Bi)?
    //           // count Fi and Bi for this crystal
    //
    //           // 2. next, is wi correlated with ratio (actually are they the same?)?
    //           // so calculate wi
    //           // divide opticals on the basis of origin crystal
    //           // count the photons detected by each detector for both origins
    //           // compute the two wi
    //           long int temp_Fi = 0;
    //           long int temp_Bi = 0;
    //           short detectorSorted[numOfCh];
    //           for(int iDet = 0; iDet < numOfCh ;iDet++)
    //           {
    //             detectorSorted[iDet] = 0;
    //           }
    //
    //           for(int oEvent = 0; oEvent < photons->size(); oEvent++)
    //           {
    //             // calc origin crystal ID
    //             int OriginCrystalID = photons->at(oEvent).OriginCrystalI * 8 +  photons->at(oEvent).OriginCrystalJ;
    //             if(OriginCrystalID == crystals[j])
    //             {
    //               if(photons->at(oEvent).ExitFace == 1)
    //                 temp_Fi++;
    //               if(photons->at(oEvent).ExitFace == 2)
    //                 temp_Bi++;
    //               //the array of detectors
    //               // std::cout << photons->at(oEvent).PositionX << " "<< (int) ((photons->at(oEvent).PositionX + 6.3) /  3.2) << std::endl;
    //               int iDetector = (int) ((photons->at(oEvent).PositionX + 6.3) /  3.2);
    //               int jDetector = (int) ((photons->at(oEvent).PositionY + 6.3) /  3.2);
    //               int detectorID = iDetector * 4 +  jDetector;
    //               detectorSorted[detectorID]++;
    //             }
    //           }
    //
    //           int pmaxID = 0;
    //           int pmax = 0;
    //           int totalDetCounts = 0;
    //           int allDetCounts = 0;
    //           double temp_floodx = 0;
    //           double temp_floody = 0;
    //
    //           for(int i = 0; i < numOfCh ; i++)
    //           {
    //             if(i == 1 | i == 2 | i == 3 | i == 5 | i == 6 | i == 7 | i == 9 | i == 10 | i == 11) //just channel 6 (where 28 and 29 are) and surrounding
    //             {
    //               totalDetCounts += detectorSorted[i];
    //
    //             }
    //
    //             allDetCounts += detectorSorted[i];
    //             temp_floodx += xmppc[i]*detectorSorted[i];
    //             temp_floody += ymppc[i]*detectorSorted[i];
    //
    //             if(detectorSorted[i] > pmax)
    //             {
    //               pmax = detectorSorted[i];
    //               pmaxID = i;
    //             }
    //           }
    //           temp_floodx = temp_floodx/( (double) allDetCounts);
    //           temp_floody = temp_floody/( (double) allDetCounts);
    //           pmax_i.push_back(pmax);
    //           ui.push_back(temp_floodx);
    //           vi.push_back(temp_floody);
    //           double temp_wi = ((double) pmax) / ((double)totalDetCounts);
    //           double temp_ratio = ((double) temp_Fi)/((double) temp_Fi+ (double) temp_Bi);
    //           Fi.push_back(temp_Fi);
    //           Bi.push_back(temp_Bi);
    //           wi.push_back(temp_wi);
    //           totalEnergyPerCrystal.push_back(temp_totalEnergyPerCrystal);
    //           ratio.push_back(temp_ratio);
    //           averageZ.push_back(temp_averageZ);
    //           normalizedZ.push_back((temp_averageZ + 7.5) / 15.0);
    //           // RatioW += (wi * totalEnergyPerCrystal) / totalEnergyDeposited;
    //
    //         }
    //
    //         bool doubleAccepted = true;
    //         for(int countCry = 0; countCry < wi.size() ;countCry++ )
    //         {
    //           if(totalEnergyPerCrystal[countCry] < 0.05) doubleAccepted = false;
    //         }
    //         if(doubleAccepted)
    //         {
    //           foundCandidate++;
    //           RatioF = ((double) Fi[0])/(( (double) Fi[0])+( (double) Fi[1]));
    //           double RatioP = ((double) pmax_i[0])/(((double) pmax_i[0])+( (double) pmax_i[1]));
    //           double RatioU = (totalEnergyPerCrystal[0] * ui[0] + totalEnergyPerCrystal[1] * ui[1] )/ totalEnergyDeposited;
    //           double RatioV = (totalEnergyPerCrystal[0] * vi[0] + totalEnergyPerCrystal[1] * vi[1] )/ totalEnergyDeposited;
    //           for(int countCry = 0; countCry < wi.size() ;countCry++ )
    //           {
    //             averageZvsRatio->Fill(ratio[countCry],averageZ[countCry]);
    //             normalizedZvsRatio->Fill(ratio[countCry],normalizedZ[countCry]);
    //             wivsRatio->Fill(ratio[countCry],wi[countCry]);
    //             averageZvsWi->Fill(wi[countCry],averageZ[countCry]);
    //             RatioW += (wi[countCry] * totalEnergyPerCrystal[countCry]) / totalEnergyDeposited;
    //             // std::cout << pmax_i[countCry] << " " << Fi[countCry] << std::endl;
    //             pmaxiVsRatioF->Fill(Fi[countCry],pmax_i[countCry]);
    //           }
    //           wMeasuredvsRatioW->Fill(RatioW,Measured_w);
    //           kMeasuredvsRatioF->Fill(RatioF,Measured_k);
    //           kMeasuredvsPmaxi->Fill(RatioP,Measured_k);
    //           uMeasuredVsRatioU->Fill(RatioU,floodx);
    //           vMeasuredVsRatioV->Fill(RatioV,floody);
    //         }
    //
    //       }
    //     }
    //   }
    // }
    //
    // if(totalEnergyDeposited > 0.510)
    // {
    //   if(crystals.size() == 1) singleCounter++;
    //   if(crystals.size() == 2) doubleCounter++;
    //   if(crystals.size() == 3) tripleCounter++;
    //   if(crystals.size() > 3)  multipleCounter++;
    // }
    // std::cout << "Event " << iEvent << " - Dep in crystals = " << crystals.size() << std::endl;




    //calculate the average x and y deposition position
    // double averageX = 0;
    // double averageY = 0;
    // double energyColumnSum = 0;
    // double energyRowSum = 0;
    //
    // for(int i = 0; i < nCrystals ; i++)
    // {
    //   for (int j = 0 ; j < pEdep[i]->size() ; j++)
    //   {
    //     energyRowSum    += px[i]->at(j) * pEdep[i]->at(j);
    //     energyColumnSum += py[i]->at(j) * pEdep[i]->at(j);
    //   }
    // }
    // averageX = energyRowSum / totalEnergyDeposited;
    // averageY = energyColumnSum / totalEnergyDeposited ;
    //
    // positions->Fill(averageX,averageY);
    //t1->Fill();


    // ExtendedTimeTag = 1e-9;
    // DeltaTimeTag = 1e-9;
    //
    // NumbOfInteractions = 0;
    // CrystalsHit = 0;
    //
    //
    // for(int i = 0; i < numOfCh ; i++)
    // {
    //   //convert to ADC channels, as if it was data from a digitizer
    //   //mppc gain = 1.25e6
    //   //adc channel binning 156e-15 C
    //   double adcCh = detector[i]*1.25e6*1.6e-19/156e-15;
    //   charge[i*2] = (Short_t) adcCh;
    // }
    //
    // RealX = RealY = RealZ = 0;
    //
    // // calculate a weigthed energy deposition in x,y,z
    // for(int i = 0; i < numOfCry ; i++) //first total energy deposited
    // {
    //   NumbOfInteractions += px[i]->size();
    //   if(px[i]->size()) CrystalsHit++;
    //   for(int j = 0; j < px[i]->size(); j++)
    //   {
    //     RealX += (px[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
    //   }
    //   for(int j = 0; j < px[i]->size(); j++)
    //   {
    //     RealY += (py[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
    //   }
    //   for(int j = 0; j < px[i]->size(); j++)
    //   {
    //     RealZ += (pz[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
    //   }
    // }
    //
    // if(NumbOfInteractions > 0) // discard events with no energy deposition (they would never trigger the detectors anyway..)
    // {
    //   t1->Fill();
    // }

    counter++;

    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }


  }

  // std::cout << std::endl;
  // std::cout << "1 cry [511 KeV deposition] events = "    << singleCounter << std::endl;
  // std::cout << "2 cry [511 KeV deposition] events = "    << doubleCounter << std::endl;
  // std::cout << "3 cry [511 KeV deposition] events = "    << tripleCounter << std::endl;
  // std::cout << "Multi cry [511 KeV deposition] events = "<< multipleCounter << std::endl;
  // std::cout << "Candidates = "<< foundCandidate << std::endl;
  //
  // TF1 *line = new TF1("line","x",-7,7);
  // line->SetLineColor(kRed);
  // TF1 *line2 = new TF1("line2","x",0,2000);
  // line2->SetLineColor(kRed);
  //
  //
  //
  // std::string outFile = "FileOut.root";
  // TFile* fOut = new TFile(outFile.c_str(),"recreate");
  //
  // TCanvas *C_flood = new TCanvas("C_flood","C_flood",800,800);
  // flood->Draw("COLZ");
  // C_flood->Write();
  // TCanvas *C_averageZvsRatio = new TCanvas("C_averageZvsRatio","C_averageZvsRatio",800,800);
  // C_averageZvsRatio->cd();
  // averageZvsRatio->Draw();
  // C_averageZvsRatio->Write();
  //
  // TCanvas *C_normalizedZvsRatio = new TCanvas("C_normalizedZvsRatio","C_normalizedZvsRatio",800,800);
  // C_normalizedZvsRatio->cd();
  // normalizedZvsRatio->Draw();
  // C_normalizedZvsRatio->Write();
  //
  // TCanvas *C_averageZvsWi = new TCanvas("C_averageZvsWi","C_averageZvsWi",800,800);
  // C_averageZvsWi->cd();
  // averageZvsWi->Draw();
  // C_averageZvsWi->Write();
  //
  // TCanvas *C_wivsRatio = new TCanvas("C_wivsRatio","C_wivsRatio",800,800);
  // C_wivsRatio->cd();
  // wivsRatio->Draw();
  // line->Draw("same");
  // C_wivsRatio->Write();
  //
  // TCanvas *C_wMeasuredvsRatioW = new TCanvas("C_wMeasuredvsRatioW","C_wMeasuredvsRatioW",800,800);
  // C_wMeasuredvsRatioW->cd();
  // wMeasuredvsRatioW->Draw();
  // line->Draw("same");
  // C_wMeasuredvsRatioW->Write();
  //
  // TCanvas *C_uMeasuredVsRatioU = new TCanvas("C_uMeasuredVsRatioU","C_uMeasuredVsRatioU",800,800);
  // C_uMeasuredVsRatioU->cd();
  // uMeasuredVsRatioU->Draw();
  // line->Draw("same");
  // C_uMeasuredVsRatioU->Write();
  //
  // TCanvas *C_vMeasuredVsRatioV = new TCanvas("C_vMeasuredVsRatioV","C_vMeasuredVsRatioV",800,800);
  // C_vMeasuredVsRatioV->cd();
  // vMeasuredVsRatioV->Draw();
  // line->Draw("same");
  // C_vMeasuredVsRatioV->Write();
  //
  //
  //
  //
  // TCanvas *C_kMeasuredvsRatioF = new TCanvas("C_kMeasuredvsRatioF","C_kMeasuredvsRatioF",800,800);
  // C_kMeasuredvsRatioF->cd();
  // kMeasuredvsRatioF->Draw();
  // line->Draw("same");
  // C_kMeasuredvsRatioF->Write();
  //
  // TCanvas *C_kMeasuredvsPmaxi = new TCanvas("C_kMeasuredvsPmaxi","C_kMeasuredvsPmaxi",800,800);
  // C_kMeasuredvsPmaxi->cd();
  // kMeasuredvsPmaxi->Draw();
  // line->Draw("same");
  // C_kMeasuredvsPmaxi->Write();
  //
  // TCanvas *C_pmaxiVsRatioF = new TCanvas("C_pmaxiVsRatioF","C_pmaxiVsRatioF",800,800);
  // C_pmaxiVsRatioF->cd();
  // pmaxiVsRatioF->Draw();
  // line2->Draw("same");
  // C_pmaxiVsRatioF->Write();


  // t1->Write();

  //   f1->Close();

  // std::string outFile = "analysis_OUT.root";
  // TFile* fOut = new TFile(outFile.c_str(),"recreate");

  // flood->Write();
  // positions->Write();
  //   f1->Close();

  TFile* fOut = new TFile(outFileName.c_str(),"recreate");
  t1->Write();
  fOut->Close();



  return 0;
}
