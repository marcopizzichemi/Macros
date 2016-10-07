// compile with
// g++ -o ../build/simRead_struct simRead_struct.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/
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

  //energy deposition events
  // these are 1 array of std::vector<float> per TTree entry. Explanation:
  // one TTree entry is 1 gamma shot, but
  // since the energy deposition events per gamma shot can be whatever, from 0 to N
  // we saved the data in std::vectors since they have variable size
  // if in a particular sim event the gamma deposited energy in N different points
  // the std::vectors are effectively arrays of length N
  // holding info of each energy deposition event, in this way
  // pEdep[N] = E0,E1,E2,E3,...,EN
  // px[N]    = x0,x1,x2,x3,...,xN
  // py[N]    = y0,y1,y2,y3,...,yN
  // pz[N]    = z0,z1,z2,z3,...,zN
  // where En is the energy deposited by the gamma in the n-th interaction, xn the x position of this interaction, and so on.
  // Just to complicate life (not really, in fact to make it simpler) the std::vectors are already sorted per crystal here, so
  // connecting them to the TBranches is easier. Remember that in the output from the simulation, the vectors above are TTree entries called
  // cryN
  // cryNPosXEnDep
  // cryNPosYEnDep
  // cryNPosZEnDep
  // here we need to store these values in some variables once we read them while looping on the TTree, and we don't know
  // a priori the number of crystals. So instead of a simple pEdep std::vector, we need an array of pEdep std::vectors, with
  // dynamic length (and same holds for px, py and pz). Hence, We don't declare a std::vector, but a pointer to std::vector, which
  // allows to make an array of std::vectors run time.
  // Finally, ROOT comes in the way and forces us not to pass a std::vector to the SetBranchAddress function, but a f...ing **ptr because
  // of template classes. All very complicated but in the end this is how it work:
  // in the event loop, after GetEvent(i) the variable pEdep[i] is a std::vector for the i-th crystal, that has N elements, one per each
  // energy deposition occurred in that crystal.
  // std::vector<float> **pEdep;         // for each energy deposition event, the amount of energy deposited
  // std::vector<float> **px;            // for each energy deposition event, the x position
  // std::vector<float> **py;            // for each energy deposition event, the y position
  // std::vector<float> **pz;            // for each energy deposition event, the z position
  // // create the arrays
  // pEdep = new std::vector<float>* [numOfCry];
  // px    = new std::vector<float>* [numOfCry];
  // py    = new std::vector<float>* [numOfCry];
  // pz    = new std::vector<float>* [numOfCry];
  // // inizialize to 0... or you'll have bad surprises
  // for (int i = 0 ; i < numOfCry ; i++)
  // {
  //   pEdep[i] = 0;
  //   px[i] = 0;
  //   py[i] = 0;
  //   pz[i] = 0;
  // }

  //NEW energy deposition, each gamma 511 event has a std::vector of struct (type enDep) with all the data of each energy deposition
  std::vector<enDep> *energyDeposition = 0;

  //Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];

  //hits on detectors events
  //these are 1 std::vector<float> per TTree entry
  //For each gamma shot, a great number of opticals is generated in the matrix
  //some of them enter the detectors, pass the PDE check and are "detected" (i.e. saved). So
  //each one of these std::vector has length M, where M is the number of optical photons
  //saved for that gamma event. Notice that this M is DIFFERENT from N of the energy deposition part.
  //data is then arranged as above, so for example
  // pOpticalX[M]               = x0,x1,x2,x3,...,xN
  // pOpticalPreMomentumX[M]    = m0,m1,m2,m3,...,mN
  // pOpticalPostMomentumX[M]   = v0,v1,v2,v3,...,vN
  // pGlobalTime[M]             = t0,t1,t2,t3,...,tN
  // and so on.

  //NEW optical photons. for each gamma 511 event, every optical photon detected is a struct of type optPhot. a std::vector<optPhot> is saved for each gamma 511
  std::vector<optPhot> *photons = 0;


  //Position of optical photon when entering detector
  // std::vector<float> *pOpticalX;            // x position of optical photon when entering detector
  // std::vector<float> *pOpticalY;            // y position of optical photon when entering detector
  // std::vector<float> *pOpticalZ;            // z position of optical photon when entering detector
  // //inizialize...
  // pOpticalX = 0;
  // pOpticalY = 0;
  // pOpticalZ = 0;
  // //Momentum before and after entering the detector
  // //before
  // std::vector<float> *pOpticalPreMomentumX;            // x component of momentum unitary vector of optical photon before entering detector
  // std::vector<float> *pOpticalPreMomentumY;            // y component of momentum unitary vector of optical photon before entering detector
  // std::vector<float> *pOpticalPreMomentumZ;            // z component of momentum unitary vector of optical photon before entering detector
  // //inizialize...
  // pOpticalPreMomentumX = 0;
  // pOpticalPreMomentumY = 0;
  // pOpticalPreMomentumZ = 0;
  // //after
  // std::vector<float> *pOpticalPostMomentumX;            // x component of momentum unitary vector of optical photon after entering detector
  // std::vector<float> *pOpticalPostMomentumY;            // y component of momentum unitary vector of optical photon after entering detector
  // std::vector<float> *pOpticalPostMomentumZ;            // z component of momentum unitary vector of optical photon after entering detector
  // //inizialize...
  // pOpticalPostMomentumX = 0;
  // pOpticalPostMomentumY = 0;
  // pOpticalPostMomentumZ = 0;
  // //Global time (i.e. from emission of primary gamma for this event) of detection for the optical
  // std::vector<float> *pGlobalTime;
  // pGlobalTime = 0;//inizialize...
  // //Type of optical photon (0 = scintillation, 1 = Cherenkov, 2 = other (should not exist))
  // std::vector<int> *pPhotonType;
  // pPhotonType = 0;//inizialize...
  // //Energy of optical photon (in eV)
  // std::vector<float> *pPhotonEnergy;
  // pPhotonEnergy = 0;//inizialize...


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
  // for (int i = 0 ; i < numOfCry ; i++)
  // {
  //   std::stringstream snames;
  //   snames << "cry" << i;
  //   tree->SetBranchAddress(snames.str().c_str(),&pEdep[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosXEnDep";
  //   tree->SetBranchAddress(snames.str().c_str(),&px[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosYEnDep";
  //   tree->SetBranchAddress(snames.str().c_str(),&py[i]);
  //   snames.str("");
  //   snames<< "cry" << i << "PosZEnDep";
  //   tree->SetBranchAddress(snames.str().c_str(),&pz[i]);
  // }
  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  // tree->SetBranchAddress("PositionX",&pOpticalX);
  // tree->SetBranchAddress("PositionY",&pOpticalY);
  // tree->SetBranchAddress("PositionZ",&pOpticalZ);
  //
  // tree->SetBranchAddress("PreMomentumX",&pOpticalPreMomentumX);
  // tree->SetBranchAddress("PreMomentumY",&pOpticalPreMomentumY);
  // tree->SetBranchAddress("PreMomentumZ",&pOpticalPreMomentumZ);
  // tree->SetBranchAddress("PostMomentumX",&pOpticalPostMomentumX);
  // tree->SetBranchAddress("PostMomentumY",&pOpticalPostMomentumY);
  // tree->SetBranchAddress("PostMomentumZ",&pOpticalPostMomentumZ);
  //
  // tree->SetBranchAddress("GlobalTime",&pGlobalTime);
  // tree->SetBranchAddress("PhotonType",&pPhotonType);
  // tree->SetBranchAddress("PhotonEnergy",&pPhotonEnergy);



  //output ttree
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
  Double_t xmppc[16]={-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};
  Double_t ymppc[16]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};

  // wrong geometry
  //   Double_t xmppc[32]={-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0,-4.65,0,-1.55,0,1.55,0,4.65,0};
  //   Double_t ymppc[32]={-4.65,0,-4.65,0,-4.65,0,-4.65,0,-1.55,0,-1.55,0,-1.55,0,-1.55,0,1.55,0,1.55,0,1.55,0,1.55,0,4.65,0,4.65,0,4.65,0,4.65,0};


  TH2F *flood = new TH2F("FloodHisto","FloodHisto",1000,-7,7,1000,-7,7);
  flood->GetXaxis()->SetTitle("X [mm]");
  flood->GetYaxis()->SetTitle("Y [mm]");
  flood->GetZaxis()->SetTitle("N");
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
  double columsum = 0;
  double rowsum = 0;
  double total = 0;
  double floodx = 0;
  double floody = 0;
  double floodz = 0;


  // TBranch *b_floodx = tree->Branch("FloodX",&floodx,"floodx/F");
  // TBranch *b_floody = tree->Branch("FloodY",&floody,"floody/F");
  // TBranch *b_floodz = tree->Branch("FloodZ",&floodz,"floodz/F");

  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;

  TH2F *averageZvsRatio = new TH2F("AverageZ vs Ratio","AverageZ vs Ratio",100,0,1,100,-8,8);
  averageZvsRatio->GetXaxis()->SetTitle("Ratio (Fi/(Fi+Bi))");
  averageZvsRatio->GetYaxis()->SetTitle("AverageZ [mm]");

  TH2F *wivsRatio = new TH2F("wi vs Ratio","wi vs Ratio",100,0,1,100,0,1);
  wivsRatio->GetXaxis()->SetTitle("Ratio (Fi/(Fi+Bi))");
  wivsRatio->GetYaxis()->SetTitle("Wi");

  TH2F *wMeasuredvsRatioW = new TH2F("wMeasuredvsRatioW","wMeasuredvsRatioW",100,0,1,100,0,1);
  wMeasuredvsRatioW->GetXaxis()->SetTitle("RatioW Sum_i(Wi*Ei/Etot)");
  wMeasuredvsRatioW->GetYaxis()->SetTitle("Measured W");

  TH2F *averageZvsWi = new TH2F("AverageZ vs Wi","AverageZ vs Wi",100,0,1,100,-8,8);
  averageZvsWi->GetXaxis()->SetTitle("Wi");
  averageZvsWi->GetYaxis()->SetTitle("AverageZ [mm]");

  // for(int i =0 ; i < 2; i++)
  // {
  //   std::stringstream name;
  //   name << "AverageZ vs. (Fi/(Fi+Bi)) - Crystal " << i
  //   averageZvsRatio[i] = new TH2F()
  // }
  long int foundCandidate = 0;
  long int singleCounter = 0;
  long int doubleCounter = 0;
  long int tripleCounter = 0;
  long int multipleCounter = 0;
  for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  {
    tree->GetEvent(iEvent);
    std::vector<int> crystals;
    // std::vector<enDep> EventDepositions;
    for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
    {
      //read the crystal where energy was deposited
      int cry = energyDeposition->at(eEvent).CrystalID;
      //loop in the crystals found
      //look for the same id
      bool sameID = false;
      for(int j = 0 ; j < crystals.size(); j++)
      {
        if(crystals[j] == cry) sameID = true;
      }
      if(!sameID) crystals.push_back(cry);
    }

    //now take only events where energy was deposited in 2 crystals and total energy deposited is 511KeV
    if((crystals[0] == 27 && crystals[1] == 28) | (crystals[0] == 28 && crystals[1] == 27))
    // if(true)
    {
      if(totalEnergyDeposited > 0.510)
      {
        if(crystals.size() == 2)
        {
          //consider only non-lateral crystals
          bool isCandidate = true;
          for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and check is not on borders
          {
            bool check_b = energyDeposition->at(eEvent).CrystalI < 2 | energyDeposition->at(eEvent).CrystalI > 5 | energyDeposition->at(eEvent).CrystalJ < 2 | energyDeposition->at(eEvent).CrystalJ > 5;
            if(check_b)
            {
              isCandidate = false;
            }
          }
          if(isCandidate)
          {
            //calculate average Z of production in the 2 crystals
            // foundCandidate++;

            //calculate the measured W
            int Measured_pmaxID = 0;
            int Measured_pmax = 0;
            int Measured_pSecond = 0;
            int Measured_totalDetCounts = 0;
            std::vector<int> vToSort;
            for(int i = 0; i < numOfCh ; i++)
            {
              Measured_totalDetCounts+=detector[i];
              vToSort.push_back(detector[i]);
            }
            std::sort (vToSort.begin(),vToSort.end());
            Measured_pmax = vToSort[vToSort.size()-1];
            Measured_pSecond = vToSort[vToSort.size()-2];

            // {
            //   Measured_totalDetCounts += detector[i];
            //   if(detector[i] > Measured_pmax1)
            //   {
            //     Measured_pmax1 = detector[i];
            //     Measured_pmaxID = i;
            //   }
            // }



            // double Measured_w = ((double) Measured_pmax) / ((double)Measured_totalDetCounts);
            double Measured_w = ((double) Measured_pmax + (double)Measured_pSecond) / (2.0*(double)Measured_totalDetCounts);
            double RatioW = 0.0;


            std::vector<double> totalEnergyPerCrystal;
            std::vector<double> averageZ;
            std::vector<double> wi;
            std::vector<double> ratio;

            for(int j = 0 ; j < crystals.size() ; j++) // we are in a specific crystal
            {
              double temp_totalEnergyPerCrystal = 0.0;
              for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++)
              {

                if(energyDeposition->at(eEvent).CrystalID == crystals[j])
                temp_totalEnergyPerCrystal += energyDeposition->at(eEvent).EnergyDeposited;
              }

              double temp_averageZ = 0.0;
              for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++)
              {
                if(energyDeposition->at(eEvent).CrystalID == crystals[j])
                temp_averageZ += (energyDeposition->at(eEvent).DepositionZ *  energyDeposition->at(eEvent).EnergyDeposited)/temp_totalEnergyPerCrystal;
              }



              // now on opticals
              // 1.is averageZ correlated to Fi/(Fi+Bi)?
              // count Fi and Bi for this crystal

              // 2. next, is wi correlated with ratio (actually are they the same?)?
              // so calculate wi
              // divide opticals on the basis of origin crystal
              // count the photons detected by each detector for both origins
              // compute the two wi
              long int Fi = 0;
              long int Bi = 0;
              short detectorSorted[numOfCh];
              for(int iDet = 0; iDet < numOfCh ;iDet++)
              {
                detectorSorted[iDet] = 0;
              }

              for(int oEvent = 0; oEvent < photons->size(); oEvent++) //run on energy depositions and find in how many crystals energy was deposited
              {
                // calc origin crystal ID
                int OriginCrystalID = photons->at(oEvent).OriginCrystalI * 8 +  photons->at(oEvent).OriginCrystalJ;
                if(OriginCrystalID == crystals[j])
                {
                  if(photons->at(oEvent).ExitFace == 1)
                  Fi++;
                  if(photons->at(oEvent).ExitFace == 2)
                  Bi++;
                  //the array of detectors
                  // std::cout << photons->at(oEvent).PositionX << " "<< (int) ((photons->at(oEvent).PositionX + 6.3) /  3.2) << std::endl;
                  int iDetector = (int) ((photons->at(oEvent).PositionX + 6.3) /  3.2);
                  int jDetector = (int) ((photons->at(oEvent).PositionY + 6.3) /  3.2);
                  int detectorID = iDetector * 4 +  jDetector;
                  detectorSorted[detectorID]++;
                }
              }

              int pmaxID = 0;
              int pmax = 0;
              int totalDetCounts = 0;
              for(int i = 0; i < numOfCh ; i++)
              {
                totalDetCounts += detectorSorted[i];
                if(detectorSorted[i] > pmax)
                {
                  pmax = detectorSorted[i];
                  pmaxID = i;
                }
              }
              double temp_wi = ((double) pmax) / ((double)totalDetCounts);


              double temp_ratio = ((double) Fi)/((double) Fi+ (double) Bi);
              wi.push_back(temp_wi);
              totalEnergyPerCrystal.push_back(temp_totalEnergyPerCrystal);
              ratio.push_back(temp_ratio);
              averageZ.push_back(temp_averageZ);

              // RatioW += (wi * totalEnergyPerCrystal) / totalEnergyDeposited;




            }

            bool doubleAccepted = true;
            for(int countCry = 0; countCry < wi.size() ;countCry++ )
            {
              if(totalEnergyPerCrystal[countCry] < 0.05) doubleAccepted = false;
            }
            if(doubleAccepted)
            {
              foundCandidate++;
              for(int countCry = 0; countCry < wi.size() ;countCry++ )
              {
                averageZvsRatio->Fill(ratio[countCry],averageZ[countCry]);
                wivsRatio->Fill(ratio[countCry],wi[countCry]);
                averageZvsWi->Fill(wi[countCry],averageZ[countCry]);
                RatioW += (wi[countCry] * totalEnergyPerCrystal[countCry]) / totalEnergyDeposited;
              }
              wMeasuredvsRatioW->Fill(RatioW,Measured_w);
            }

          }
        }
      }
    }







    if(totalEnergyDeposited > 0.510)
    {
      if(crystals.size() == 1) singleCounter++;
      if(crystals.size() == 2) doubleCounter++;
      if(crystals.size() == 3) tripleCounter++;
      if(crystals.size() > 3)  multipleCounter++;
    }
    // std::cout << "Event " << iEvent << " - Dep in crystals = " << crystals.size() << std::endl;

    // columsum = 0;
    // rowsum = 0;
    // total = 0;
    // floodx = 0;
    // floody = 0;
    // floodz = 0;
    //
    //
    //
    // //first clean the array
    // for(int i = 0; i < 16 ; i++)
    // {
    //   charge[i] = 0;
    // }
    //
    // //then fill it with the detector data
    // for(int i = 0; i < numOfCh ; i++)
    // {
    //   charge[i] = detector[i];
    // }
    //
    // //calculate weighted energy and fill 2d histo
    // for(int i = 0; i < 16 ; i++)
    // {
    //   columsum += charge[i]*ymppc[i];
    //   rowsum += charge[i]*xmppc[i];
    //   total += charge[i];
    // }
    // floodx = rowsum/total;
    // floody = columsum/total;
    // b_floodx->Fill();
    //
    // flood->Fill(floodx,floody);


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

  std::cout << std::endl;
  std::cout << "1 cry [511 KeV deposition] events = "    << singleCounter << std::endl;
  std::cout << "2 cry [511 KeV deposition] events = "    << doubleCounter << std::endl;
  std::cout << "3 cry [511 KeV deposition] events = "    << tripleCounter << std::endl;
  std::cout << "Multi cry [511 KeV deposition] events = "<< multipleCounter << std::endl;
  std::cout << "Candidates = "<< foundCandidate << std::endl;

  std::string outFile = "FileOut.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  averageZvsRatio->Write();
  averageZvsWi->Write();
  wivsRatio->Write();
  wMeasuredvsRatioW->Write();
  // t1->Write();

  //   f1->Close();

  // std::string outFile = "analysis_OUT.root";
  // TFile* fOut = new TFile(outFile.c_str(),"recreate");

  // flood->Write();
  // positions->Write();
  //   f1->Close();

  fOut->Close();



  return 0;
}
