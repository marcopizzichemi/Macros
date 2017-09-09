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
#include "TRandom3.h"

#include "../code/struct.hh"

// bool compareByTime(const enDep &a,const enDep  &b)
// {
//     return a.DepositionTime < b.DepositionTime;
// }

struct sipm_t
{
  float detx;
  float dety;
  short int **spad;
  int counts;
};


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


  bool saturation = false;

  //-------------------
  // Input Files or Flags
  //-------------------
  TChain *tree =  new TChain("tree"); // read input files
  for (int i = 1 ; i < argc ; i++)
  {
    std::string argument(argv[i]);

    if(argument.compare("--saturation") == 0)
    {
      std::cout << "SiPM saturation will be taken into account " << std::endl;
      saturation = true;
    }
    else
    {
      std::cout << "Adding file " << argv[i] << std::endl;
      tree->Add(argv[i]);
    }
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
  // prepare also a vector for optical photons

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
  //output ttree
  long long int DeltaTimeTag,ExtendedTimeTag;

  Short_t *charge; //adc type
  charge = new Short_t[numOfCh];
  Float_t RealX,RealY,RealZ;
  Short_t CrystalsHit;
  Short_t NumbOfInteractions;
  Float_t TotalEnergyDeposited_out;

  TTree* t1 = new TTree("adc","adc");

  t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
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
  t1->Branch("TotalEnergyDeposited",&TotalEnergyDeposited_out,"TotalEnergyDeposited/F");

  //saturation part
  //create an array of spads
  int n_spad_x = 60;
  int n_spad_y = 60;
  int n_dead_spad_x = 4;
  int n_dead_spad_y = 4;
  double qe = 0.3;

  float xmppc[16] = {-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  float ymppc[16] = {-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};

  double detector_pitch = 3.2;
  double det_size_x = 3.0;
  double det_size_y = 3.0;
  double spad_size_x = det_size_x / n_spad_x;
  double spad_size_y = det_size_y / n_spad_y;

  sipm_t* sipm;
  sipm = new sipm_t[numOfCh];
  TRandom3 *rand = new TRandom3(0);
  //initialize spads
  //they are in the same order of ch0, ch1, etc..
  //FIXME not general at all!!!
  for(int i = 0 ; i < numOfCh ; i++)
  {
    //set position of he sipm
    sipm[i].detx = xmppc[i];
    sipm[i].dety = ymppc[i];
    //create the 2d array of spads
    sipm[i].spad = new short int*[n_spad_x];
    for(int j = 0; j < n_spad_x; j++)
    {
      sipm[i].spad[j] = new short int[n_spad_y];
    }
    //fill the array of spads with 0s
    for(int iSpad = 0; iSpad < n_spad_x; iSpad++)
    {
      for(int jSpad = 0; jSpad < n_spad_y; jSpad++)
      {
        sipm[i].spad[iSpad][jSpad] = 0;
      }
    }
    //set count to 0
    sipm[i].counts = 0;
  }


  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;

  for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  {
    tree->GetEvent(iEvent);

    ExtendedTimeTag = 1e-9;
    DeltaTimeTag = 1e-9;

    NumbOfInteractions = 0;
    CrystalsHit = 0;


    for(int iPhot = 0; iPhot < photons->size(); iPhot++) // run on all opticals
    {
      // find which sipm was hit
      for(int iSipm = 0; iSipm < numOfCh ; iSipm++)
      {
        if((photons->at(iPhot).PositionX > (sipm[iSipm].detx - det_size_x/2.0 )) && (photons->at(iPhot).PositionX < (sipm[iSipm].detx + det_size_x/2.0 )) )
        {
          if((photons->at(iPhot).PositionY > (sipm[iSipm].dety - det_size_y/2.0 )) && (photons->at(iPhot).PositionY < (sipm[iSipm].dety + det_size_y/2.0 )) )
          {
            if(saturation)
            {
              // find which spad was hit
              // bring sipm hit to start in 0,0
              float hitx = photons->at(iPhot).PositionX - sipm[iSipm].detx + (det_size_x/2.0);
              float hity = photons->at(iPhot).PositionY - sipm[iSipm].dety + (det_size_y/2.0);
              // find spad I and J
              int hiti = (int) (hitx / spad_size_x);
              int hitj = (int) (hity / spad_size_y);
              // std::cout << photons->at(iPhot).PositionX << "\t"
              //           << photons->at(iPhot).PositionY << "\t"
              //           << sipm[iSipm].detx << "\t"
              //           << sipm[iSipm].dety << "\t"
              //           << hitx << "\t"
              //           << hity << "\t"
              //           << hiti << "\t"
              //           << hitj << "\t"
              //           << std::endl;
              // ignore the NxN central spads
              // 0-27 (28-29-20-31) 32-59
              if( (hiti > ( n_spad_x/2 - n_dead_spad_x/2 - 1 ) ) &&
              (hiti < ( n_spad_x/2 + n_dead_spad_x/2 - 1 ) ) &&
              (hitj > ( n_spad_y/2 - n_dead_spad_y/2 - 1 ) ) &&
              (hitj < ( n_spad_y/2 + n_dead_spad_y/2 - 1 ) ) )
              {
                // do nothing, this part of the sipm is not active
              }
              else // increment the counts of the sipm, if the spad was not hit yet
              {
                //HACK to avoid seg fault when the optical photon is exactly on the border (which makes the hiti of hitj being exatly 60 for example)
                if(hiti == n_spad_x) hiti = hiti -1;
                if(hitj == n_spad_y) hitj = hitj -1;

                //quantum efficiency test

                double numb = rand->Uniform(1.0);

                if(numb < qe)
                {
                  if(sipm[iSipm].spad[hiti][hitj] == 0)
                  {
                    sipm[iSipm].counts++;
                    sipm[iSipm].spad[hiti][hitj] = 1;
                  }
                  else
                  {
                    //ignore the hit, the spad has already fired and the optical photon is lost
                  }
                }
              }
            }
            else // just qe test for each photon and sipm
            {
              // TRandom3 *rand = new TRandom3(0);
              double numb = rand->Uniform(1.0);
              if(numb < qe)
              {
                sipm[iSipm].counts++;
              }
            }
          }
        }
      }
    }
    // fill the charge vector
    for(int i = 0; i < numOfCh ; i++)
    {
      charge[i] = (Short_t) sipm[i].counts;
    }
    // re-initialize the sipms counters
    for(int i = 0 ; i < numOfCh ; i++)
    {
      //fill the array of spads with 0s
      for(int iSpad = 0; iSpad < n_spad_x; iSpad++)
      {
        for(int jSpad = 0; jSpad < n_spad_y; jSpad++)
        {
          sipm[i].spad[iSpad][jSpad] = 0;
        }
      }
      //set count to 0
      sipm[i].counts = 0;
    }




    RealX = RealY = RealZ = 0;

    NumbOfInteractions = energyDeposition->size();

    std::vector<int> crystals;
    for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++)// run on energy depositions for this gamma event
    {
      // -- counting the crystals where energy was deposited in this event
      //read the crystal where energy was deposited
      int cry = energyDeposition->at(eEvent).CrystalID;
      //loop in the crystals found
      //look for the same id
      bool sameID = false;
      for(int j = 0 ; j < crystals.size(); j++)
      {
        if(crystals[j] == cry) sameID = true;
      }
      if(!sameID) crystals.push_back(cry);  // add the crystal if it was not already counted as hit

      // -- calculate the average coordinate of energy deposition

	      // RealX += (px[i]->at(j) * pEdep[i]->at(j))/totalEnergyDeposited;
      RealX += (energyDeposition->at(eEvent).DepositionX * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
      RealY += (energyDeposition->at(eEvent).DepositionY * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
      RealZ += (energyDeposition->at(eEvent).DepositionZ * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
    }
    TotalEnergyDeposited_out = totalEnergyDeposited;
    CrystalsHit = crystals.size();

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


  }
  std::cout << std::endl;
  std::cout << "Writing output to file "<< outFileName << std::endl;



  TFile* fOut = new TFile(outFileName.c_str(),"recreate");
  t1->Write();

  fOut->Close();

  //free memory
  for(int i = 0 ; i < numOfCh ; i++)
  {
    for(int j = 0; j < n_spad_x; j++)
    {
      delete sipm[i].spad[j];
    }
    delete sipm[i].spad;
  }
  delete detector;
  delete charge;
  delete sipm;
  delete rand;
  return 0;
}
