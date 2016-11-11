// compile with
// g++ -o ../build/simRead_base simRead_base.cpp `root-config --cflags --glibs`
// syntax
// simRead_base `ls out*`

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

int main (int argc, char** argv)
{
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors

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
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;


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
  std::vector<float> **pEdep;         // for each energy deposition event, the amount of energy deposited
  std::vector<float> **px;            // for each energy deposition event, the x position
  std::vector<float> **py;            // for each energy deposition event, the y position
  std::vector<float> **pz;            // for each energy deposition event, the z position
  // create the arrays
  pEdep = new std::vector<float>* [numOfCry];
  px    = new std::vector<float>* [numOfCry];
  py    = new std::vector<float>* [numOfCry];
  pz    = new std::vector<float>* [numOfCry];
  // inizialize to 0... or you'll have bad surprises
  for (int i = 0 ; i < numOfCry ; i++)
  {
    pEdep[i] = 0;
    px[i] = 0;
    py[i] = 0;
    pz[i] = 0;
  }

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

  //Position of optical photon when entering detector
  std::vector<float> *pOpticalX;            // x position of optical photon when entering detector
  std::vector<float> *pOpticalY;            // y position of optical photon when entering detector
  std::vector<float> *pOpticalZ;            // z position of optical photon when entering detector
  //inizialize...
  pOpticalX = 0;
  pOpticalY = 0;
  pOpticalZ = 0;
  //Momentum before and after entering the detector
  //before
  std::vector<float> *pOpticalPreMomentumX;            // x component of momentum unitary vector of optical photon before entering detector
  std::vector<float> *pOpticalPreMomentumY;            // y component of momentum unitary vector of optical photon before entering detector
  std::vector<float> *pOpticalPreMomentumZ;            // z component of momentum unitary vector of optical photon before entering detector
  //inizialize...
  pOpticalPreMomentumX = 0;
  pOpticalPreMomentumY = 0;
  pOpticalPreMomentumZ = 0;
  //after
  std::vector<float> *pOpticalPostMomentumX;            // x component of momentum unitary vector of optical photon after entering detector
  std::vector<float> *pOpticalPostMomentumY;            // y component of momentum unitary vector of optical photon after entering detector
  std::vector<float> *pOpticalPostMomentumZ;            // z component of momentum unitary vector of optical photon after entering detector
  //inizialize...
  pOpticalPostMomentumX = 0;
  pOpticalPostMomentumY = 0;
  pOpticalPostMomentumZ = 0;
  //Global time (i.e. from emission of primary gamma for this event) of detection for the optical
  std::vector<float> *pGlobalTime;
  pGlobalTime = 0;//inizialize...
  //Type of optical photon (0 = scintillation, 1 = Cherenkov, 2 = other (should not exist))
  std::vector<int> *pPhotonType;
  pPhotonType = 0;//inizialize...
  //Energy of optical photon (in eV)
  std::vector<float> *pPhotonEnergy;
  pPhotonEnergy = 0;//inizialize...


  //------------------------
  // Set Branch Addresses
  //------------------------
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
  tree->SetBranchAddress("PositionX",&pOpticalX);
  tree->SetBranchAddress("PositionY",&pOpticalY);
  tree->SetBranchAddress("PositionZ",&pOpticalZ);

  tree->SetBranchAddress("PreMomentumX",&pOpticalPreMomentumX);
  tree->SetBranchAddress("PreMomentumY",&pOpticalPreMomentumY);
  tree->SetBranchAddress("PreMomentumZ",&pOpticalPreMomentumZ);
  tree->SetBranchAddress("PostMomentumX",&pOpticalPostMomentumX);
  tree->SetBranchAddress("PostMomentumY",&pOpticalPostMomentumY);
  tree->SetBranchAddress("PostMomentumZ",&pOpticalPostMomentumZ);

  tree->SetBranchAddress("GlobalTime",&pGlobalTime);
  tree->SetBranchAddress("PhotonType",&pPhotonType);
  tree->SetBranchAddress("PhotonEnergy",&pPhotonEnergy);


  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  for(int iEvent = 0; iEvent < 1 ; iEvent++)
  {

    tree->GetEvent(iEvent);
    int CrystalsHit = 0;      // counter of how many crystals were hit by this gamma

    float* energyPerCrystal;  // total energy deposited in each crystal
    energyPerCrystal = new float[numOfCry];
    for(int i = 0; i < numOfCry ; i++)
    {
      energyPerCrystal[i] = 0;
    }

    // some example of global variable
    std::cout << "Event Numb \t\t\t = "              << iEvent << std::endl;
    std::cout << "Total energy deposited \t\t = "  << totalEnergyDeposited << std::endl;
    std::cout << "Total Numb of Opticals produced\t = "  << NumOptPhotons+NumCherenkovPhotons << std::endl;

    // one example of energy deposition related std::vectors
    for(int i = 0; i < numOfCry ; i++) //crystal by crystal
    {
      if(px[i]->size()) //if the i-th crystal had energy deposition (--> if(0) => false)
      {
        CrystalsHit++; //count the crystal as hit
        for(int j = 0; j < px[i]->size(); j++) //run on the hits in this crystal
        {
          energyPerCrystal[i] += pEdep[i]->at(j); //add the energy deposited
        }
        //output to user
        std::cout << "Crystal[" << i <<"]\t\t\t = "  << energyPerCrystal[i] << std::endl;
      }
    }
    std::cout << "Crystals Hit\t\t\t = " << CrystalsHit << std::endl;

    // one example of optical photon related std::vectors
    // get the total number of optical detected
    // the size of any optical photon related std::vector is M (see above)
    std::cout << "Number of Opticals detected\t = " << pGlobalTime->size() << std::endl;
    for(int i = 0; i < pGlobalTime->size() ; i++) //crystal by crystal
    {
      std::cout << pGlobalTime->at(i) << std::endl;
    }

    std::cout << "-----------------------------" << std::endl;
    std::cout << std::endl;





    //progress feedback to the user...
    counter++;
    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
    }
  }
  std::cout << std::endl;
  TFile *foooooo = new TFile("prova.root","recreate");
  return 0;
}
