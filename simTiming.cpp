// compile with
// g++ -o ../build/simTiming simTiming.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/
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
#include <getopt.h>

#include "TPaveLabel.h"
#include "TPaveText.h"
#include "TLegend.h"
#include "THStack.h"
#include "TSystemDirectory.h"

#include "../code/struct.hh"

// bool compareByTime(const enDep &a,const enDep  &b)
// {
//     return a.DepositionTime < b.DepositionTime;
// }


struct Slice_t
{
  TH1F* pulse;
  double zMin;
  double zMax;
  long int events;
};

struct sipm_t
{
  float detx;
  float dety;
  short int **spad;
  UShort_t counts;
  std::vector<Float_t> listOfTimestamps;
  Float_t timestamp;
  std::vector<Slice_t> slice;
};



bool compare_by_GlobalTime(const optPhot a, const optPhot b)
{
  return a.GlobalTime < b.GlobalTime;
}

bool compare_by_DepositionTime(const enDep a, const enDep b)
{
  return a.DepositionTime < b.DepositionTime;
}


void usage()
{
  std::cout << "\t\t" << "[ -i <input file prefix>    prefix of input files name] " << std::endl
            << "\t\t" << "[ -o <output file>    name of output file] " << std::endl
            << "\t\t" << "[ --photons <N>       average time on first N photons - default = 5]" << std::endl
            << "\t\t" << "[ --saturation        flag to use saturation of mppc - default false ] " << std::endl
            << "\t\t" << "[ --cherenkov         flag to include cherenkov - default false ] " << std::endl
            << "\t\t" << "[ --nmppcx <N>        number of mppc in x             - default = 4]" << std::endl
            << "\t\t" << "[ --nmppcy <N>        number of mppc in y             - default = 4]" << std::endl
            << "\t\t" << "[ --pitchx <N>        distance between center of mppcs , in x [mm] - default = 3.2]" << std::endl
            << "\t\t" << "[ --pitchy <N>        distance between center of mppcs , in y [mm] - default = 3.2]" << std::endl
            << "\t\t" << "[ --qe <N>            quantum efficiency of mppcs  - default = 0.3]" << std::endl
            << "\t\t" << "[ --sptr <N>          sigma sptr of the detector   - default = 0.087]" << std::endl
            << "\t\t" << "[ --doiSlices <N>     number of slices in doi - default = 10]" << std::endl
            << "\t\t" << "[ --length <N>        crystal length in mm - default = 15]" << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{

  if(argc < 2) // check input from command line
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
    usage();
    return 1;
  }



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
  bool cherenkov = false;
  int numb_of_phot_for_time_average = 5;
  std::string inputFileName = "";
  std::string outputFileName = "";
  bool inputGiven = false;
  bool outputGiven = false;
  TChain *tree =  new TChain("tree"); // read input files
  // int nmodulex = 1;
  // int nmoduley = 1;
  int nmppcx = 4;
  int nmppcy = 4;
  float pitchx = 3.2;
  float pitchy = 3.2;
  double qe = 0.3;
  double sigmaSPTR = 0.087;
  int doiSlices = 10;
  double length = 15.0;
  std::string filePrefix = "";
  // int ncrystalsx = 2;
  // int ncrystalsy = 2;

  static struct option longOptions[] =
  {
      { "photons", required_argument, 0, 0 },
      { "saturation", no_argument, 0, 0 },
      { "nmppcx", required_argument, 0, 0 },
      { "nmppcy", required_argument, 0, 0 },
      { "pitchx", required_argument, 0, 0 },
      { "pitchy", required_argument, 0, 0 },
      { "qe", required_argument, 0, 0 },
      { "sptr", required_argument, 0, 0 },
      { "doiSlices", required_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "cherenkov", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			filePrefix = (char *)optarg;
      // std::cout << "Adding file " << inputFileName << std::endl;
      // tree->Add(inputFileName.c_str());
      inputGiven = true;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
      outputGiven = true;
    }
		else if (c == 0 && optionIndex == 0){
      numb_of_phot_for_time_average = atoi((char *)optarg);
      std::cout << "Time average on first " << numb_of_phot_for_time_average << " photons"<< std::endl;
    }
    else if (c == 0 && optionIndex == 1){
      std::cout << "SiPM saturation will be taken into account " << std::endl;
      saturation = true;
    }
    else if (c == 0 && optionIndex == 2){
      nmppcx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3){
      nmppcy = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      pitchx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      pitchx = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      qe = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      sigmaSPTR = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      doiSlices = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      length = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      std::cout << "Including Cherenkov photons " << std::endl;
      cherenkov = true;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(!inputGiven | !outputGiven)
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
		usage();
		return 1;
  }


  //open all files that begin with prefix, add them to TChain
  TSystemDirectory dir("./","./" );
  TList *files = dir.GetListOfFiles();
  int nfiles = files->GetEntries();
  std::vector<std::string> filesName;
  // fill a vector with the file names
  for(int i = 0 ; i < nfiles ; i++){
    filesName.push_back(files->At(i)->GetName());
  }

  for(unsigned int i = 0 ; i < filesName.size() ; i++)
  {
    if (!filesName[i].compare(0, filePrefix.size(), filePrefix))
    {
      std::cout << "Adding file " << filesName[i] << std::endl;
      tree->Add(filesName[i].c_str());
    }
  }

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

  // std::string outFileName = "treeout.root"; //+ std::string(argv[1]);
  //output ttree
  // long long int DeltaTimeTag,ExtendedTimeTag;
  //
  // UShort_t *charge; //adc type
  // charge = new UShort_t[numOfCh];
  // Float_t *timestamp;
  // timestamp = new Float_t[numOfCh];
  // UShort_t taggingCharge = 1;
  // Float_t taggingTimeStamp = 0;
  // Float_t RealX,RealY,RealZ;
  // Short_t CrystalsHit;
  // Short_t NumbOfInteractions;
  // Float_t TotalEnergyDeposited_out;
  //
  // TTree* t1 = new TTree("adc","adc");
  //
  // t1->Branch("ExtendedTimeTag",&ExtendedTimeTag,"ExtendedTimeTag/l"); 	//absolute time tag of the event
  // t1->Branch("DeltaTimeTag",&DeltaTimeTag,"DeltaTimeTag/l"); 			//delta time from previous event
  // //branches of the channels data
  // std::stringstream snames,stypes;
  // for (int i = 0 ; i < numOfCh ; i++)
  // {
  //   //empty the stringstreams
  //
  //   charge[i] = 0;
  //   timestamp[i] = 0;
  //   snames << "ch" << i;
  //   stypes << "ch" << i << "/s";
  //   t1->Branch(snames.str().c_str(),&charge[i],stypes.str().c_str());
  //   snames.str("");
  //   stypes.str("");
  //   snames << "t" << i;
  //   stypes << "t" << i << "/F";
  //   t1->Branch(snames.str().c_str(),&timestamp[i],stypes.str().c_str());
  //   snames.str("");
  //   stypes.str("");
  // }
  // //create a fake additional channel, faking an external tagging crystal, for ModuleCalibration
  // snames << "ch" << numOfCh;
  // stypes << "ch" << numOfCh << "/s";
  // t1->Branch(snames.str().c_str(),&taggingCharge,stypes.str().c_str());
  // snames.str("");
  // stypes.str("");
  // snames << "t" << numOfCh;
  // stypes << "t" << numOfCh << "/F";
  // t1->Branch(snames.str().c_str(),&taggingTimeStamp,stypes.str().c_str());
  // snames.str("");
  // stypes.str("");
  // t1->Branch("RealX",&RealX,"RealX/F");
  // t1->Branch("RealY",&RealY,"RealY/F");
  // t1->Branch("RealZ",&RealZ,"RealZ/F");
  // t1->Branch("CrystalsHit",&CrystalsHit,"CrystalsHit/S");
  // t1->Branch("NumbOfInteractions",&NumbOfInteractions,"NumbOfInteractions/S");
  // t1->Branch("TotalEnergyDeposited",&TotalEnergyDeposited_out,"TotalEnergyDeposited/F");

  //saturation part
  //create an array of spads
  int n_spad_x = 60;
  int n_spad_y = 60;
  int n_dead_spad_x = 4;
  int n_dead_spad_y = 4;


  float *xmppc;
  float *ymppc;
  xmppc = new float[numOfCh];
  ymppc = new float[numOfCh];

  for(int i = 0; i < nmppcx; i++)
  {
    for(int j = 0 ; j < nmppcy;j++)
    {
      xmppc[i*nmppcy+j] = (i * pitchx) - (pitchx*nmppcx/2.0) + (pitchx/2.0);
      ymppc[i*nmppcy+j] = (j * pitchy) - (pitchy*nmppcy/2.0) + (pitchy/2.0);
    }
  }
  // float xmppc[16] = {-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  // float ymppc[16] = {-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};

  // std::cout << "xmppc = {";
  // for(int i = 0 ; i < numOfCh ; i++)
  // {
  //   std::cout << xmppc[i]<< ",";
  // }
  // std::cout << "}" << std::endl;
  //
  // std::cout << "ymppc = {";
  // for(int i = 0 ; i < numOfCh ; i++)
  // {
  //   std::cout << ymppc[i]<< ",";
  // }
  // std::cout << "}"<< std::endl;

  // float xmppc[64] = {-11.2,-11.2,-11.2,-11.2,-11.2,-11.2,-11.2,-11.2,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-8.0,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8,4.8,4.8,4.8,4.8,8.0,8.0,8.0,8.0,8.0,8.0,8.0,8.0,11.2,11.2,11.2,11.2,11.2,11.2,11.2,11.2};
  // float ymppc[64] = {-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2,-11.2,-8.0,-4.8,-1.6,1.6,4.8,8.0,11.2};
  // float *xmppc;
  // float *ymppc;
  // xmppc = new float[8*8];
  // ymppc = new float[8*8];
  // for(int iMppc = 0; iMppc < 8 ; iMppc++)
  // {
  //   for(int dd = 0; dd < 8 ; dd++)
  //   {
  //     xmppc[iMppc+dd] =  iMppc * 3.2;
  //   }
  // }
  // double detector_pitch = 3.2;
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

    // std::vector<Slice_t> slice;

    for(int iDoi = 0; iDoi < doiSlices ; iDoi++)
    {
      Slice_t tempSlice;
      tempSlice.zMin = -(length/2.0) + iDoi*(length/doiSlices);
      tempSlice.zMax = -(length/2.0) + (iDoi+1)*(length/doiSlices);
      std::stringstream Hname;
      Hname << "Ch "<< i << " Pulse z " << tempSlice.zMin << " to " << tempSlice.zMax ;
      tempSlice.pulse = new TH1F(Hname.str().c_str(),Hname.str().c_str(),1000,0,2.5);
      tempSlice.events = 0;
      sipm[i].slice.push_back(tempSlice);
    }
  }
  // for(int iSipm = 0; iSipm < numOfCh ; iSipm++)
  // {
  //   std::cout << iSipm << " "
  //             << sipm[iSipm].detx - det_size_x/2.0 << " "
  //             << sipm[iSipm].detx + det_size_x/2.0 << " "
  //             << sipm[iSipm].dety - det_size_y/2.0 << " "
  //             << sipm[iSipm].dety + det_size_y/2.0 << " "
  //             << std::endl;
  // }

  //Prepare histograms of "pulses"
  // slice the doi coordinate doiSlices times, for each build and average "pulse"




  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;

  for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  {
    tree->GetEvent(iEvent);

    // ExtendedTimeTag = 1e-9;
    // DeltaTimeTag = 1e-9;

    int NumbOfInteractions = 0;
    int CrystalsHit = 0;

    // sort the optical photons in time
    std::sort(photons->begin(), photons->end(), compare_by_GlobalTime );
    //also sort the energy depositions in time
    std::sort(energyDeposition->begin(), energyDeposition->end(), compare_by_DepositionTime );







    //------------------------//
    // ENERGY DEPOSITIONS     //
    //------------------------//
    double RealX = 0;
    double RealY = 0;
    double RealZ = 0;
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
      RealX += (energyDeposition->at(eEvent).DepositionX * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
      RealY += (energyDeposition->at(eEvent).DepositionY * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
      RealZ += (energyDeposition->at(eEvent).DepositionZ * energyDeposition->at(eEvent).EnergyDeposited)/totalEnergyDeposited;
    }


    // TotalEnergyDeposited_out = totalEnergyDeposited;
    CrystalsHit = crystals.size();



    // if(NumbOfInteractions > 0) // discard events with no energy deposition (they would never trigger the detectors anyway..)
    // {
    //   t1->Fill();
    // }

    if(totalEnergyDeposited > 0.5 && CrystalsHit == 1 && crystals[0] == 36)
    {
      for(int iSipm = 0; iSipm < numOfCh ; iSipm++)
      {
        for(unsigned int iSlice = 0 ; iSlice < sipm[iSipm].slice.size(); iSlice++)
        {
          if(RealZ > sipm[iSipm].slice[iSlice].zMin && RealZ < sipm[iSipm].slice[iSlice].zMax)
          {
            // sipm[iSipm].slice[iSlice].pulse->Fill(photons->at(iPhot).GlobalTime);
            sipm[iSipm].slice[iSlice].events++;
          }
        }
      }
    }
    //------------------------//
    // OPTICAL PHOTONS        //
    //------------------------//
    for(int iPhot = 0; iPhot < photons->size(); iPhot++) // run on all opticals
    {
      //skip cherenkov unless they are included
      if(photons->at(iPhot).PhotonType != 1 || cherenkov)
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
                    if(sipm[iSipm].spad[hiti][hitj] == 0) // if this spad was not hit yet
                    {
                      sipm[iSipm].counts++;
                      sipm[iSipm].spad[hiti][hitj] = 1;
                      sipm[iSipm].listOfTimestamps.push_back((Float_t) photons->at(iPhot).GlobalTime); //add its time stamp to the sipm
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
                if(numb < qe)  // accepted photon
                {
                  sipm[iSipm].counts++;
                  sipm[iSipm].listOfTimestamps.push_back((Float_t) photons->at(iPhot).GlobalTime); //add its time stamp to the sipm
                  if(totalEnergyDeposited > 0.5 && CrystalsHit == 1 && crystals[0] == 36)
                  {
                    for(unsigned int iSlice = 0 ; iSlice < sipm[iSipm].slice.size(); iSlice++)
                    {
                      if(RealZ > sipm[iSipm].slice[iSlice].zMin && RealZ < sipm[iSipm].slice[iSlice].zMax)
                      {
                        sipm[iSipm].slice[iSlice].pulse->Fill(photons->at(iPhot).GlobalTime);
                        // sipm[iSipm].slice[iSlice].events++;
                      }
                    }
                  }
                }
              }

            }
          }
        }
      }

    }
    // calculate the global sipm parameters
    for(int i = 0; i < numOfCh ; i++)
    {
      // fill the charge vector
      // charge[i] = (UShort_t) sipm[i].counts;
      // calculate the sipm timestamp from average of first N timestamps
      sipm[i].timestamp = 0.0;
      int effectiveN = numb_of_phot_for_time_average;
      if(numb_of_phot_for_time_average > sipm[i].listOfTimestamps.size())
        effectiveN = sipm[i].listOfTimestamps.size();
      for(int j = 0 ; j < effectiveN; j++)
      {
        sipm[i].timestamp +=  (Float_t) ((gRandom->Gaus(sipm[i].listOfTimestamps[j],sigmaSPTR) / effectiveN)*1e-9); // default smearing at 0.087, and convert to seconds
      }
      // timestamp[i] = (Float_t) sipm[i].timestamp;
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
      //clear time stamps;
      sipm[i].listOfTimestamps.clear();
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
  std::cout << "Writing output to file "<< outputFileName << std::endl;



  TFile* fOut = new TFile(outputFileName.c_str(),"recreate");
  // t1->Write();
  for(int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream sStack;
    sStack << "Pulses of Ch " << i;
    THStack *hs = new THStack(sStack.str().c_str(),sStack.str().c_str());
    TCanvas *canvas = new TCanvas(sStack.str().c_str(),sStack.str().c_str(),1200,800);

    sStack.str("");
    sStack << "Normalized Pulses of Ch " << i;
    THStack *hsNorm = new THStack(sStack.str().c_str(),sStack.str().c_str());
    TCanvas *canvasNorm = new TCanvas(sStack.str().c_str(),sStack.str().c_str(),1200,800);

    TLegend *legend = new TLegend(0.15,0.62,0.30,0.89,"");
    legend->SetFillStyle(0);

    TLegend *legendNorm = new TLegend(0.15,0.62,0.30,0.89,"");
    legendNorm->SetFillStyle(0);
    // sname.str("");
    // sname << "No correction        = " << round(1e12*crystal[iCry].simpleCTR->GetRMS()*2.355*TMath::Sqrt(2.0)) << "ps";
    int color = 1; // avoid white
    for(unsigned int iSlice = 0 ; iSlice < sipm[i].slice.size(); iSlice++)
    {
      std::stringstream sDoi;
      sDoi << "z = " <<(sipm[i].slice[iSlice].zMax + sipm[i].slice[iSlice].zMin) /2.0 << "mm";
      if(color == 5)
      {
        color++; // no yellow!!
      }
      sipm[i].slice[iSlice].pulse->SetLineColor(color);
      color++;
      //normalize the histograms to the number of events in the doi slice - mandatory to average
      if(sipm[i].slice[iSlice].events != 0)
      {
        sipm[i].slice[iSlice].pulse->Scale(1.0/sipm[i].slice[iSlice].events);
      }
      //clone the histo
      TH1F *clone = (TH1F*) sipm[i].slice[iSlice].pulse->Clone();
      hs->Add(clone);
      legend->AddEntry(clone,sDoi.str().c_str(),"l");

      //normalize histograms to max
      int maxBin = sipm[i].slice[iSlice].pulse->GetMaximumBin();
      double maxValue = sipm[i].slice[iSlice].pulse->GetBinContent(maxBin);
      sipm[i].slice[iSlice].pulse->Scale(1.0/maxValue);
      hsNorm->Add(sipm[i].slice[iSlice].pulse);
      legendNorm->AddEntry(sipm[i].slice[iSlice].pulse,sDoi.str().c_str(),"l");
      sipm[i].slice[iSlice].pulse->Write();
    }
    canvas->cd();
    hs->Draw("nostack");
    legend->Draw();
    canvas->Write();


    canvasNorm->cd();
    hsNorm->Draw("nostack");
    legendNorm->Draw();
    canvasNorm->Write();
  }
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
  // delete charge;
  // delete timestamp;
  // delete sipm;
  delete rand;
  return 0;
}
