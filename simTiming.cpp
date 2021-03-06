// compile with
// g++ -o ../build/simTiming simTiming.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/

//---------------------------------//
//                                 //
// simTiming                       //
//                                 //
//---------------------------------//

// Program to analyze timing in a scintillator matrix
// 1. Input is the output of a g4matrix simulation. The input is passed by specifying the prefix common to all simulation ouput file, via the flag -i
// 2. The crystal hit is split in N parts, and pulses/histos are built for each one
// 3. An impulse is built for each part for each event
//     a. by time of arrival of photons
//     b. with sum of single spad impulses
// 4. Calculates time stamp of each event
//     a. by average of first M photons (M input by user)
//     b. by crossing of fixed threshold (thr input by user)
// 5. Outputs global CTR and CTR for each part 

// syntax explained by running without args:

// Usage: simTiming
	// 	[ -i <input file prefix>    prefix of input files name]
	// 	[ -o <output file>          name of output file]
	// 	[ --photons <N>             average time on first N photons                  - default = 5]
	// 	[ --saturation              flag to use saturation of mppc                   - default = false]
	// 	[ --nmppcx <N>              number of mppc in x direction                    - default = 4]
	// 	[ --nmppcy <N>              number of mppc in y direction                    - default = 4]
	// 	[ --pitchx <N>              distance x between center of mppcs [mm]          - default = 3.2]
	// 	[ --pitchy <N>              distance y between center of mppcs [mm]          - default = 3.2]
	// 	[ --qeScint <N>             quantum eff. for scintillation photons[0-1]      - default = 0.3]
	// 	[ --qeCher <N>              quantum eff. for cherenkov photons [0-1]         - default = 0.2]
	// 	[ --sptr <N>                sigma sptr of the detector [ps]                  - default = 0, no smearing]
	// 	[ --doiSlices <N>           number of slices in doi                          - default = 10]
	// 	[ --length <N>              crystal length in mm                             - default = 15]
	// 	[ --timeBinSize <N>         size of individual time bin [ps]                 - default = 5]
	// 	[ --spadPulse               flag to use read spad pulses                     - default = false]
	// 	[ --spadRise <N>            spad pulse rise time [ps]                        - default = 200]
	// 	[ --spadDecay <N>           spad pulse decay time [ps]                       - default = 10000]
	// 	[ --pulseEnd <N>            length of pluses [ps]                            - default = 10000]


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
#include "TMultiGraph.h"
#include "TGraph.h"

// #include <omp.h>

#include "../code/struct.hh"

// bool compareByTime(const enDep &a,const enDep  &b)
// {
//     return a.DepositionTime < b.DepositionTime;
// }


struct Slice_t
{
  TH1F* pulse;
  TH1F* ctr;
  double zMin;
  double zMax;
  long int events;
  // std::vector<Float_t> listOfTimestamps;
};

struct sipm_t
{
  TH1F* ctr;
  float detx;
  float dety;
  short int **spad;
  UShort_t counts;
  std::vector<Float_t> listOfTimestamps;
  Float_t timestamp;
  std::vector<Slice_t> slice;
};

double singleSpad(double t, double t0, double tr,double td)
{
  double As = 1.0;
  double b = td/tr;
  double Sa = As*0.1;
  double A = (1.0/(2*3.1415*Sa))*exp(-(As)/(2.0*Sa));
  double C = (A)/(pow(b,(1.0/(1.0-b))) - pow(b,(1.0/((1.0/b)-1.0))));
  double s = C * (exp( -(t-t0)/(td) ) - exp( -(t-t0)/(tr) ) );
  if(t>=t0)
  {
    return s;
  }
  else
  {
    return 0;
  }
}



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
            << "\t\t" << "[ -o <output file>          name of output file] " << std::endl
            << "\t\t" << "[ --photons <N>             if spadPulse is not used: average time of first N (integer) photons                               - default = 5]" << std::endl
            << "\t\t" << "[                           if spadPulse is used: trigger threshold, expressed in multiples (double) of single photon pulses  - default = 5]" << std::endl
            << "\t\t" << "[ --saturation              flag to use saturation of mppc                   - default = false] " << std::endl
            // << "\t\t" << "[ --cherenkov         flag to include cherenkov               - default false ] " << std::endl
            << "\t\t" << "[ --nmppcx <N>              number of mppc in x direction                    - default = 4]" << std::endl
            << "\t\t" << "[ --nmppcy <N>              number of mppc in y direction                    - default = 4]" << std::endl
            << "\t\t" << "[ --pitchx <N>              distance x between center of mppcs [mm]          - default = 3.2]" << std::endl
            << "\t\t" << "[ --pitchy <N>              distance y between center of mppcs [mm]          - default = 3.2]" << std::endl
            << "\t\t" << "[ --qeScint <N>             quantum eff. for scintillation photons[0-1]      - default = 0.3]" << std::endl
            << "\t\t" << "[ --qeCher <N>              quantum eff. for cherenkov photons [0-1]         - default = 0.2]" << std::endl
            << "\t\t" << "[ --sptr <N>                sigma sptr of the detector [ps]                  - default = 0, no smearing]" << std::endl
            << "\t\t" << "[ --doiSlices <N>           number of slices in doi                          - default = 10]" << std::endl
            << "\t\t" << "[ --length <N>              crystal length in mm                             - default = 15]" << std::endl
            << "\t\t" << "[ --timeBinSize <N>         size of individual time bin [ps]                 - default = 5]" << std::endl
            << "\t\t" << "[ --spadPulse               flag to use read spad pulses                     - default = false] " << std::endl
            << "\t\t" << "[ --spadRise <N>            spad pulse rise time [ps]                        - default = 200]" << std::endl
            << "\t\t" << "[ --spadDecay <N>           spad pulse decay time [ps]                       - default = 10000]" << std::endl
            << "\t\t" << "[ --pulseEnd <N>            length of pluses [ps]                            - default = 10000]" << std::endl
            << "\t\t" << std::endl;
}

// # ifndef __CINT__
int main (int argc, char** argv)
{

  if(argc < 2) // check input from command line
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
    usage();
    return 1;
  }


  // ROOT::EnableThreadSafety();
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
  double numb_of_phot_for_time_average = 5.0;
  std::string inputFileName = "";
  std::string outputFileName = "";
  bool inputGiven = false;
  bool outputGiven = false;
  bool spadPulse = false;
  TChain *tree =  new TChain("tree"); // read input files
  // int nmodulex = 1;
  // int nmoduley = 1;
  int nmppcx = 4;
  int nmppcy = 4;
  float pitchx = 3.2;
  float pitchy = 3.2;
  double qeScint = 0.3;
  double qeCher = 0.2;
  double sigmaSPTR = 0.0;
  bool smearTimeStamp = false;
  int doiSlices = 10;
  double length = 15.0;
  double timeBinSize = 50.0; //ps
  std::string filePrefix = "";
  double pulseStart = 0.0;  // in ps
  double pulseEnd = 10000.0;   // in ps
  double spadRise = 200.0; // 200 ps
  double spadDecay = 10000.0; // 10 ns
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
      { "qeScint", required_argument, 0, 0 },
      { "sptr", required_argument, 0, 0 },
      { "doiSlices", required_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "qeCher", required_argument, 0, 0 },
      { "timeBinSize", required_argument, 0, 0 },
      { "spadPulse", no_argument, 0, 0 },
      { "spadRise", required_argument, 0, 0 },
      { "spadDecay", required_argument, 0, 0 },
      { "pulseEnd", required_argument, 0, 0 },
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
      numb_of_phot_for_time_average = atof((char *)optarg);
      // std::cout << "Time average on first " << numb_of_phot_for_time_average << " photons"<< std::endl;
    }
    else if (c == 0 && optionIndex == 1){
      // std::cout << "SiPM saturation will be taken into account " << std::endl;
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
      qeScint = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      sigmaSPTR = atof((char *)optarg);
      if(sigmaSPTR > 0)
      {
        smearTimeStamp = true;
      }
    }
    else if (c == 0 && optionIndex == 8){
      doiSlices = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      length = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      // std::cout << "Including Cherenkov photons " << std::endl;
      qeCher = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      timeBinSize = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      spadPulse = true;
    }
    else if (c == 0 && optionIndex == 13){
      spadRise = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      spadDecay = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      pulseEnd = atof((char *)optarg);
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

  // if the real spad pulses is not used, the timestamp will be calculated as average of first N photons,
  // but the the N has to be integer
  if(!spadPulse)
  {
    numb_of_phot_for_time_average = round(numb_of_phot_for_time_average);
  }

  int pulseBins = (int) round((pulseEnd-pulseStart)/(timeBinSize)); // /1000 because input timeBinSize is ps

  // feedback to user
  std::cout << "|----------------------------------------------|" << std::endl;
  std::cout << "|               USER PARAMETERS                |" << std::endl;
  std::cout << "|----------------------------------------------|" << std::endl;
  std::cout << std::endl;
  std::cout << "Input file prefix               = " << filePrefix << std::endl;
  std::cout << "Output file                     = " << outputFileName << std::endl;
  std::cout << "Time average on first           = " << numb_of_phot_for_time_average << " photons"<< std::endl;
  std::cout << "nmppcx                          = " << nmppcx << std::endl;
  std::cout << "nmppcy                          = " << nmppcy << std::endl;
  std::cout << "pitchx                          = " << pitchx << std::endl;
  std::cout << "pitchy                          = " << pitchy << std::endl;
  std::cout << "PDE scintillation ph            = " << qeScint<< std::endl;
  std::cout << "PDE Cherenkov ph                = " << qeCher << std::endl;
  std::cout << "SPTR [ps] sigma                 = " << sigmaSPTR << std::endl;
  std::cout << "Number of slices in z           = " << doiSlices << std::endl;
  std::cout << "Crystal length [mm]             = " << length << std::endl;
  std::cout << "Size of time bins [ps]          = " << timeBinSize << std::endl;
  if(saturation)
  {
    std::cout << "--> SiPM saturation will be taken into account " << std::endl;
  }
  if(spadPulse)
  {
    std::cout << "--> Pulse generated using spad pulses: " << std::endl;
    std::cout << "Spad pulse rise time [ps]       = " << spadRise << std::endl;
    std::cout << "Spad pulse decay time [ps]      = " << spadDecay << std::endl;
  }
  std::cout << std::endl;
  std::cout << "|----------------------------------------------|" << std::endl;
  std::cout << std::endl;

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
  // #pragma omp parallel for
  for(int i = 0 ; i < numOfCh ; i++)
  {
    //
    std::stringstream Htitle;
    Htitle << "Global CTR channel " << i;
    sipm[i].ctr = new TH1F(Htitle.str().c_str(),Htitle.str().c_str(),1000,pulseStart*1e-12,pulseEnd*1e-12); // in seconds
    sipm[i].ctr->GetXaxis()->SetTitle("Time [s]");
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

      //pulse
      std::stringstream Hname;
      Hname << "Ch "<< i << " Pulse z " << tempSlice.zMin << " to " << tempSlice.zMax ;


      tempSlice.pulse = new TH1F(Hname.str().c_str(),Hname.str().c_str(),pulseBins,pulseStart/1000.0,pulseEnd/1000.0); // in ns
      std::stringstream Htitle;
      Htitle << "Time [ns]";
      tempSlice.pulse->GetXaxis()->SetTitle(Htitle.str().c_str());
      Htitle.str("");
      Htitle << "Photons / " << timeBinSize << "ps";
      tempSlice.pulse->GetYaxis()->SetTitle(Htitle.str().c_str());

      //ctr
      Hname.str("");
      Hname << "Ch "<< i << " CTR z " << tempSlice.zMin << " to " << tempSlice.zMax ;
      tempSlice.ctr = new TH1F(Hname.str().c_str(),Hname.str().c_str(),1000,pulseStart*1e-12,pulseEnd*1e-12);
      Htitle.str("");
      Htitle << "Time [s]";
      tempSlice.ctr->GetXaxis()->SetTitle(Htitle.str().c_str());
      Htitle.str("");
      Htitle << "N";
      tempSlice.ctr->GetYaxis()->SetTitle(Htitle.str().c_str());

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

  //find max of single spad pulse
  std::vector<double> x,y;
  // double tEnd = 10.0;
  // double tBegin = 0.0;
  // int divisions = 10000;
  // double t0 = 1.0;
  // double tr = 0.2;
  // double td = 10;
  double maxSingleSpad = 0.0;
  for(int i = 0 ; i < pulseBins ; i++)
  {
    double valueX = i*((pulseEnd - pulseStart)/pulseBins)/1000.0;   //in ns
    double valueY = singleSpad(valueX,0.0,spadRise/1000.0,spadDecay/1000.0);
    if(valueY > maxSingleSpad)
    {
      maxSingleSpad = valueY;
    }
    // x.push_back(valueX);
    // y.push_back(valueY);
  }
  // std::cout << max << std::endl;


  //----------------------------------------//
  //             LOOP ON EVENTS             //
  //----------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  // #pragma omp parallel for
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
    // #pragma omp parallel for
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

    bool keepTheEvent = false;    // whether to keep or not the event, in terms of energyDeposition, single crystal hit and which crystal was hit
    int phSlice = 0;             // which slice of the chosen crystal was hit

    // count the "photoelectric" depositions in the slices, for the chosen crystal
    if(totalEnergyDeposited > 0.5 && CrystalsHit == 1 && crystals[0] == 36)
    {
      keepTheEvent = true;
      for(int iSipm = 0; iSipm < numOfCh ; iSipm++)
      {
        for(unsigned int iSlice = 0 ; iSlice < sipm[iSipm].slice.size(); iSlice++)
        {
          if(RealZ > sipm[iSipm].slice[iSlice].zMin && RealZ < sipm[iSipm].slice[iSlice].zMax)
          {
            phSlice = iSlice;
            sipm[iSipm].slice[iSlice].events++;
          }
        }
      }
    }
    //------------------------//
    // OPTICAL PHOTONS        //
    //------------------------//
    if(keepTheEvent)
    {
      // #pragma omp parallel for
      for(int iPhot = 0; iPhot < photons->size(); iPhot++) // run on all opticals
      {
        bool keepThePhoton = false;  // whether to keep the photon or not (for pde and saturation tests, so if the photons has triggers an avalanche)
        int phSipm = 0;              // which Sipm was hit by the optical photon
        // Float_t timeOfArrival = photons->at(iPhot).GlobalTime;
        // if(timeOfArrival > 2.0) break;

        // set PDE to the appropriate value, depending on photon type
        // if PDE = 0 for a given photon type, those photons will be discarded
        double qe;
        if(photons->at(iPhot).PhotonType == 1)// cherenkov photons
        {
          qe = qeCher;
        }
        else // all other photons are scintillation photons
        {
          qe = qeScint;
        }

        // find which sipm was hit and decide if keepThePhoton
        for(int iSipm = 0; iSipm < numOfCh ; iSipm++) // run on sipms
        {
          if((photons->at(iPhot).PositionX > (sipm[iSipm].detx - det_size_x/2.0 )) && (photons->at(iPhot).PositionX < (sipm[iSipm].detx + det_size_x/2.0 )) ) // check if x of photon is inside sensitive area of this sipm
          {
            if((photons->at(iPhot).PositionY > (sipm[iSipm].dety - det_size_y/2.0 )) && (photons->at(iPhot).PositionY < (sipm[iSipm].dety + det_size_y/2.0 )) )//check if y of photon is inside sensitive area of this sipm
            {
              phSipm = iSipm; // get Sipm index
              if(saturation)  // consider saturation if user wants so
              {
                // find which spad was hit
                // bring sipm hit to start in 0,0
                float hitx = photons->at(iPhot).PositionX - sipm[iSipm].detx + (det_size_x/2.0);
                float hity = photons->at(iPhot).PositionY - sipm[iSipm].dety + (det_size_y/2.0);
                // find spad I and J
                int hiti = (int) (hitx / spad_size_x);
                int hitj = (int) (hity / spad_size_y);
                // check if the photon has hit the dead area of the sipm (central part, for TSV)
                if( (hiti > ( n_spad_x/2 - n_dead_spad_x/2 - 1 ) ) &&
                (hiti < ( n_spad_x/2 + n_dead_spad_x/2 - 1 ) ) &&
                (hitj > ( n_spad_y/2 - n_dead_spad_y/2 - 1 ) ) &&
                (hitj < ( n_spad_y/2 + n_dead_spad_y/2 - 1 ) )   )
                {
                  // do nothing, this part of the sipm is not active
                }
                else //photon on active spad (but it could have been hit already)
                {
                  //HACK to avoid seg fault when the optical photon is exactly on the border (which makes the hiti of hitj being exatly 60 for example)
                  if(hiti == n_spad_x) hiti = hiti -1;
                  if(hitj == n_spad_y) hitj = hitj -1;

                  //quantum efficiency test
                  if(sipm[iSipm].spad[hiti][hitj] == 0) // if this spad was not hit yet
                  {
                    // qe test
                    double numb = rand->Uniform(1.0);
                    if(numb < qe) // this test involves already the photon nature (scintillation/cherenkov)
                    {
                      keepThePhoton = true; //keep the photon
                      sipm[iSipm].spad[hiti][hitj] = 1; // set the spad as dead
                    }
                  }
                  else // spad was already hit, so don't keepThePhoton
                  {
                    //ignore the hit, the spad has already fired and the optical photon is lost
                  }
                }
              }
              else // perform just qe test for this photon
              {
                double numb = rand->Uniform(1.0);
                if(numb < qe)  // accepted photon
                {
                  keepThePhoton = true; //keep the photon
                }
              }
            }
          }
        }

        // so if the photon has started an avalanche
        if(keepThePhoton)
        {
          sipm[phSipm].counts++;  // increment counts on the SiPM that was hit
          Float_t RealTimeStamp = photons->at(iPhot).GlobalTime;  // get time of arrival of the photon
          if(smearTimeStamp)  // smear time stamp if user says so, according to sigmaSPTR
          {
            RealTimeStamp = gRandom->Gaus(RealTimeStamp,sigmaSPTR/1000.0); // RealTimeStamp in ns, sigmaSPTR in ps, converted here
          }
          sipm[phSipm].listOfTimestamps.push_back(RealTimeStamp); // add to the list of timestamps for the sipm hit
          if(totalEnergyDeposited > 0.5 && CrystalsHit == 1 && crystals[0] == 36) // check if energy deposited is in photopeak and the crystal is the chosen one (36)
          {
            if(spadPulse) // sum to the SiPM pulse this spad pulse
            {
              for(int iBin = 0; iBin < sipm[phSipm].slice[phSlice].pulse->GetNbinsX(); iBin++) //run on all bins
              {
                double tBin = sipm[phSipm].slice[phSlice].pulse->GetBinCenter(iBin+1);
                if( tBin >= RealTimeStamp) //nothing happens before RealTimeStamp
                {
                  // sipm[phSipm].slice[phSlice].pulse->Fill(tBin,1.0*( exp(-(tBin - RealTimeStamp)/(spadDecay/1000.0)) * (1.0 - exp(-(tBin - RealTimeStamp)/(spadRise/1000.0)) ) ) );  //input of spadDecay and spadRise is in ps, all times are in ns so these two are converted here
                  sipm[phSipm].slice[phSlice].pulse->Fill(tBin,singleSpad(tBin,RealTimeStamp,spadRise/1000.0,spadDecay/1000.0)/maxSingleSpad );  //input of spadDecay and spadRise is in ps, all times are in ns so these two are converted here
                }
              }
            }
            else // simply build the pulse with the timestamps
            {
              sipm[phSipm].slice[phSlice].pulse->Fill(RealTimeStamp); // fill the pulse histo for the slice of this sipm
              // sipm[phSipm].slice[phSlice].listOfTimestamps.push_back(RealTimeStamp); // fill the listOfTimestamps for the slice of this sipm)
            }
          }
        }
      }
      // end of run on all photons of this gamma event

      // calculate the global (and per slice) sipm parameters
      for(int i = 0; i < numOfCh ; i++)
      {
        sipm[i].timestamp = 0.0;


        if(spadPulse)
        {
          double threshold = numb_of_phot_for_time_average;
          for(int iBin = 0; iBin < sipm[i].slice[phSlice].pulse->GetNbinsX()-1; iBin++) //run on all bins
          {
            double tBin = sipm[i].slice[phSlice].pulse->GetBinCenter(iBin+1);
            double vBin = sipm[i].slice[phSlice].pulse->GetBinContent(iBin+1);
            double tBinNext = sipm[i].slice[phSlice].pulse->GetBinCenter(iBin+2);
            double vBinNext = sipm[i].slice[phSlice].pulse->GetBinContent(iBin+2);
            if(vBin < threshold && vBinNext > threshold)
            {
              sipm[i].timestamp = 1e-9*(tBin + tBinNext)/2.0;
              sipm[i].ctr->Fill(sipm[i].timestamp);
              sipm[i].slice[phSlice].ctr->Fill(sipm[i].timestamp);
              break;
            }
          }

        }
        else
        {
          // calculate the sipm timestamp from average of first N timestamps
          // first, the timestamps need to be ordered again (because of gaussian smearing, if applied)
          std::sort(sipm[i].listOfTimestamps.begin(),sipm[i].listOfTimestamps.end());
          int effectiveN = numb_of_phot_for_time_average;
          if(numb_of_phot_for_time_average > sipm[i].listOfTimestamps.size())
            effectiveN = sipm[i].listOfTimestamps.size();
          for(int j = 0 ; j < effectiveN; j++)
          {
            sipm[i].timestamp +=  (Float_t) ( (sipm[i].listOfTimestamps[j] / effectiveN)*1e-9); // already smeared if user says so, and converted now to seconds
          }
          sipm[i].ctr->Fill(sipm[i].timestamp);
          sipm[i].slice[phSlice].ctr->Fill(sipm[i].timestamp);
        }



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

    sStack.str("");
    sStack << "Graph Pulses of Ch " << i;
    TMultiGraph *mg = new TMultiGraph();
    TCanvas *canvasMg = new TCanvas(sStack.str().c_str(),sStack.str().c_str(),1200,800);

    TLegend *legend = new TLegend(0.15,0.62,0.30,0.89,"");
    legend->SetFillStyle(0);

    TLegend *legendNorm = new TLegend(0.15,0.62,0.30,0.89,"");
    legendNorm->SetFillStyle(0);

    sipm[i].ctr->Write();

    // sname.str("");
    // sname << "No correction        = " << round(1e12*crystal[iCry].simpleCTR->GetRMS()*2.355*TMath::Sqrt(2.0)) << "ps";
    int color = 1; // avoid white


    for(unsigned int iSlice = 0 ; iSlice < sipm[i].slice.size(); iSlice++)
    {
      sipm[i].slice[iSlice].ctr->Write();
      std::stringstream sDoi;
      sDoi << "z = " <<(sipm[i].slice[iSlice].zMax + sipm[i].slice[iSlice].zMin) /2.0 << "mm";
      if(color == 5)
      {
        color++; // no yellow!!
      }


      sipm[i].slice[iSlice].pulse->SetLineColor(color);

      //normalize the histograms to the number of events in the doi slice - mandatory to average
      if(sipm[i].slice[iSlice].events != 0)
      {
        sipm[i].slice[iSlice].pulse->Scale(1.0/sipm[i].slice[iSlice].events);
      }

      //rescale pulse as multiple of single spad pulses
      //and get vectors
      std::vector<double> valueX;
      std::vector<double> valueY;
      for(int iBin = 0; iBin < sipm[i].slice[iSlice].pulse->GetNbinsX(); iBin++) //run on all bins
      {
        double Xvalue = sipm[i].slice[iSlice].pulse->GetBinCenter(iBin+1);
        double Yvalue = sipm[i].slice[iSlice].pulse->GetBinContent(iBin+1);
        valueX.push_back(Xvalue);
        valueY.push_back(Yvalue);
        // sipm[i].slice[iSlice].pulse->SetBinContent(iBin+1,Yvalue/maxSingleSpad);
      }
      TGraph *gr = new TGraph(sipm[i].slice[iSlice].pulse->GetNbinsX(),&valueX[0],&valueY[0]);
      gr->SetLineColor(color);
      gr->SetLineWidth(2);

      mg->Add(gr,"l");
      //clone the histo
      TH1F *clone = (TH1F*) sipm[i].slice[iSlice].pulse->Clone();
      hs->Add(clone);
      legend->AddEntry(clone,sDoi.str().c_str(),"l");

      //normalize histograms to max, by fitting to line
      // TF1 *line = new TF1("line","[0]",1.75,2.5); //FIXME hardcoded
      // int maxBin = sipm[i].slice[iSlice].pulse->GetMaximumBin();
      // double maxValue = sipm[i].slice[iSlice].pulse->GetBinContent(maxBin);
      // sipm[i].slice[iSlice].pulse->Fit(line,"RQN");
      // sipm[i].slice[iSlice].pulse->Scale(1.0/line->GetParameter(0));
      // hsNorm->Add(sipm[i].slice[iSlice].pulse);
      // legendNorm->AddEntry(sipm[i].slice[iSlice].pulse,sDoi.str().c_str(),"l");
      // sipm[i].slice[iSlice].pulse->Write();
      color++;
    }
    canvasMg->cd();
    mg->Draw("a");
    legend->Draw();
    canvasMg->Write();

    canvas->cd();
    hs->Draw("nostack");
    legend->Draw();
    canvas->Write();


    // canvasNorm->cd();
    // hsNorm->Draw("nostack");
    // legendNorm->Draw();
    // canvasNorm->Write();
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
// # endif
