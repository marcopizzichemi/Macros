// compile with
// g++ -o ../build/timeAnalysis timeAnalysis.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

// small program to extract timing calibration and data

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"
#include "TNamed.h"
#include "TPaveLabel.h"
#include "THStack.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort

#include <sys/types.h>
#include <dirent.h>

// typedef std::vector<std::string> stringvec;
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}


struct Crystal_t
{
  int number;
  int detectorChannel;
  int timingChannel;
  std::vector<int> relevantForW;
  std::vector<int> delayTimingChannels;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TGraph *wz;
  TH1F *simpleCTR;
  TH1F *centralCTR;
  TH1F *allCTR;
  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TTreeFormula *Formula;
  std::vector<double> z;
  TGraph* tw_correction;
  TGraph* rms_tw_correction;
  std::vector<TGraph*> delay;
  std::vector<TGraph*> rms_delay;
  const char* path;
  bool accepted;
  // TCanvas
};

struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};


/*** find width at half max ***/
float ComputeFWHM(TH1F* histo)
{
   float max = histo->GetMaximum();
   float halfMax = max / 2.0;
   int binMin = histo->FindFirstBinAbove(halfMax);
   int binMax = histo->FindLastBinAbove(halfMax);
   float down = histo->GetBinCenter(binMin);
   float up   = histo->GetBinCenter(binMax);
   float ret = up -down;
   return ret;
}

void extractCTR(TH1F* histo,double fitMin,double fitMax, int divs, double* res)
{
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  // std::cout << f1min << " " << f1max << std::endl;
  TF1* f1  = new TF1("f1","crystalball");
  f1->SetLineColor(kBlack);
  f1->SetParameters(histo->GetMaximum(),histo->GetMean(),histo->GetRMS(),1,3);
  histo->Fit(f1,"Q","",fitMin,fitMax);
  double min,max,min10,max10;
  // int divs = 3000;
  double step = (f1max-f1min)/divs;
  // is [0] the max of the function???
  for(int i = 0 ; i < divs ; i++)
  {
    if( (f1->Eval(f1min + i*step) < f1->GetParameter(0)/2.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetParameter(0)/2.0) )
    {
      min = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > f1->GetParameter(0)/2.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetParameter(0)/2.0) )
    {
      max = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) < f1->GetParameter(0)/10.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetParameter(0)/10.0) )
    {
      min10 = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > f1->GetParameter(0)/10.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetParameter(0)/10.0) )
    {
      max10 = f1min + (i+0.5)*step;
    }
  }
  res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(70e-12,2));
  res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((70e-12/2.355)*4.29,2));
}

/*** find effective sigma ***/
void FindSmallestInterval(float* ret, TH1F* histo, const float&
fraction, const bool& verbosity)
{
   float integralMax = fraction * histo->Integral();

   int N = histo -> GetNbinsX();
   std::vector<float> binCenters(N);
   std::vector<float> binContents(N);
   std::vector<float> binIntegrals(N);
   for(int bin1 = 0; bin1 < N; ++bin1)
   {
     binCenters[bin1] = histo->GetBinCenter(bin1+1);
     binContents[bin1] = histo->GetBinContent(bin1+1);

     for(int bin2 = 0; bin2 <= bin1; ++bin2)
       binIntegrals[bin1] += binContents[bin2];
   }

   float min = 0.;
   float max = 0.;
   float delta = 999999.;
   for(int bin1 = 0; bin1 < N; ++bin1)
   {
     for(int bin2 = bin1+1; bin2 < N; ++bin2)
     {
       if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;

       float tmpMin = histo -> GetBinCenter(bin1);
       float tmpMax = histo -> GetBinCenter(bin2);

       if( (tmpMax-tmpMin) < delta )
       {
         delta = (tmpMax - tmpMin);
         min = tmpMin;
         max = tmpMax;
       }

       break;
     }
   }

   TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
   for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
   {
     if( smallHisto->GetBinCenter(bin) < min )
       smallHisto -> SetBinContent(bin,0);

     if( smallHisto->GetBinCenter(bin) > max )
       smallHisto -> SetBinContent(bin,0);
   }
   smallHisto -> SetFillColor(kYellow);

   float mean = smallHisto -> GetMean();
   float meanErr = smallHisto -> GetMeanError();

   ret[0] = mean;
   ret[1] = meanErr;
   ret[2] = min;
   ret[3] = max;
}


bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
{
  return a.number < b.number;
}

void usage()
{
  std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
            << "\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
            // << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
            << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
            << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
            << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
            << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 70e-12"  << std::endl
            << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
            << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            << "\t\t" << "--smooth <value>                                 - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
            << "\t\t" << "--fitMin <value>                                 - min of crystalball fit in sec - default = -9.4e-9"  << std::endl
            << "\t\t" << "--fitMax <value>                                 - max of crystalball fit in sec - default = -8.1e-9"  << std::endl
            << "\t\t" << "--divs <value>                                   - n of divisions when looking for FWHM - default = 10000"  << std::endl
            << "\t\t" << "--bins <value>                                   - n of bins in summary CTR histograms - deafult 40"  << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  std::string inputFileName = "";
  std::string outputFileName = "";
  std::string calibrationFileName = "";
  // std::string coincidenceCalibrationFileName = "";

  bool simulation = false;
  Float_t length = 15.0; //mm
  Float_t doiFraction = 0.5;
  Float_t tagFwhm = 70.0e-12; //s
  Float_t rmsLow = 1.75;
  Float_t rmsHigh = 1.75;
  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s
  Float_t fitMin = -9.4e-9; // s
  Float_t fitMax = -8.1e-9; //s
  int divs       = 10000;
  int histoBins = 500;
  int smooth = 0; //
  int bins = 40;
  double minCTR = 100;
  double maxCTR = 500;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "simulation", no_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "doiFraction", required_argument, 0, 0 },
      // { "coincidence", required_argument, 0, 0 },
      { "tagFwhm", required_argument, 0, 0 },
      { "rmsLow", required_argument, 0, 0 },
      { "rmsHigh", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "smooth", required_argument, 0, 0 },
      { "fitMin", required_argument, 0, 0 },
      { "fitMax", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "bins", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:c:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 'c'){
      calibrationFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      calibrationFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      std::cout << "Dataset from simulation " << std::endl;
      simulation = true;
    }
    else if (c == 0 && optionIndex == 4){
      length = atof((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 5){
      doiFraction = atof((char *)optarg);;
    }
    // else if (c == 0 && optionIndex == 6){
    //   coincidenceCalibrationFileName = (char *)optarg;
    // }
    else if (c == 0 && optionIndex == 6){
      tagFwhm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      rmsLow = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      rmsHigh = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      smooth = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      fitMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      fitMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      divs = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 16){
      bins = atoi((char *)optarg);
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;

  // read file in dir
  std::vector<std::string> v;
  read_directory(".", v);
  // std::copy(v.begin(), v.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;

  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(v[i]);
    }
  }



  // TFile *treeFile = new TFile(inputFileName.c_str());
  // TTree* tree = (TTree*) treeFile->Get("adc");

  //----------------------------------------------------------//
  //  Get TChain of input files                               //
  //----------------------------------------------------------//
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }

  // Add several TTreeFormula to the list;
  TList formulas;

  //variables for the input TChain
  // ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  // ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous event
  // Int_t        *ChainAdcChannel;
  // Short_t      *ChainDesktopAdcChannel;                              // input TChain data for desktop digitizers - data is int_16
  // UShort_t     *ChainVMEadcChannel;                                  // input TChain data for VME digitizers - data is uint_16
  // Float_t      *ChainTimeStamp;
  // Float_t      *TDCBinning;
  // // Short_t      *ChainPetirocChannel;                                 //FIXME temporary data type of petiroc charge input - ask
  // Float_t       RealX;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       RealY;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       RealZ;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       simTaggingCharge;
  // Float_t       simTaggingTime;
  // Short_t       CrystalsHit;                                         // "real" number of crystals hit in the event (from simulation data)
  // Short_t       NumbOfInteractions;                                  // "real" number of interaction (energy depositions) in the event (from simulation data)
  //
  // //branches for the input TChain
  // TBranch      *bChainExtendedTimeTag;                               // branches for above data
  // TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  // TBranch     **bChainAdcChannel;                                    // branches for above data
  // TBranch     **bChainTimeStamp;
  // TBranch      *bRealX;                                              // branches for above data
  // TBranch      *bRealY;                                              // branches for above data
  // TBranch      *bRealZ;                                              // branches for above data
  // TBranch      *bsimTaggingCharge;                                              // branches for above data
  // TBranch      *bsimTaggingTime;                                              // branches for above data
  // TBranch      *bCrystalsHit;                                        // branches for above data
  // TBranch      *bNumbOfInteractions;                                 // branches for above data
  // TBranch      *bTotalCryEnergy;                                     //
  //
  // if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  // {
  //   for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
  //   {
  //     std::cout << "Adding file " << argv[i] << std::endl;
  //     tree->Add(argv[i]);
  //   }
  // }
  // else // the config file was indeed the default one
  // {
  //   for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
  //   {
  //     std::cout << "Adding file " << argv[i] << std::endl;
  //     tree->Add(argv[i]);
  //   }
  // }

  std::vector<int> detector_channels;

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
  std::string ch_prefix("ch");
  std::string t_prefix("t");

  // std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
    // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
    // numOfCry++;
  }
  //the string "cry" appears 4 times per crystal..
  // numOfCry = numOfCry / 4;
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  // std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;


  // first, create the adc channels variables and branches
  // ChainAdcChannel        = new Int_t [numOfCh];
  // ChainDesktopAdcChannel = new Short_t [numOfCh]; // input from ADC desktop
  // ChainVMEadcChannel     = new UShort_t [numOfCh]; // input from VME
  // ChainTimeStamp         = new Float_t[numOfCh];
  // // TDCBinning             = new Float_t[numOfCh];
  // // DigitizerChannelOn     = new bool[adcChannels];
  // bChainAdcChannel       = new TBranch* [numOfCh];
  // bChainTimeStamp        = new TBranch* [numOfCh];
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  UShort_t      *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;

  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfCh];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfCh];

  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  // if(usingRealSimData)
  // {
  //   tree->SetBranchAddress("RealX", &RealX, &bRealX);
  //   tree->SetBranchAddress("RealY", &RealY, &bRealY);
  //   tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
  //   // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
  //   // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
  //   tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
  //   tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
  //   // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
  // }
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    // stype << "ch" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    // stype.str("");

    sname << "t" << detector_channels[i];
    // stype << "t" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
    // stype.str("");
  }



  // find the crystals with complete calibration data
  // this means looking into TWO files:
  // 1. the standard calibration file, but produced with active timing part
  //    this is produced NOT in coincidence and will give all but two plot
  //    for the crystals calibration. the only two missing are the only two that HAVE
  //    to be measured in coincidence, i.e. the tw_correction plots (and it's RMS) and the tag photopeak cut
  // 2. a configuration run performed in coincidence, from which the tw_correction are taken
  //
  // afterwards the crystals will be accepted if they are present in both calibration files

  std::vector<Crystal_t> crystal;
  std::vector<detector_t> detectorSaturation;

  // STANDARD CALIBRATION FILE
  TFile *calibrationFile = new TFile(calibrationFileName.c_str());
  calibrationFile->cd("Module 0.0");
  TList *listModule = gDirectory->GetListOfKeys();
  int nKeysMod = listModule->GetEntries();
  std::vector<std::string> keysModName;
  // fill a vector with the leaves names
  std::string mppc_prefix("MPPC");
  for(int i = 0 ; i < nKeysMod ; i++){
    keysModName.push_back(listModule->At(i)->GetName());
  }

  TCut *taggingPhotopeakCut;
  int taggingCrystalTimingChannel;
  std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  std::string taggingCrystalTimingChannel_prefix("taggingCrystalTimingChannel");
  // std::string det_prefix("channels");
  // std::string saturation_prefix("saturation");

  std::vector<int> *pChannels;
  std::vector<float> *pSaturation;
  std::vector<float> *pPedestal;
  gDirectory->GetObject("channels",pChannels);
  gDirectory->GetObject("saturation",pSaturation);
  gDirectory->GetObject("pedestal",pPedestal);

  std::vector<int> DetChannels = pChannels[0];
  std::vector<float> saturation = pSaturation[0];
  std::vector<float> pedestal = pPedestal[0];

  for(unsigned int iSat = 0; iSat < DetChannels.size(); iSat++)
  {
    detector_t tempDetector;
    tempDetector.digitizerChannel = DetChannels[iSat];
    tempDetector.saturation = saturation[iSat];
    tempDetector.pedestal = pedestal[iSat];
    detectorSaturation.push_back(tempDetector);
  }

  std::vector<std::string> MPPCfolders;
  for(unsigned int i = 0 ; i < keysModName.size() ; i++)
  {
    if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
    {
      MPPCfolders.push_back(keysModName[i]);
    }
    if(!keysModName[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
    {
//       std::cout << keysCryName[i] << std::endl;
      taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
//       if(cut)
//         temp_crystal.CrystalCut = cut;
    }
    if(!keysModName[i].compare(0,taggingCrystalTimingChannel_prefix.size(),taggingCrystalTimingChannel_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      taggingCrystalTimingChannel = atoi(snameCh.str().c_str());
    }

    // if(!keysModName[i].compare(0,det_prefix.size(),det_prefix)) // find tcut
    // {
    //   std::stringstream snameCh;
    //   std::vector<int> *v;
    //   gDirectory->GetObject("channelsNumRelevantForW",v);
    //   temp_crystal.relevantForW = v[0];
    // }
    //
    // if(!keysModName[i].compare(0,saturation_prefix.size(),saturation_prefix)) // find tcut
    // {
    //   std::stringstream snameCh;
    //   std::vector<int> *v;
    //   gDirectory->GetObject("channelsNumRelevantForW",v);
    //   temp_crystal.relevantForW = v[0];
    // }


  }

  std::stringstream sformulaname;
  sformulaname << "FormulaTag";
  TCut taggingPhotopeakCutName;
  taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // std::cout << "FormulaTag ------------- \n" << taggingPhotopeakCutName << std::endl;
  TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);
  formulas.Add(FormulaTag);



  for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
  {
    // std::cout << MPPCfolders[iMppc] << std::endl;
    gDirectory->cd(MPPCfolders[iMppc].c_str());
    TList *listMppc = gDirectory->GetListOfKeys();
    int nKeysMppc = listMppc->GetEntries();
    std::vector<std::string> keysMppcName;
    // fill a vector with the leaves names
    std::string crystal_prefix("Crystal");
    // std::string det_prefix("digitizerChannel");
    // std::string saturation_prefix("saturation");
    for(int i = 0 ; i < nKeysMppc ; i++){
      keysMppcName.push_back(listMppc->At(i)->GetName());
    }

    std::vector<std::string> CrystalFolders;

    for(unsigned int i = 0 ; i < keysMppcName.size() ; i++)
    {
      if (!keysMppcName[i].compare(0, crystal_prefix.size(), crystal_prefix))
      {
        CrystalFolders.push_back(keysMppcName[i]);
      }
      // if (!keysMppcName[i].compare(0, det_prefix.size(), det_prefix))
      // {
      //   std::stringstream snameCh;
      //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
      //   tempDetector.digitizerChannel = atoi(snameCh.str().c_str());
      // }
      // if (!keysMppcName[i].compare(0, saturation_prefix.size(), saturation_prefix))
      // {
      //   std::stringstream snameCh;
      //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
      //   tempDetector.saturation = atof(snameCh.str().c_str());
      // }
    }
    // detectorSaturation.push_back(tempDetector);

    for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
    {
      //  std::cout << CrystalFolders[iCry] << std::endl;
       gDirectory->cd(CrystalFolders[iCry].c_str());

       Crystal_t temp_crystal;
       temp_crystal.number = -1;
       temp_crystal.CrystalCut = NULL;
       temp_crystal.CrystalCutWithoutCutG = NULL;
       temp_crystal.PhotopeakEnergyCut = NULL;
       temp_crystal.calibrationGraph = NULL;
       temp_crystal.simpleCTR = NULL;
       temp_crystal.centralCTR = NULL;
       temp_crystal.allCTR = NULL;
       temp_crystal.simpleCTR_norm = NULL;
       temp_crystal.centralCTR_norm = NULL;
       temp_crystal.allCTR_norm = NULL;
       temp_crystal.wz = NULL;
       temp_crystal.accepted = true;
       temp_crystal.tw_correction = NULL;


       //get crystal number
       temp_crystal.number = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());

       TList *listCry = gDirectory->GetListOfKeys();
       int nKeysCry = listCry->GetEntries();
       std::vector<std::string> keysCryName;
       if(nKeysCry) //if directory not empty
       {

         for(int i = 0 ; i < nKeysCry ; i++){
           keysCryName.push_back(listCry->At(i)->GetName());
         }
         std::string CalibName;
         std::string CutName;
         std::vector<std::string> cutgNames;
         std::string cutG_prefix("cutg");
         std::string calibration_prefix("Calibration");
         std::string crystalCut_prefix("CrystalCut");
         std::string crystalCutWithoutCutG_prefix("CrystalCutWithoutCutG");
         std::string photopeakEnergyCut_prefix("PhotopeakEnergyCut");
         std::string channel_prefix("digitizerChannel");
         std::string w_channels_prefix("channelsNumRelevantForW");
         std::string timing_channel_prefix("timingChannel");
         std::string wz_prefix("w(z)");
        //  std::string correction_prefix("Correction Graph");
        //  std::string correction_rms_prefix("RMS Correction Graphs");
         for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
         {
           if(!keysCryName[i].compare(0,calibration_prefix.size(),calibration_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.calibrationGraph = calibGraph;
           }

           if(!keysCryName[i].compare(0,wz_prefix.size(),wz_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.wz = calibGraph;
           }

           if(!keysCryName[i].compare(0,crystalCut_prefix.size(),crystalCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCut = cut;
           }
           if(!keysCryName[i].compare(0,crystalCutWithoutCutG_prefix.size(),crystalCutWithoutCutG_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCutWithoutCutG = cut;
           }
           if(!keysCryName[i].compare(0,photopeakEnergyCut_prefix.size(),photopeakEnergyCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.PhotopeakEnergyCut = cut;
           }

           if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCutG* cutg = (TCutG*) gDirectory->Get(keysCryName[i].c_str());

             temp_crystal.cutg.push_back(cutg);
           }

           if(!keysCryName[i].compare(0,channel_prefix.size(),channel_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.detectorChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }




           if(!keysCryName[i].compare(0,timing_channel_prefix.size(),timing_channel_prefix)) // find timing channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.timingChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,w_channels_prefix.size(),w_channels_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             std::vector<int> *v;
             gDirectory->GetObject("channelsNumRelevantForW",v);
             // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.relevantForW = v[0];
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }
         }

         bool dirExists;
         std::stringstream sname;
         dirExists = gDirectory->cd("TimeCorrection");
         if(dirExists)
         {
           sname.str("");
           sname << "Central correction - Crystal " << temp_crystal.number;
           temp_crystal.centralCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           sname.str("");
           sname << "Full correction - Crystal " << temp_crystal.number;
           temp_crystal.allCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           sname.str("");

           temp_crystal.path = gDirectory->GetPath();
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
             std::string graph_delay_prefix = "Graph Delay ch_";
             std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
               {
                 TGraph *calibGraph = NULL;
                 calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                   temp_crystal.tw_correction = calibGraph;
               }

               if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
               {
                TGraph *calibGraph = NULL;
                calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
                if(calibGraph)
                   temp_crystal.rms_tw_correction = calibGraph;
               }

               if(!keysTcorrName[i].compare(0,graph_delay_prefix.size(),graph_delay_prefix))
               {
                 TGraph *calibGraph = NULL;
                 calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                   temp_crystal.delay.push_back(calibGraph);
               }

               if(!keysTcorrName[i].compare(0,rms_graph_delay_prefix.size(),rms_graph_delay_prefix))
               {
                 TGraph *calibGraph = NULL;
                 calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                   temp_crystal.rms_delay.push_back(calibGraph);
               }

               if(!keysTcorrName[i].compare(0,delay_timing_ch_prefix.size(),delay_timing_ch_prefix))
               {
                 //  std::cout << keysCryName[i] << std::endl;
                 std::stringstream snameCh;
                 std::vector<int> *v;
                 gDirectory->GetObject("delayTimingChannels",v);
                 // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                 //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                //  istringstream()
                 temp_crystal.delayTimingChannels = v[0];
                 //  std::cout <<temp_crystal.detectorChannel << std::endl;
                //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                          //  << temp_crystal.detectorChannel << std::endl;
               }
             }
           }
           gDirectory->cd("..");
         }

         sname << "No correction - Crystal " << temp_crystal.number;
         temp_crystal.simpleCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);

         TCut globalCut; // the cut for the formula
         globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }

         sname << "Formula" << temp_crystal.number;
         TTreeFormula* Formula = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         formulas.Add(Formula);
         temp_crystal.Formula = Formula;
         sname.str("");

         if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut && (temp_crystal.cutg.size() == 2))
         {

           crystal.push_back(temp_crystal);
         }

       }
       gDirectory->cd("..");
    }
    calibrationFile->cd("Module 0.0");
  }

  // //DEBUG
  // for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  // {
  //   std::cout << detectorSaturation[iSat].digitizerChannel << " " << detectorSaturation[iSat].saturation << " " << detectorSaturation[iSat].pedestal << std::endl;
  // }

  // COINCIDENCE CALIBRATION FILE
  //
  // TFile *coincidenceCalibrationFile = new TFile(coincidenceCalibrationFileName.c_str());
  // coincidenceCalibrationFile->cd("Module 0.0");
  // TList *listModuleCoinc = gDirectory->GetListOfKeys();
  // int nKeysModCoinc = listModuleCoinc->GetEntries();
  // std::vector<std::string> keysModNameCoinc;
  // // fill a vector with the leaves names
  // // std::string mppc_prefixCoinc("MPPC");
  // for(int i = 0 ; i < nKeysModCoinc ; i++){
  //   keysModNameCoinc.push_back(listModuleCoinc->At(i)->GetName());
  // }

//   TCut *taggingPhotopeakCut;
//   std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  // std::vector<std::string> MPPCfoldersCoinc;
//   for(unsigned int i = 0 ; i < keysModNameCoinc.size() ; i++)
//   {
    // if (!keysModNameCoinc[i].compare(0, mppc_prefixCoinc.size(), mppc_prefixCoinc))
    // {
    //   MPPCfoldersCoinc.push_back(keysModNameCoinc[i]);
    // }
//     if(!keysModNameCoinc[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
//     {
     //  std::cout << keysCryName[i] << std::endl;
//       taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
      // if(cut)
      //   temp_crystal.CrystalCut = cut;
//     }
//   }

  // std::stringstream sformulaname;
  // sformulaname << "FormulaTag";
  // TCut taggingPhotopeakCutName;
  // taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);

  //look only into crystals found in the other calibration file
//   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
//   {
//     bool dirFound = gDirectory->cd(crystal[iCry].path);
//     if(dirFound)
//     {
//       crystal[iCry].accepted = true;
//       TList *listTcorr = gDirectory->GetListOfKeys();
//       int nKeysTcorr = listTcorr->GetEntries();
//       std::vector<std::string> keysTcorrName;
//       if(nKeysTcorr) //if directory not empty
//       {
//         for(int i = 0 ; i < nKeysTcorr ; i++){
//           keysTcorrName.push_back(listTcorr->At(i)->GetName());
//         }
//
//         std::string deltaWGraph_prefix = "DeltaW Graph";
//         std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
//         // std::string graph_delay_prefix = "Graph Delay ch_";
//         // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
//
//         for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
//         {
//           if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
//           {
//             TGraph *calibGraph = NULL;
//             calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//             if(calibGraph)
//               crystal[iCry].tw_correction = calibGraph;
//           }
//
//           if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
//           {
//            TGraph *calibGraph = NULL;
//            calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//            if(calibGraph)
//               crystal[iCry].rms_tw_correction = calibGraph;
//           }
//         }
//       }
//     }
//   }



  // // //BEGIN of DEBUG
  // for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  // {
  //   std::cout << crystal[i].number << "\t"
  //             // << crystal[i].cut->GetTitle() << "\t"
  //             << crystal[i].calibrationGraph->GetName() << "\t";
  //   for(unsigned int j = 0 ; j < crystal[i].cutg.size(); j++)
  //   {
  //     std::cout << crystal[i].cutg[j]->GetName() << "\t";
  //   }
  //   std::cout << std::endl;
  // }
  // // //END of DEBUG


  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
    }
  }

  //notify TTreeFormula(s) to TChain
  tree->SetNotify(&formulas);

  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;
  long int goodEvents = 0;
  long int counter = 0;
  // for (long long int i=0;i<1000000;i++)
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
        {
          if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
          {

            goodEvents++;

            //temp commented
            Float_t centralcorrection = 0.0;
            Float_t zeroCorrection    = 0.0;
            //no corr
            crystal[iCry].simpleCTR->Fill(timeStamp[crystal[iCry].timingChannel] -
                                          timeStamp[taggingCrystalTimingChannel]);

            if(crystal[iCry].tw_correction)
            {

              //calculate FloodZ...
              Float_t FloodZ;
              float centralChargeOriginal;
              float centralSaturation;
              float centralPedestal;
              Float_t division = 0.0;

              centralChargeOriginal = charge[crystal[iCry].detectorChannel];
              for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
              {
                if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                {
                  centralSaturation = detectorSaturation[iSat].saturation;
                  centralPedestal = detectorSaturation[iSat].pedestal;
                }
              }
              float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

              for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
              {
                float originalCh = charge[crystal[iCry].relevantForW[iW]];
                float saturationCh;
                float pedestalCorr;
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                  {
                    saturationCh = detectorSaturation[iSat].saturation;
                    pedestalCorr = detectorSaturation[iSat].pedestal;
                  }
                }
                division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
              }

              FloodZ = centralChargeCorr / division;

              // central corr
              // std::string deltaWGraph_prefix = "DeltaW Graph";
              // std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";


              Float_t averageTimeStamp = 0.0;
              Float_t totalWeight = 0.0;
              // averageTimeStamp += timeStamp[crystal[iCry].detectorChannel];
              centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
              crystal[iCry].centralCTR->Fill((timeStamp[crystal[iCry].timingChannel] + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel]);


              if(crystal[iCry].delay.size())
              {
                // full corr v0
                // first evalutate all ts as if they where read by the central crystal, i.e. correct all the neighboring channels
                // using the delta
                //
                Float_t weight = 0.0;
                weight = pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),-2);
                averageTimeStamp += weight * timeStamp[crystal[iCry].timingChannel];
                totalWeight += weight;
                //
                //full corr
                // zeroCorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
                //
                // Float_t weight = 0.0;
                // weight = pow(sqrt(pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),2)+pow(crystal[iCry].rms_tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)),2)),-2);
                //
                // averageTimeStamp += weight*(timeStamp[crystal[iCry].timingChannel]- zeroCorrection);
                // totalWeight += weight;
                // // std::cout << i << " " << iCry << " " <<  crystal[iCry].number << "\n";
                for(unsigned int iGraph = 0; iGraph < crystal[iCry].delay.size();iGraph++)
                {

                  std::string graphName = crystal[iCry].delay[iGraph]->GetName();
                  std::string rmsName = crystal[iCry].rms_delay[iGraph]->GetName();

                  std::size_t foundGraph = graphName.find_last_of("_");
                  std::size_t foundRms   = rmsName.find_last_of("_");

                  std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
                  std::string tChannelStringFromRms   = rmsName.substr(foundRms  +1);

                  // std::stringstream sname;
                  // sname << "Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
                  // std::string graph_delay_prefix = sname.str();
                  // std::cout << tChannelStringFromGraph << std::endl;
                  // std::cout << tChannelStringFromRms << std::endl;

                  // sname.str("");
                  // sname << "RMS Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
                  // std::string rms_graph_delay_prefix = sname.str();
                  // std::cout << tChannelStringFromRms << std::endl;
                  // sname.str("");






                  int graphCh = atoi(tChannelStringFromGraph.c_str() );
                  // std::cout << graphCh  << "\n";


                  int rmsCh   = atoi( tChannelStringFromRms.c_str() );
                  if(graphCh != rmsCh)
                  {
                    std::cout << "ERROR! TGraphs of delay and rms are from different timing channels!!!!" << std::endl;
                    break;
                  }
                  // std::cout << rmsCh  << "\n";
                  weight = pow(crystal[iCry].rms_delay[iGraph]->Eval(FloodZ),-2);
                  totalWeight += weight;
                  averageTimeStamp += weight*(timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ));
                  // std::cout << timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ) << "\t";
                }
                averageTimeStamp = averageTimeStamp/totalWeight;
                // then correct the average of ts for central correction
                crystal[iCry].allCTR->Fill((averageTimeStamp + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel]);
                // crystal[iCry].allCTR->Fill(averageTimeStamp + zeroCorrection  - timeStamp[taggingCrystalTimingChannel]);
              }
            }
            // end of temp commented
          }
        }
      }
      // delete Formula;
    }
    counter++;

    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }
  std::cout << "Good events = " << goodEvents << std::endl;
  std::sort(crystal.begin(), crystal.end(), compareByNumber);


  TH1F* noCorr = new TH1F("No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* centralCorr = new TH1F("Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* fullCorr = new TH1F("Full Correction","Full Correction",bins,minCTR,maxCTR);
  // std::vector<TH1F*> histograms;

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {

    std::stringstream sname;



    if(smooth)
    {
      if(crystal[iCry].simpleCTR)
      {
        crystal[iCry].simpleCTR ->Smooth(smooth);
      }
      if(crystal[iCry].centralCTR)
      {
        crystal[iCry].centralCTR ->Smooth(smooth);
      }
      if(crystal[iCry].allCTR)
      {
        crystal[iCry].allCTR ->Smooth(smooth);
      }
    }


    Float_t realBasicCTRfwhm,realBasicCTRfwtm ;
    Float_t realCentralCTRfwhm,realCentralCTRfwtm;
    Float_t realAllCTRfwhm,realAllCTRfwtm;
    double ret[2];

    // std::cout << "BASIC CTRs --------------------" << std::endl;

    if(crystal[iCry].simpleCTR)
    {
      crystal[iCry].simpleCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].simpleCTR->SetFillStyle(3001);
      crystal[iCry].simpleCTR->SetFillColor(kGreen);
      crystal[iCry].simpleCTR->SetLineColor(kGreen);
      crystal[iCry].simpleCTR->SetStats(0);
      crystal[iCry].simpleCTR_norm = (TH1F*) crystal[iCry].simpleCTR->Clone();
      extractCTR(crystal[iCry].simpleCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].simpleCTR->GetName() << "\t\t\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realBasicCTRfwhm = ret[0]*1e12;
      realBasicCTRfwtm = ret[1]*1e12;
      noCorr->Fill(ret[0]*1e12);

      crystal[iCry].simpleCTR->Write();
      crystal[iCry].simpleCTR_norm->Scale(1.0/crystal[iCry].simpleCTR_norm->GetMaximum());
    }

    if(crystal[iCry].centralCTR)
    {
      crystal[iCry].centralCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].centralCTR->SetFillStyle(3001);
      crystal[iCry].centralCTR->SetFillColor(kBlue);
      crystal[iCry].centralCTR->SetLineColor(kBlue);
      crystal[iCry].centralCTR->SetStats(0);
      crystal[iCry].centralCTR_norm = (TH1F*) crystal[iCry].centralCTR->Clone();
      extractCTR(crystal[iCry].centralCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].centralCTR->GetName() << "\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realCentralCTRfwhm = ret[0]*1e12;
      realCentralCTRfwtm = ret[1]*1e12;
      centralCorr->Fill(ret[0]*1e12);

      crystal[iCry].centralCTR->Write();
      // crystal[iCry].centralCTR->Scale(1.0/crystal[iCry].centralCTR->GetMaximum());
      crystal[iCry].centralCTR_norm->Scale(1.0/crystal[iCry].centralCTR_norm->GetMaximum());
    }

    if(crystal[iCry].allCTR)
    {
      crystal[iCry].allCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].allCTR->SetFillStyle(3001);
      crystal[iCry].allCTR->SetFillColor(kRed);
      crystal[iCry].allCTR->SetLineColor(kRed);
      crystal[iCry].allCTR->SetStats(0);
      crystal[iCry].allCTR->Write();
      crystal[iCry].allCTR_norm = (TH1F*) crystal[iCry].allCTR->Clone();
      extractCTR(crystal[iCry].allCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].allCTR->GetName() << "\t\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realAllCTRfwhm = ret[0]*1e12;
      realAllCTRfwtm = ret[1]*1e12;
      fullCorr->Fill(ret[0]*1e12);

      crystal[iCry].allCTR_norm->Scale(1.0/crystal[iCry].allCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());
    }

    sname.str("");

    sname << "Summary - Crystal " << crystal[iCry].number;
    TCanvas* c_summary = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    c_summary->cd();
    THStack *hs = new THStack("hs","");
    hs->Add(crystal[iCry].simpleCTR_norm);
    hs->Add(crystal[iCry].centralCTR_norm);
    hs->Add(crystal[iCry].allCTR_norm);


    // std::cout << "Crystal " << crystal[iCry].number << std::endl;
    // std::cout << "CTR FWHM [ps] " << std::endl;
    hs->Draw("hist nostack");
    sname.str("");
    sname << "CTR - Crystal " << crystal[iCry].number << " - width in FWHM";
    hs->SetTitle(sname.str().c_str());
    hs->GetXaxis()->SetTitle("Time [s]");
    hs->GetXaxis()->SetTitleOffset(1);
    hs->GetXaxis()->SetTitleSize(0.045);
    hs->GetXaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetLabelSize(0.045);
    TLegend *legend = new TLegend(0.54,0.62,0.89,0.89,"");
    legend->SetFillStyle(0);
    if(crystal[iCry].simpleCTR)
    {
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].simpleCTR,sname.str().c_str(),"f");
      // std::cout << "No correction       = "<< realBasicCTRfwhm   << " ps" << std::endl;
    }
    if(crystal[iCry].centralCTR)
    {
      sname.str("");
      sname << "Central correction = " << realCentralCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].centralCTR,sname.str().c_str(),"f");
      // std::cout << "Central correction  = "<< realCentralCTRfwhm << " ps" << std::endl;
    }
    if(crystal[iCry].allCTR)
    {
      sname.str("");
      sname << "Full correction       = " << realAllCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].allCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }

    sname.str("");
    legend->Draw();
    gStyle->SetOptTitle(0);
    TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","brndc");
    title->Draw();
    // std::cout << std::endl;

    c_summary->Write();


    TH1F* cloneBasic;
    TH1F* cloneCentral;
    TH1F* cloneAll;
    THStack *cloneHs = (THStack*) hs->Clone();
    TLegend *legend1 = new TLegend(0.15,0.69,0.49,0.89,"");
    legend1->SetFillStyle(0);
    sname.str("");
    sname << "Multi - Crystal " << crystal[iCry].number;
    TCanvas* c_multi = new TCanvas(sname.str().c_str(),sname.str().c_str(),1800,1400);
    c_multi->Divide(2,2);

    if(crystal[iCry].simpleCTR_norm)
    {
      cloneBasic   = (TH1F*) crystal[iCry].simpleCTR->Clone();
      c_multi->cd(1);
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend1->AddEntry(cloneBasic,sname.str().c_str(),"f");
      cloneBasic->Draw();
      legend1->Draw();
    }
    if(crystal[iCry].centralCTR_norm)
    {
      cloneCentral = (TH1F*) crystal[iCry].centralCTR->Clone();
      c_multi->cd(2);
      TLegend *legend2 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend2->SetFillStyle(0);
      sname.str("");
      sname << "Central correction   = " << realCentralCTRfwhm << "ps";
      legend2->AddEntry(cloneCentral,sname.str().c_str(),"f");
      cloneCentral->Draw();
      legend2->Draw();
    }
    if(crystal[iCry].allCTR_norm)
    {
      cloneAll     = (TH1F*) crystal[iCry].allCTR->Clone();
      c_multi->cd(3);
      TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend3->SetFillStyle(0);
      sname.str("");
      sname << "Full correction      = " << realAllCTRfwhm << "ps";
      legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
      cloneAll->Draw();
      legend3->Draw();
    }







    c_multi->cd(4);
    c_multi->cd(4)->SetGrid();
    cloneHs->Draw("hist nostack");
    c_multi->Write();

  }
  noCorr->Write();
  centralCorr->Write();
  fullCorr->Write();
  // treeFile->Close();
  calibrationFile->Close();
  outputFile->Close();
  return 0;
}
