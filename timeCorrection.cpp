// compile with
// g++ -o ../build/timeCorrection timeCorrection.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer


// g++ -o ../build/timeCorrection timeCorrection.cpp -pthread -std=c++11 -m64 -I/home/marco/Programs/Root/root-6.10.08/include -L/home/marco/Programs/Root/root-6.10.08/lib -lGui -lCore -lRIO -lNet -lHist -lGraf -lGraf3d -lGpad -lTree -lTreePlayer -lRint -lPostscript -lMatrix -lPhysics -lMathCore -lThread -lMultiProc -pthread -lm -ldl -rdynamic -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

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


struct Crystal_t
{
  int number;
  int detectorChannel;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TGraph *wz;
  TH1F *simpleCTR;
  TH1F *centralCTR;
  TH1F *allCTR;
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

// class of real z positions of tag points (after check of alignment)
// class inputZ_t
// {
// public:
//   int pointsFromDoi;
//   int i;
//   int j;
//   std::vector<double> z;
//   inputZ_t(int a){ pointsFromDoi = a;};
//   void clear(){z.clear();};
//
//   friend std::istream& operator>>(std::istream& input, inputZ_t& s)
//   {
//     input >> s.i; //read i
//     input >> s.j;           //read
//     for(int p = 0; p < s.pointsFromDoi; p++)
//     {
//       double zValue;
//       input >> zValue;
//       s.z.push_back(zValue);
//     }
//     return input;
//   }
//
// };

// // define a function with 3 parameters
// Double_t crystalball(Double_t *x,Double_t *par) {
//       Double_t arg = 0;
//       if (par[2]!=0) arg = (x[0] - par[1])/par[2];
//       Double_t fitval = par[0]*TMath::Exp(-0.5*arg*arg);
//       return fitval;
// }


// void ComputeWithCrystalBall(TH1F* histo,float *res)
// {
//   // res[0] is fwhm
//   // res[1] is fwtm
//   double f1min = -10e-9;
//   double f1max = -8e-9;
//   TF1* f1  = new TF1("f1","crystalball",f1min,f1max);
//   f1->SetParameters(280,histo->GetMean(),histo->GetRMS(),1,3);
//   histo->Fit(f1,"","",-9.4e-9,-8.42-9);
//   double min = 1;
//   double max = -1;
//   double min10 = 1;
//   double max10 = -1;
//
//   int divs = 1000;
//   double step = (f1max-f1min)/divs;
//   for(int i = 0 ; i < divs ; i++)
//   {
//     if( (f1->Eval(f1min + i*step) < f1->GetMaximum()/2.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetMaximum()/2.0) )
//     {
//       min = f1min + (i+0.5)*step;
//     }
//     if( (f1->Eval(f1min + i*step) > f1->GetMaximum()/2.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetMaximum()/2.0) )
//     {
//       max = f1min + (i+0.5)*step;
//     }
//     if( (f1->Eval(f1min + i*step) < f1->GetMaximum()/10.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetMaximum()/10.0) )
//     {
//       min10 = f1min + (i+0.5)*step;
//     }
//     if( (f1->Eval(f1min + i*step) > f1->GetMaximum()/10.0) && (f1->Eval(f1min + (i+1)*step) < f1->GetMaximum()/10.0) )
//     {
//       max10 = f1min + (i+0.5)*step;
//     }
//   }
//
//   res[0] = 2.355*sqrt(2)*sqrt(pow((max-min)/2.355,2)-pow(70e-12/2.355,2));
//   res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((70e-12/2.355)*4.29,2));
// }

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
  f1->SetParameters(histo->GetMaximum(),histo->GetMean(),histo->GetRMS(),1,3);
  histo->Fit(f1,"Q","",fitMin,fitMax);
  double min,max,min10,max10;

  // int divs = 3000;
  double step = (f1max-f1min)/divs;
  // is [0] the max of the function???
  for(int i = 0 ; i < divs ; i++)
  {
    // std::cout << f1->Eval(f1min + i*step)  << "\t"
              // << f1->Eval(f1min + (i+1)*step) <<  "\t"
              // << f1->GetParameter(0) << "\t"
              // << std::endl;
    if( (f1->Eval(f1min + i*step) < f1->GetParameter(0)/2.0) && (f1->Eval(f1min + (i+1)*step) > f1->GetParameter(0)/2.0) )
    {

      min = f1min + (i+0.5)*step;
      // std::cout << f1->Eval(f1min + i*step) << " "
      //           << f1->Eval(f1min + (i+1)*step) << " "
      //           << f1->GetMaximum()/2.0 << " "
      //           << min  << std::endl;
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
  // std::cout
  // << min << " "
  // << max << " "
  // << max - min << " "
  // << std::endl << "CTR FWHM [ps] = "<< 2.355*sqrt(2)*sqrt(pow((max-min)/2.355,2)-pow(70e-12/2.355,2)) << std::endl
  // << max10 - min10 << std::endl
  // << "CTR FWTM [ps] = "<< sqrt(2)*sqrt(pow((max10-min10),2)-pow((70e-12/2.355)*4.29,2))
  // << std::endl;

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
  std::cout << "\t\t" << "[-i|--input] <temp.root>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
            << "\t\t" << "<temp.root>                                        - complete dataset (analysis ttree) of run"   << std::endl
            << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
            << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
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
  std::string coincidenceCalibrationFileName = "";

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
      { "coincidence", required_argument, 0, 0 },
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
    else if (c == 0 && optionIndex == 6){
      coincidenceCalibrationFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 7){
      tagFwhm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      rmsLow = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      rmsHigh = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      smooth = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      fitMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      fitMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 16){
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

  std::cout << "length * doiFraction = " << length * doiFraction << std::endl;



  TFile *treeFile = new TFile(inputFileName.c_str());
  TTree* tree = (TTree*) treeFile->Get("adc");
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

  std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
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
  }




  for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
  {
    // std::cout << MPPCfolders[iMppc] << std::endl;
    gDirectory->cd(MPPCfolders[iMppc].c_str());
    TList *listMppc = gDirectory->GetListOfKeys();
    int nKeysMppc = listMppc->GetEntries();
    std::vector<std::string> keysMppcName;
    // fill a vector with the leaves names
    std::string crystal_prefix("Crystal");
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
    }
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
         std::string channel_prefix("ChNum");
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


  // COINCIDENCE CALIBRATION FILE
  //
  TFile *coincidenceCalibrationFile = new TFile(coincidenceCalibrationFileName.c_str());
  coincidenceCalibrationFile->cd("Module 0.0");
  TList *listModuleCoinc = gDirectory->GetListOfKeys();
  int nKeysModCoinc = listModuleCoinc->GetEntries();
  std::vector<std::string> keysModNameCoinc;
  // fill a vector with the leaves names
  // std::string mppc_prefixCoinc("MPPC");
  for(int i = 0 ; i < nKeysModCoinc ; i++){
    keysModNameCoinc.push_back(listModuleCoinc->At(i)->GetName());
  }

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

  std::stringstream sformulaname;
  sformulaname << "FormulaTag";
  TCut taggingPhotopeakCutName;
  taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);

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

  Float_t FloodX,FloodY,FloodZ;
  Float_t ZPosition;
  TBranch      *bFloodX;
  TBranch      *bFloodY;
  TBranch      *bFloodZ;
  TBranch      *bZPosition;
  int inputChannels = 32;
  Float_t      *charge;
  Float_t      Tagging;
  Float_t      TaggingTimeStamp;
  TBranch      *bTagging;
  TBranch      *bTaggingTimeStamp;
  TBranch      **bCharge;
  Float_t      *timeStamp;
  TBranch      **btimeStamp;
  charge = new Float_t[inputChannels];
  timeStamp = new Float_t[inputChannels];
  bCharge = new TBranch*[inputChannels];
  btimeStamp = new TBranch*[inputChannels];
  int triggerChannel;
  TBranch      *bTriggerChannel;

  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    std::stringstream sname,stype;
    sname << "ch" << detector_channels[i];
    // stype << "ch" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    stype.str("");

    sname << "t" << detector_channels[i];
    // stype << "t" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
    stype.str("");
  }
  tree->SetBranchAddress("TriggerChannel",&triggerChannel,&bTriggerChannel);
  tree->SetBranchAddress("FloodX", &FloodX, &bFloodX);
  tree->SetBranchAddress("FloodY", &FloodY, &bFloodY);
  tree->SetBranchAddress("FloodZ", &FloodZ, &bFloodZ);
  tree->SetBranchAddress("ZPosition", &ZPosition, &bZPosition);
  tree->SetBranchAddress("Tagging", &Tagging, &bTagging);
  tree->SetBranchAddress("TaggingTimeStamp", &TaggingTimeStamp, &bTaggingTimeStamp);

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
            Float_t centralcorrection = 0.0;
            Float_t zeroCorrection    = 0.0;
            //no corr
            crystal[iCry].simpleCTR->Fill(timeStamp[crystal[iCry].detectorChannel] - TaggingTimeStamp);

            if(crystal[iCry].tw_correction)
            {
              // central corr
              // std::string deltaWGraph_prefix = "DeltaW Graph";
              // std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
              std::string graph_delay_prefix = "Graph Delay ch_";
              std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
              Float_t averageTimeStamp = 0.0;
              Float_t totalWeight = 0.0;
              // averageTimeStamp += timeStamp[crystal[iCry].detectorChannel];
              centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
              crystal[iCry].centralCTR->Fill((timeStamp[crystal[iCry].detectorChannel] + (centralcorrection)) - TaggingTimeStamp);


              if(crystal[iCry].delay.size())
              {
                //full corr
                zeroCorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*0)) - crystal[iCry].tw_correction->Eval(FloodZ);

                Float_t weight = 0.0;
                weight = pow(sqrt(pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),2)+pow(crystal[iCry].rms_tw_correction->Eval(crystal[iCry].wz->Eval(length*0)),2)),-2);

                averageTimeStamp += weight*(timeStamp[crystal[iCry].detectorChannel]- zeroCorrection);
                totalWeight += weight;
                // std::cout << i << " " << iCry << " " <<  crystal[iCry].number << "\n";
                for(unsigned int iGraph = 0; iGraph < crystal[iCry].delay.size();iGraph++)
                {
                  std::string graphName = crystal[iCry].delay[iGraph]->GetName();
                  int graphCh = atoi( graphName.substr( graph_delay_prefix.size(), graphName.size() - graph_delay_prefix.size() ).c_str() );
                  // std::cout << graphCh  << "\t";

                  std::string rmsName = crystal[iCry].rms_delay[iGraph]->GetName();
                  int rmsCh = atoi( rmsName.substr( rms_graph_delay_prefix.size(), rmsName.size() - rms_graph_delay_prefix.size() ).c_str() );
                  // std::cout << rmsCh  << "\t";
                  weight = pow(crystal[iCry].rms_delay[iGraph]->Eval(FloodZ),-2);
                  totalWeight += weight;
                  averageTimeStamp += weight*(timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ));
                  // std::cout << timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ) << "\t";
                }
                averageTimeStamp = averageTimeStamp/totalWeight;

                crystal[iCry].allCTR->Fill(averageTimeStamp + centralcorrection  - TaggingTimeStamp);
              }
            }
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
      extractCTR(crystal[iCry].simpleCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].simpleCTR->GetName() << "\t\t\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realBasicCTRfwhm = ret[0]*1e12;
      realBasicCTRfwtm = ret[1]*1e12;
      noCorr->Fill(ret[0]*1e12);
      crystal[iCry].simpleCTR->SetFillStyle(3001);
      crystal[iCry].simpleCTR->SetFillColor(kGreen);
      crystal[iCry].simpleCTR->SetLineColor(kGreen);
      crystal[iCry].simpleCTR->SetStats(0);
      crystal[iCry].simpleCTR->Write();
      crystal[iCry].simpleCTR->Scale(1.0/crystal[iCry].simpleCTR->GetMaximum());
    }

    if(crystal[iCry].centralCTR)
    {
      crystal[iCry].centralCTR->GetXaxis()->SetTitle("Time [s]");
      extractCTR(crystal[iCry].centralCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].centralCTR->GetName() << "\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realCentralCTRfwhm = ret[0]*1e12;
      realCentralCTRfwtm = ret[1]*1e12;
      centralCorr->Fill(ret[0]*1e12);
      crystal[iCry].centralCTR->SetFillStyle(3001);
      crystal[iCry].centralCTR->SetFillColor(kBlue);
      crystal[iCry].centralCTR->SetLineColor(kBlue);
      crystal[iCry].centralCTR->SetStats(0);
      crystal[iCry].centralCTR->Write();
      crystal[iCry].centralCTR->Scale(1.0/crystal[iCry].centralCTR->GetMaximum());
    }

    if(crystal[iCry].allCTR)
    {
      crystal[iCry].allCTR->GetXaxis()->SetTitle("Time [s]");
      extractCTR(crystal[iCry].allCTR,fitMin,fitMax,divs,ret);
      std::cout << crystal[iCry].allCTR->GetName() << "\t\t";
      std::cout << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << std::endl;
      realAllCTRfwhm = ret[0]*1e12;
      realAllCTRfwtm = ret[1]*1e12;
      fullCorr->Fill(ret[0]*1e12);
      crystal[iCry].allCTR->SetFillStyle(3001);
      crystal[iCry].allCTR->SetFillColor(kRed);
      crystal[iCry].allCTR->SetLineColor(kRed);
      crystal[iCry].allCTR->SetStats(0);
      crystal[iCry].allCTR->Write();
      crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());
    }

    sname.str("");

    sname << "Summary - Crystal " << crystal[iCry].number;
    TCanvas* c_summary = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    c_summary->cd();
    THStack *hs = new THStack("hs","");
    hs->Add(crystal[iCry].simpleCTR);
    hs->Add(crystal[iCry].centralCTR);
    hs->Add(crystal[iCry].allCTR);


    // std::cout << "Crystal " << crystal[iCry].number << std::endl;
    // std::cout << "CTR FWHM [ps] " << std::endl;
    hs->Draw("nostack");
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

    if(crystal[iCry].simpleCTR)
    {
      cloneBasic   = (TH1F*) crystal[iCry].simpleCTR->Clone();
      c_multi->cd(1);
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend1->AddEntry(cloneBasic,sname.str().c_str(),"f");
      cloneBasic->Draw();
      legend1->Draw();
    }
    if(crystal[iCry].centralCTR)
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
    if(crystal[iCry].allCTR)
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
    cloneHs->Draw("nostack");
    c_multi->Write();

  }
  noCorr->Write();
  centralCorr->Write();
  fullCorr->Write();
  treeFile->Close();
  calibrationFile->Close();
  outputFile->Close();
  return 0;
}
