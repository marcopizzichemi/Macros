// compile with
// g++ -o ../build/amplitudeWalk amplitudeWalk.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lSpectrum


// small program to extract amplitude correction impact
// procedure
// 1. process data with ModuleCalibration as usual, either all together or split by crystal or mppc
// 2. run amplitudeWalk on the dataset with the input, output and calibration file, specifying the crystal, its length and the number of divisions into which "slicing" the doi
// typical command
// amplitudeWalk -i TTree_ -o ampl34.root -c pOutput_B3.root --length 15 --divs 15 --crystal 34
// this assumes crystal 34 is on MPPC B3 (it will fail otherwise)


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


void demoDoi(TH2F* original_histoADCvsDt,TH2F* original_histo1_5,TCanvas* canvas,double min,double max, double length,double minADC,double maxADC,double minT,double maxT,int rebin_all,int rebin_slice,TCanvas* single_canvases[4])
{

  // TFile *inputFile = new TFile(fileName.c_str());
  // inputFile->cd();
  TH2F* histoADCvsDt = (TH2F*) original_histoADCvsDt->Clone();
  TH2F* histo1_5     = (TH2F*) original_histo1_5->Clone();
  // TH2F* histoADCvsDt = (TH2F*) gDirectory->Get("histoADCvsDt");
  std::stringstream sname;
  // sname << "histoSat_" << min << "_" << max;
  // TH2F* histo1_5     = (TH2F*) gDirectory->Get(sname.str().c_str());
  // sname.str("");

  TH2F* clone_histoADCvsDt = (TH2F*) histoADCvsDt->Clone();
  TH2F* clone_histo1_5     = (TH2F*) histo1_5->Clone();

  clone_histoADCvsDt->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histo1_5    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histoADCvsDt->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histo1_5    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");

  clone_histoADCvsDt->SetTitle("All crystal -> DOI + Amplitude effect");

  sname << "Fixed DOI = "<< length-min << "mm to "<< length-max<< "mm -> Only Amplitude effect";
  clone_histo1_5    ->SetTitle(sname.str().c_str());
  sname.str("");

  clone_histoADCvsDt->Rebin2D(rebin_all,rebin_all);
  clone_histo1_5->Rebin2D(rebin_all,rebin_all);

  histoADCvsDt->Rebin2D(rebin_slice,rebin_slice);
  histo1_5->Rebin2D(rebin_slice,rebin_slice);
  histoADCvsDt->ProfileX();
  histo1_5->ProfileX();

  sname << histoADCvsDt->GetName() << "_pfx";
  TProfile *histoADCvsDt_pfx = (TProfile*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  sname << histo1_5->GetName() << "_pfx";
  TProfile *histo1_5_pfx = (TProfile*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  TProfile* clone_histoADCvsDt_pfx = (TProfile*) histoADCvsDt_pfx->Clone();
  TProfile* clone_histo1_5_pfx     = (TProfile*) histo1_5_pfx->Clone();


  clone_histoADCvsDt_pfx->GetYaxis()->SetRangeUser(minT,maxT);
  clone_histo1_5_pfx    ->GetYaxis()->SetRangeUser(minT,maxT);
  clone_histoADCvsDt_pfx->SetLineColor(kBlue);
  clone_histo1_5_pfx    ->SetLineColor(kRed);
  clone_histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  TLegend *legend1 = new TLegend(0.54,0.62,0.89,0.89,"");
  legend1->SetFillStyle(0);
  legend1->AddEntry(clone_histoADCvsDt_pfx,"DOI + Amplitude","l");
  legend1->AddEntry(clone_histo1_5_pfx,"Only Amplitude","l");

  TLegend *legend2 = new TLegend(0.54,0.62,0.89,0.89,"");
  legend2->SetFillStyle(0);
  legend2->AddEntry(histoADCvsDt_pfx,"DOI + Amplitude","l");
  legend2->AddEntry(histo1_5_pfx,"Only Amplitude","l");

  histoADCvsDt_pfx->GetYaxis()->SetRangeUser(minT,maxT);
  histo1_5_pfx    ->GetYaxis()->SetRangeUser(minT,maxT);
  histoADCvsDt_pfx->SetLineColor(kBlue);
  histo1_5_pfx    ->SetLineColor(kRed);
  histoADCvsDt_pfx->GetXaxis()->SetRangeUser(minADC,maxADC);
  histo1_5_pfx->GetXaxis()    ->SetRangeUser(minADC,maxADC);
  histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  // fitting
  // TF1 *pol1_doi = new TF1("pol1_doi","pol1");
  // TF1 *pol1_ampl = new TF1("pol1_ampl","pol1");
  // pol1_doi ->SetLineColor(kBlue);
  // pol1_ampl->SetLineColor(kRed);
  // histoADCvsDt_pfx->Fit(pol1_doi);
  // histo1_5_pfx->Fit(pol1_ampl);

  // TCanvas* canvas = new TCanvas("canvas","canvas",1600,1300);
  canvas->Divide(2,2);
  canvas->cd(1);
  clone_histoADCvsDt->Draw("COLZ");
  single_canvases[0]->cd();
  clone_histoADCvsDt->Draw("COLZ");

  canvas->cd(2);
  clone_histo1_5->Draw("COLZ");
  single_canvases[1]->cd();
  clone_histo1_5->Draw("COLZ");

  canvas->cd(3);
  clone_histoADCvsDt_pfx->SetTitle("All Amplitude range");
  clone_histoADCvsDt_pfx->Draw();
  clone_histo1_5_pfx->Draw("same");
  legend1->Draw("same");
  single_canvases[2]->cd();
  clone_histoADCvsDt_pfx->SetTitle("All Amplitude range");
  clone_histoADCvsDt_pfx->Draw();
  clone_histo1_5_pfx->Draw("same");
  legend1->Draw("same");

  canvas->cd(4);
  histoADCvsDt_pfx->SetTitle("511 KeV Amplitude range");
  histoADCvsDt_pfx->Draw();
  histo1_5_pfx->Draw("same");
  legend2->Draw("same");
  single_canvases[3]->cd();
  histoADCvsDt_pfx->SetTitle("511 KeV Amplitude range");
  histoADCvsDt_pfx->Draw();
  histo1_5_pfx->Draw("same");
  legend2->Draw("same");
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
  TH1F *poliCorrCTR;
  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TH1F *poliCorrCTR_norm;
  TTreeFormula *Formula;
  TTreeFormula *FormulaPhotopeak;
  std::vector<double> z;
  TGraph* tw_correction;
  TGraph* rms_tw_correction;
  std::vector<TGraph*> delay;
  std::vector<TGraph*> rms_delay;
  const char* path;
  bool accepted;
  bool polishedCorrection;
  std::vector<int>    tChannelsForPolishedCorrection;
  std::vector<double> meanForPolishedCorrection;
  std::vector<double> fwhmForPolishedCorrection;
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

void extractCTR(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double tagFwhm, double* res)
{
  //first, dummy gaussian fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  histo->Fit(gaussDummy,"QN");

  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  // std::cout << f1min << " " << f1max << std::endl;
  TF1* f1  = new TF1("f1","crystalball");
  f1->SetLineColor(kBlack);
  f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }
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
  res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
  res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));
  delete cTemp;
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
            << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
            << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            << "\t\t" << "--divs <value>                                     - num of divisions of doi length"  << std::endl
            << "\t\t" << "--crystal <value>                                  - num of crystal invesitgated"  << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for spectra - default = 500"  << std::endl
            << "\t\t" << "--adcMinSat <value>                                - min of adc histograms - default = 50000"  << std::endl
            << "\t\t" << "--adcMaxSat <value>                                - max of adc histograms - default = 90000"  << std::endl
            << "\t\t" << "--adcMinNoSat <value>                              - min of adc histograms - default = 34000"  << std::endl
            << "\t\t" << "--adcMaxNoSat <value>                              - max of adc histograms - default = 46000"  << std::endl
            << "\t\t" << "--rebin_all <value>                                - rebin value for all ADCvsCTR spectra - default = 10"  << std::endl
            << "\t\t" << "--rebin_slice <value>                              - rebin value for slice ADCvsCTR spectra - default = 30"  << std::endl




            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  gStyle->SetOptStat(0);
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

  int divs       = 10000;
  int histoBins = 500;
  int steps = 15;
  int crystalNum = 0;
  double length = 15.0; //mm
  int rebin_all = 10;
  int rebin_slice = 30;

  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s

  double adcMinSat   = 50000;
  double adcMaxSat   = 90000;
  double adcMinNoSat = 34000;
  double adcMaxNoSat = 46000;

  // int smooth = 0; //
  // int bins = 40;
  // double minCTR = 100;
  // double maxCTR = 500;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "simulation", no_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "crystal", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "adcMinSat", required_argument, 0, 0 },
      { "adcMaxSat", required_argument, 0, 0 },
      { "adcMinNoSat", required_argument, 0, 0 },
      { "adcMaxNoSat", required_argument, 0, 0 },
      { "rebin_all", required_argument, 0, 0 },
      { "rebin_slice", required_argument, 0, 0 },
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
      steps = atoi((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 6){
      crystalNum = atoi((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 7){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      adcMinSat = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      adcMaxSat = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      adcMinNoSat = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      adcMaxNoSat = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      rebin_all = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      rebin_slice = atoi((char *)optarg);
    }

		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  // std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;

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


  //prepare output text file
  // std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  // textFileName += ".txt";
  // // std::cout << textFileName << std::endl;
  //
  // std::ofstream textfile;
  // textfile.open (textFileName.c_str(),std::ofstream::out);




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
       temp_crystal.poliCorrCTR = NULL;
       temp_crystal.simpleCTR_norm = NULL;
       temp_crystal.centralCTR_norm = NULL;
       temp_crystal.allCTR_norm = NULL;
       temp_crystal.poliCorrCTR_norm = NULL;
       temp_crystal.wz = NULL;
       temp_crystal.accepted = true;
       temp_crystal.tw_correction = NULL;
       temp_crystal.polishedCorrection = false;


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
         std::string t_channels_for_poli_prefix("tChannelsForPolishedCorrection");
         std::string mean_for_poli_prefix("meanForPolishedCorrection");
         std::string fwhm_for_poli_prefix("fwhmForPolishedCorrection");
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

           if(!keysCryName[i].compare(0,t_channels_for_poli_prefix.size(),t_channels_for_poli_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrection",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrection = v[0];
           }

           if(!keysCryName[i].compare(0,mean_for_poli_prefix.size(),mean_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("meanForPolishedCorrection",v);
             temp_crystal.meanForPolishedCorrection = v[0];
           }

           if(!keysCryName[i].compare(0,fwhm_for_poli_prefix.size(),fwhm_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("fwhmForPolishedCorrection",v);
             temp_crystal.fwhmForPolishedCorrection = v[0];
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
         sname.str("");
         sname << "Polished correction - Crystal " << temp_crystal.number;
         temp_crystal.poliCorrCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
         sname.str("");
         TCut globalCut; // the cut for the formula
         globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         // globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }
         sname.str("");

         sname << "Formula" << temp_crystal.number;
         TTreeFormula* Formula = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         formulas.Add(Formula);
         temp_crystal.Formula = Formula;
         sname.str("");

         TCut photopeakCutFormula;
         photopeakCutFormula += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         sname << "FormulaPhotopeak" << temp_crystal.number;
         TTreeFormula* FormulaPhotopeak = new TTreeFormula(sname.str().c_str(),photopeakCutFormula,tree);
         formulas.Add(FormulaPhotopeak);
         temp_crystal.FormulaPhotopeak = FormulaPhotopeak;
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

  bool foundCrystal = false;
  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
      if(crystalNum == crystal[i].number)
      {
        foundCrystal = true;
      }
    }
  }

  if(!foundCrystal)
  {
    std::cout << "Chosen crystal (number " << crystalNum << ") not found in calibration data!! Aborting" << std::endl;
    return 1;
  }

  //notify TTreeFormula(s) to TChain
  tree->SetNotify(&formulas);

  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  std::cout << "Total number of events = " << nevent << std::endl;
  long int goodEvents = 0;
  long int counter = 0;

  std::vector<TH2F*> histoNoSat;
  // = new TH2F("histoNoSat","histoNoSat",1000,0,65000,1000,-9.5e-9,-7.5e-9);
  std::vector<TH2F*> histoSat;
  // = new TH2F("histoSat","histoSat",1000,0,100000,1000,-9.5e-9,-7.5e-9);

  TH2F *histoADCvsW = new TH2F("histoADCvsW","histoADCvsW",1000,adcMinSat,adcMaxSat,1000,0,1);
  TH2F *histoADCvsDt = new TH2F("histoADCvsDt","histoADCvsDt",1000,adcMinSat,adcMaxSat,1000,histoMin,histoMax);
  TH2F *histoADCvsDtNoSat = new TH2F("histoADCvsDtNoSat","histoADCvsDtNoSat",1000,adcMinNoSat,adcMaxNoSat,1000,histoMin,histoMax);

  TH1F* basicCTR             = new TH1F("basicCTR","basicCTR",histoBins,histoMin,histoMax) ;
  TH1F* singleChargeSpectrum_NoSat = new TH1F("singleChargeSpectrum_NoSat","singleChargeSpectrum_NoSat",histoBins,0,adcMaxNoSat) ;
  TH1F* singleChargeSpectrum_Sat = new TH1F("singleChargeSpectrum_Sat","singleChargeSpectrum_Sat",histoBins,9,adcMaxSat) ;

  // TH1F *standardCTR = new TH1F("standardCTR","standardCTR",1000,-9.5e-9,-7.5e-9);
  // TH1F *amplCorrCTR = new TH1F("amplCorrCTR","amplCorrCTR",1000,-9.5e-9,-7.5e-9);

  // TH1F *standardCTR_w = new TH1F("standardCTR_w","standardCTR_w",1000,-9.5e-9,-7.5e-9);
  // TH1F *amplCorrCTR_w = new TH1F("amplCorrCTR_w","amplCorrCTR_w",1000,-9.5e-9,-7.5e-9);

  // TH1F *standardCTR_corr2 = new TH1F("standardCTR_corr2","standardCTR_corr2",1000,-9.5e-9,-7.5e-9);
  // TH1F *amplCorrCTR_corr2 = new TH1F("amplCorrCTR_corr2","amplCorrCTR_corr2",1000,-9.5e-9,-7.5e-9);

  // int steps = 5;
  // length = 15
  for(int iStep = 0; iStep < steps; iStep++)
  {
    std::stringstream sname;
    sname << "histoSat_" << (length/steps)*iStep << "_" << (length/steps)*(iStep+1);
    TH2F* tempHisto = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,adcMinSat,adcMaxSat,1000,histoMin,histoMax);
    histoSat.push_back(tempHisto);

    sname.str("");
    sname << "histoNoSat_" << (length/steps)*iStep << "_" << (length/steps)*(iStep+1);
    TH2F* tempHisto2 = new TH2F(sname.str().c_str(),sname.str().c_str(),1000,adcMinNoSat,adcMaxNoSat,1000,histoMin,histoMax);
    histoNoSat.push_back(tempHisto2);


  }

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
          if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal without photopeak
          {
            if(crystal[iCry].number == crystalNum)  //temp
            {
              if(crystal[iCry].tw_correction)
              {
                goodEvents++;
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

                histoADCvsW->Fill(centralChargeCorr,FloodZ);
                double ctr = timeStamp[crystal[iCry].timingChannel] - timeStamp[taggingCrystalTimingChannel];


                basicCTR->Fill(ctr);
                singleChargeSpectrum_NoSat ->Fill(centralChargeOriginal);
                singleChargeSpectrum_Sat   ->Fill(centralChargeCorr);

                // find wlimits
                // choose 1 mm steps


                for(int iW = 0; iW < steps; iW++)
                {
                  if( (FloodZ <=   (crystal[iCry].wz->Eval( ((length/steps)*iW) )))  &&  (FloodZ >=   (crystal[iCry].wz->Eval( ((length/steps)*(iW+1)) )))  )
                  {
                    histoSat[iW]->Fill(centralChargeCorr,ctr);
                    histoNoSat[iW]->Fill(centralChargeOriginal,ctr);
                  }

                }

                histoADCvsDt->Fill(centralChargeCorr,ctr);
                histoADCvsDtNoSat->Fill(centralChargeOriginal,ctr);
              }
            }



            // //temp commented
            // Float_t centralcorrection = 0.0;
            // Float_t zeroCorrection    = 0.0;
            // //no corr
            // crystal[iCry].simpleCTR->Fill(timeStamp[crystal[iCry].timingChannel] -
            //                               timeStamp[taggingCrystalTimingChannel]);
            //
            // if(crystal[iCry].tw_correction)
            // {
            //
            //   //calculate FloodZ...
            //   Float_t FloodZ;
            //   float centralChargeOriginal;
            //   float centralSaturation;
            //   float centralPedestal;
            //   Float_t division = 0.0;
            //
            //   centralChargeOriginal = charge[crystal[iCry].detectorChannel];
            //   for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
            //   {
            //     if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
            //     {
            //       centralSaturation = detectorSaturation[iSat].saturation;
            //       centralPedestal = detectorSaturation[iSat].pedestal;
            //     }
            //   }
            //   float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
            //
            //   for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
            //   {
            //     float originalCh = charge[crystal[iCry].relevantForW[iW]];
            //     float saturationCh;
            //     float pedestalCorr;
            //     for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
            //     {
            //       if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
            //       {
            //         saturationCh = detectorSaturation[iSat].saturation;
            //         pedestalCorr = detectorSaturation[iSat].pedestal;
            //       }
            //     }
            //     division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
            //   }
            //
            //   FloodZ = centralChargeCorr / division;
            //
            //   // central corr
            //   // std::string deltaWGraph_prefix = "DeltaW Graph";
            //   // std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
            //
            //
            //   Float_t averageTimeStamp = 0.0;
            //   Float_t totalWeight = 0.0;
            //   // averageTimeStamp += timeStamp[crystal[iCry].detectorChannel];
            //   centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
            //   crystal[iCry].centralCTR->Fill((timeStamp[crystal[iCry].timingChannel] + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel]);
            //
            //
            //   if(crystal[iCry].delay.size())
            //   {
            //     // full corr v0
            //     // first evalutate all ts as if they where read by the central crystal, i.e. correct all the neighboring channels
            //     // using the delta
            //     //
            //     Float_t weight = 0.0;
            //     weight = pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),-2);
            //     averageTimeStamp += weight * timeStamp[crystal[iCry].timingChannel];
            //     totalWeight += weight;
            //     //
            //     //full corr
            //     // zeroCorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
            //     //
            //     // Float_t weight = 0.0;
            //     // weight = pow(sqrt(pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),2)+pow(crystal[iCry].rms_tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)),2)),-2);
            //     //
            //     // averageTimeStamp += weight*(timeStamp[crystal[iCry].timingChannel]- zeroCorrection);
            //     // totalWeight += weight;
            //     // // std::cout << i << " " << iCry << " " <<  crystal[iCry].number << "\n";
            //     for(unsigned int iGraph = 0; iGraph < crystal[iCry].delay.size();iGraph++)
            //     {
            //
            //       std::string graphName = crystal[iCry].delay[iGraph]->GetName();
            //       std::string rmsName = crystal[iCry].rms_delay[iGraph]->GetName();
            //
            //       std::size_t foundGraph = graphName.find_last_of("_");
            //       std::size_t foundRms   = rmsName.find_last_of("_");
            //
            //       std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
            //       std::string tChannelStringFromRms   = rmsName.substr(foundRms  +1);
            //
            //       // std::stringstream sname;
            //       // sname << "Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
            //       // std::string graph_delay_prefix = sname.str();
            //       // std::cout << tChannelStringFromGraph << std::endl;
            //       // std::cout << tChannelStringFromRms << std::endl;
            //
            //       // sname.str("");
            //       // sname << "RMS Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
            //       // std::string rms_graph_delay_prefix = sname.str();
            //       // std::cout << tChannelStringFromRms << std::endl;
            //       // sname.str("");
            //
            //
            //
            //
            //
            //
            //       int graphCh = atoi(tChannelStringFromGraph.c_str() );
            //       // std::cout << graphCh  << "\n";
            //
            //
            //       int rmsCh   = atoi( tChannelStringFromRms.c_str() );
            //       if(graphCh != rmsCh)
            //       {
            //         std::cout << "ERROR! TGraphs of delay and rms are from different timing channels!!!!" << std::endl;
            //         break;
            //       }
            //       // std::cout << rmsCh  << "\n";
            //       weight = pow(crystal[iCry].rms_delay[iGraph]->Eval(FloodZ),-2);
            //       totalWeight += weight;
            //       averageTimeStamp += weight*(timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ));
            //       // std::cout << timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ) << "\t";
            //     }
            //     averageTimeStamp = averageTimeStamp/totalWeight;
            //     // then correct the average of ts for central correction
            //     crystal[iCry].allCTR->Fill((averageTimeStamp + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel]);
            //     // crystal[iCry].allCTR->Fill(averageTimeStamp + zeroCorrection  - timeStamp[taggingCrystalTimingChannel]);
            //   }
            // }
            //
            // if(crystal[iCry].polishedCorrection)
            // {
            //   //central time stamp
            //   // std::cout << "----------------" << std::endl;
            //   float weight = 0.0;
            //   float meanTimeStamp = 0.0;
            //   float sumWeight = 0.0;
            //   // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;
            //
            //   weight = pow(crystal[iCry].fwhmForPolishedCorrection[0],-2);
            //   float t_0 = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[0]];
            //
            //   meanTimeStamp += weight * t_0;
            //   sumWeight += weight;
            //   // std::cout << weight << " " << t_0 << " " << meanTimeStamp << " " << sumWeight << std::endl;
            //
            //   // std::cout << t_0 << "\t";
            //
            //   for(unsigned int iPoli = 1; iPoli < crystal[iCry].tChannelsForPolishedCorrection.size(); iPoli++)
            //   {
            //     // std::cout << iPoli << std::endl;
            //     // std::cout << "timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] = " << timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] << std::endl;
            //     // std::cout << "crystal[iCry].meanForPolishedCorrection[iPoli] = " << crystal[iCry].meanForPolishedCorrection[iPoli] << std::endl;
            //
            //     float t_i = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] - crystal[iCry].meanForPolishedCorrection[iPoli];
            //     float weight_i = pow(crystal[iCry].fwhmForPolishedCorrection[iPoli],-2);
            //     meanTimeStamp += weight_i * t_i;
            //     sumWeight += weight_i;
            //     // std::cout << weight_i << " " << t_i << " " << meanTimeStamp << " " << sumWeight << std::endl;
            //     // std::cout << t_i << "\t";
            //   }
            //   // std::cout << std::endl;
            //   meanTimeStamp = meanTimeStamp / sumWeight;
            //   // std::cout << std::endl
            //   //           << meanTimeStamp << "\t"
            //   //           << std::endl;
            //   crystal[iCry].poliCorrCTR->Fill(meanTimeStamp - timeStamp[taggingCrystalTimingChannel]);
            //
            // }
            //
            //
            //
            // // end of temp commented
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
  // std::cout << "Good events = " << goodEvents << std::endl;
  std::sort(crystal.begin(), crystal.end(), compareByNumber);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  // histoNoSat->Write();

  //look for limits
  // double tagFitLowerFraction = 1.8;
  // double tagFitUpperFraction = 2.0;
  double minT = basicCTR->GetMean() - 3.0*basicCTR->GetRMS();
  double maxT = basicCTR->GetMean() + 3.0*basicCTR->GetRMS();



  for(unsigned int i = 0 ; i < histoSat.size(); i++)
  {
    std::stringstream canvas_name;
    canvas_name << "With saturation correction - " << (length/steps)*i << "_" << (length/steps)*(i+1);
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(),canvas_name.str().c_str(),1600,1300);
    TCanvas* single_canvases[4];
    for(int j = 0 ; j < 4 ; j++)
    {
      std::stringstream name2;
      name2 << canvas_name.str() <<"_Canvas_" << j+1 ;
      single_canvases[j] = new TCanvas(name2.str().c_str(),name2.str().c_str(),1200,800);
    }

    demoDoi(histoADCvsDt,histoSat[i],canvas,(length/steps)*i,(length/steps)*(i+1),length,adcMinSat,adcMaxSat,minT,maxT,rebin_all,rebin_slice,single_canvases);
    histoSat[i]->Write();
    canvas->Write();
    for(int j = 0 ; j < 4 ; j++)
    {
      single_canvases[j]->Write();
    }
  }



  for(unsigned int i = 0 ; i < histoNoSat.size(); i++)
  {
    std::stringstream canvas_name;
    canvas_name << "Without saturation corretion - " << (length/steps)*i << "_" << (length/steps)*(i+1);
    TCanvas *canvas = new TCanvas(canvas_name.str().c_str(),canvas_name.str().c_str(),1600,1300);
    TCanvas* single_canvases[4];
    for(int j = 0 ; j < 4 ; j++)
    {
      std::stringstream name2;
      name2 << canvas_name.str() <<"_Canvas_" << j+1 ;
      single_canvases[j] = new TCanvas(name2.str().c_str(),name2.str().c_str(),1200,800);
    }
    demoDoi(histoADCvsDtNoSat,histoNoSat[i],canvas,(length/steps)*i,(length/steps)*(i+1),length,adcMinNoSat,adcMaxNoSat,minT,maxT,rebin_all,rebin_slice,single_canvases);
    histoNoSat[i]->Write();
    canvas->Write();
    for(int j = 0 ; j < 4 ; j++)
    {
      single_canvases[j]->Write();
    }
  }

  std::cout << std::endl;
  std::cout << "Crystal " << crystalNum << " = " << goodEvents << " events" << std::endl;

  //write data
  basicCTR->Write();
  singleChargeSpectrum_Sat->Write();
  singleChargeSpectrum_NoSat->Write();

  histoADCvsW->Write();
  histoADCvsDt->Write();
  histoADCvsDtNoSat->Write();


  calibrationFile->Close();
  outputFile->Close();

  return 0;
}
