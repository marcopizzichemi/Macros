// compile with
// g++ -o ../build/computeCTR computeCTR.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

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

}


void usage()
{
  std::cout << "\t\t" << "[-i|--input] <temp.root>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
            << "\t\t" << "<temp.root>                                        - complete dataset (analysis ttree) of run"   << std::endl
            // << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            // << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
            // << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
            // << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
            // << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            // << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            // << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
            // << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
            // << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in ps - default = 70"  << std::endl
            // << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
            // << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
            // << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in ns - default = -15"  << std::endl
            // << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in ns - default = 15"  << std::endl
            // << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            //   << "\t\t" << "--smooth <value>                                 - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
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
  double fitMin = -9.4e-9;
  double fitMax = -8.1e-9;
  int divs = 10000;


  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "fitMin", required_argument, 0, 0 },
      { "fitMax", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      fitMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 3){
      fitMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      divs = atoi((char *)optarg);
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  // if() // maybe put a check on inputFileName and outputFileName not to be empty strings?

  //FEEDBACK TO USER
  std::cout << std::endl;
  std::cout << "ANALYSIS PARAMETERS: " << std::endl
            << "Input file    = "      << inputFileName  << std::endl
            << "Output file   = "      << outputFileName << std::endl
            << "fitMin        = "      << fitMin         << std::endl
            << "fitMax        = "      << fitMax         << std::endl
            << "divs          = "      << divs           << std::endl;
  std::cout << std::endl;



  TFile *_file0 = TFile::Open(inputFileName.c_str());
  _file0->cd();

  TList *listCry = _file0->GetListOfKeys();
  std::cout << "testtttttttt" << std::endl;
  int nKeysCry = listCry->GetEntries();
  std::cout << nKeysCry << std::endl;
  std::vector<std::string> keysCryName;
  std::vector<TCanvas*> canvas;
  std::vector<TH1F*> histograms;

  int bins = 40;
  double minCTR = 100;
  double maxCTR = 500;
  TH1F* noCorr = new TH1F("No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* centralCorr = new TH1F("Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* fullCorr = new TH1F("Full Correction","Full Correction",bins,minCTR,maxCTR);

  if(nKeysCry) //if directory not empty
  {

    for(int i = 0 ; i < nKeysCry ; i++){
      keysCryName.push_back(listCry->At(i)->GetName());
    }


    std::string basic_prefix("No correction");
    std::string central_prefix("Central correction");
    std::string all_prefix("Full correction");

    // std::cout << "BASIC CTRs --------------------" << std::endl;
    for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
    {
      if(!keysCryName[i].compare(0,basic_prefix.size(),basic_prefix)) //
      {

        TH1F* histo = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
        histo->GetXaxis()->SetTitle("Time [s]");
        double ret[2];
        extractCTR(histo,fitMin,fitMax,divs,ret);
        std::cout << histo->GetName() << "\t\t\t";
        std::cout << ret[0]*1e12 << "\t"
                  << ret[1]*1e12 << std::endl;
        noCorr->Fill(ret[0]*1e12);
        histograms.push_back(histo);
      }
    }
    // std::cout << "-------------------------------" << std::endl;
    std::cout << std::endl;

    // std::cout << "CENTRAL CTRs --------------------" << std::endl;
    for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
    {
      if(!keysCryName[i].compare(0,central_prefix.size(),central_prefix)) //
      {
        TH1F* histo = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
        histo->GetXaxis()->SetTitle("Time [s]");
        double ret[2];
        extractCTR(histo,fitMin,fitMax,divs,ret);
        std::cout << histo->GetName() << "\t";
        std::cout << ret[0]*1e12 << "\t"
                  << ret[1]*1e12 << std::endl;
        centralCorr->Fill(ret[0]*1e12);
        histograms.push_back(histo);
      }
    }
    // std::cout << "-------------------------------" << std::endl;
    std::cout << std::endl;

    // std::cout << "ALL CTRs --------------------" << std::endl;
    for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
    {
      if(!keysCryName[i].compare(0,all_prefix.size(),all_prefix)) //
      {
        TH1F* histo = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
        histo->GetXaxis()->SetTitle("Time [s]");
        double ret[2];
        extractCTR(histo,fitMin,fitMax,divs,ret);
        std::cout << histo->GetName() << "\t\t";
        std::cout << ret[0]*1e12 << "\t"
                  << ret[1]*1e12 << std::endl;
        fullCorr->Fill(ret[0]*1e12);
        histograms.push_back(histo);
      }
    }
    // std::cout << "-------------------------------" << std::endl;
    std::cout << std::endl;

  }


  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  noCorr->Write();
  centralCorr->Write();
  fullCorr->Write();
  for(unsigned int i = 0; i < histograms.size() ; i++)
  {
    histograms[i]->Write();
  }
  outputFile->Close();
  _file0->Close();
}
