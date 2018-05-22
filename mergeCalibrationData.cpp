// compile with
// g++ -o ../build/mergeCalibrationData mergeCalibrationData.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer

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
#include "TFitResult.h"

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
// void read_directory(const std::string& name, std::vector<std::string> &v)
// {
//     DIR* dirp = opendir(name.c_str());
//     struct dirent * dp;
//     while ((dp = readdir(dirp)) != NULL) {
//         v.push_back(dp->d_name);
//     }
//     closedir(dirp);
// }




void usage()
{
  // std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
  //           << "\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
  //           << "\t\t" << "<output.root>                                      - output file name"   << std::endl
  //           << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
  //           // << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
  //           << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
  //           << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
  //           << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
  //           << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
  //           << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
  //           << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 70e-12"  << std::endl
  //           << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
  //           << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
  //           << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
  //           << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
  //           << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
  //           << "\t\t" << "--smooth <value>                                   - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
  //           << "\t\t" << "--fitPercMin <value>                               - time fit min is set to ((gauss fit mean) - fitPercMin*(gauss fit sigma))  - default = 5"  << std::endl
  //           << "\t\t" << "--fitPercMax <value>                               - time fit max is set to ((gauus fit mean) - fitPercMax*(gauss fit sigma))  - default = 6" << std::endl
  //           << "\t\t" << "--divs <value>                                     - n of divisions when looking for FWHM - default = 10000"  << std::endl
  //           << "\t\t" << "--bins <value>                                     - n of bins in summary CTR histograms - deafult 40"  << std::endl
  //           << "\t\t" << "--func <value>                                     - function for fitting (default = 0)"  << std::endl
  //           << "\t\t" << "                                                   - 0 = crystalball "  << std::endl
  //           << "\t\t" << "                                                   - 1 = gauss+exp "  << std::endl
  //           << "\t\t" << "--unbinned                                         - use also the unbinned method to calculate CTR - default = 0 (false)"  << std::endl
  //           << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }
  // std::string inputFileName = "";
  // std::string outputFileName = "";
  // std::string calibrationFileName = "";
  // // std::string coincidenceCalibrationFileName = "";
  //
  // bool simulation = false;
  // Float_t length = 15.0; //mm
  // Float_t doiFraction = 0.5;
  // Float_t tagFwhm = 70.0e-12; //s
  // Float_t rmsLow = 1.75;
  // Float_t rmsHigh = 1.75;
  // Float_t histoMin = -15e-9;//s
  // Float_t histoMax = 15e-9;//s
  // Float_t fitPercMin = 5;
  // Float_t fitPercMax = 6;
  // int divs       = 10000;
  // int histoBins = 500;
  // int smooth = 0; //
  // int bins = 40;
  // double minCTR = 100;
  // double maxCTR = 500;
  // int func = 0;
  // bool unbinned = false;
  //
  // // parse arguments
  // static struct option longOptions[] =
  // {
  // 		{ "input", required_argument, 0, 0 },
  //     { "output", required_argument, 0, 0 },
  //     { "calibration", required_argument, 0, 0 },
  //     { "simulation", no_argument, 0, 0 },
  //     { "length", required_argument, 0, 0 },
  //     { "doiFraction", required_argument, 0, 0 },
  //     // { "coincidence", required_argument, 0, 0 },
  //     { "tagFwhm", required_argument, 0, 0 },
  //     { "rmsLow", required_argument, 0, 0 },
  //     { "rmsHigh", required_argument, 0, 0 },
  //     { "histoMin", required_argument, 0, 0 },
  //     { "histoMax", required_argument, 0, 0 },
  //     { "histoBins", required_argument, 0, 0 },
  //     { "smooth", required_argument, 0, 0 },
  //     { "fitPercMin", required_argument, 0, 0 },
  //     { "fitPercMax", required_argument, 0, 0 },
  //     { "divs", required_argument, 0, 0 },
  //     { "bins", required_argument, 0, 0 },
  //     { "func", required_argument, 0, 0 },
  //     { "unbinned", no_argument, 0, 0 },
  // 		{ NULL, 0, 0, 0 }
  // };
  //
  // while(1) {
  // 	int optionIndex = 0;
  // 	int c = getopt_long(argc, argv, "i:o:c:", longOptions, &optionIndex);
  // 	if (c == -1) {
  // 		break;
  // 	}
  // 	if (c == 'i'){
  // 		inputFileName = (char *)optarg;
  //   }
  // 	else if (c == 'o'){
  //     outputFileName = (char *)optarg;
  //   }
  //   else if (c == 'c'){
  //     calibrationFileName = (char *)optarg;
  //   }
  // 	else if (c == 0 && optionIndex == 0){
  //     inputFileName = (char *)optarg;
  //   }
  //   else if (c == 0 && optionIndex == 1){
  //     outputFileName = (char *)optarg;
  //   }
  //   else if (c == 0 && optionIndex == 2){
  //     calibrationFileName = (char *)optarg;
  //   }
  //   else if (c == 0 && optionIndex == 3){
  //     std::cout << "Dataset from simulation " << std::endl;
  //     simulation = true;
  //   }
  //   else if (c == 0 && optionIndex == 4){
  //     length = atof((char *)optarg);;
  //   }
  //   else if (c == 0 && optionIndex == 5){
  //     doiFraction = atof((char *)optarg);;
  //   }
  //   // else if (c == 0 && optionIndex == 6){
  //   //   coincidenceCalibrationFileName = (char *)optarg;
  //   // }
  //   else if (c == 0 && optionIndex == 6){
  //     tagFwhm = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 7){
  //     rmsLow = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 8){
  //     rmsHigh = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 9){
  //     histoMin = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 10){
  //     histoMax = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 11){
  //     histoBins = atoi((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 12){
  //     smooth = atoi((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 13){
  //     fitPercMin = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 14){
  //     fitPercMax = atof((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 15){
  //     divs = atoi((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 16){
  //     bins = atoi((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 17){
  //     func = atoi((char *)optarg);
  //   }
  //   else if (c == 0 && optionIndex == 18){
  //     unbinned = true;
  //   }
  // 	else {
  //     std::cout	<< "Usage: " << argv[0] << std::endl;
  // 		usage();
  // 		return 1;
  // 	}
  // }

  // std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;

  // read file in dir
  // std::vector<std::string> v;
  // read_directory(".", v);
  // std::copy(v.begin(), v.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;

  for(unsigned int i = 1 ; i < argc ; i++)
  {
    listInputFiles.push_back(argv[i]);
  }

  for(unsigned int i = 0 ; i < listInputFiles.size() ; i++)
  {
    std::cout << listInputFiles[i] << std::endl;
  }

  // TH2F* histo2d;
  std::stringstream sname;



  TList *list2D = new TList;
  TList *list3D = new TList;
  // list->Add(h1);
  // list->Add(h2);
  // list->Add(h3);
  TH2F *spectrum2dModule;
  TH3F *spectrum3dModule;
  std::vector<TFile*> files;

  struct spectra_t
  {
    TH1F *spectrum;
    int pos;
  };

  std::vector<spectra_t> chargeSpectra;
  std::vector<spectra_t> chargeSpectraCorrected;
  std::vector<spectra_t> CTRSpectra;

  TCanvas *bigSpectra = new TCanvas("BigSpectra","BigSpectra",1200,1200);
  TCanvas *BigCTRbasic = new TCanvas("BigCTRbasic","BigCTRbasic",1200,1200);


  for(unsigned int iFile = 0 ; iFile < listInputFiles.size() ; iFile++)
  {
    TFile *_file0 = TFile::Open(listInputFiles[iFile].c_str());
    files.push_back(_file0);
  }
  // int counter = 1;
  // int counterCorr = 1;
  // int counterCTR = 1;
  for(unsigned int iFile = 0 ; iFile < listInputFiles.size() ; iFile++)
  {
    // std::cout << listInputFiles[i] << std::endl;

    files[iFile]->cd("Module 0.0");
    TCanvas *canvas = (TCanvas*) gDirectory->Get("Flood Histogram 2D");
    TH2F *histo = (TH2F*) canvas->GetPrimitive("Flood Histogram 2D - Module 0.0");
    list2D->Add(histo);
    if(iFile == 0)
    {
      spectrum2dModule = (TH2F*)histo->Clone();
    }
    TCanvas *canvas3d = (TCanvas*) gDirectory->Get("Flood Histogram 3D");
    TH3I *histo3d = (TH3I*) canvas3d->GetPrimitive("Flood Histogram 3D - Module Module 0.0");
    list3D->Add(histo3d);
    if(iFile == 0)
    {
      spectrum3dModule = (TH3F*)histo3d->Clone();
    }


    if(iFile == 0)
    {
      //get bigSpectra to find the number of ch

      TCanvas *big = (TCanvas*) gDirectory->Get("BigSpectra");
      TList* listBig = (TList*) big->GetListOfPrimitives();
      int nPads = listBig->GetEntries();
      bigSpectra->Divide((int) ceil(sqrt(nPads)),(int) ceil(sqrt(nPads)));
      BigCTRbasic->Divide(ceil((int) sqrt(nPads)),(int) ceil(sqrt(nPads)));

    }

    //look for spectra in directories
    TList *listModule = gDirectory->GetListOfKeys();
    int nKeysMod = listModule->GetEntries();
    std::vector<std::string> keysModName;
    // fill a vector with the leaves names
    std::string mppc_prefix("MPPC ");
    for(int i = 0 ; i < nKeysMod ; i++){
      keysModName.push_back(listModule->At(i)->GetName());
    }

    std::vector<std::string> MPPCfolders;
    for(unsigned int i = 0 ; i < keysModName.size() ; i++)
    {
      if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
      {
        MPPCfolders.push_back(keysModName[i]);
      }
    }
    // std::cout << MPPCfolders.size() << std::endl;
    for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
    {
      // std::cout << MPPCfolders[iMppc] << std::endl;
      gDirectory->cd(MPPCfolders[iMppc].c_str());
      TList *listMppc = gDirectory->GetListOfKeys();
      int nKeysMppc = listMppc->GetEntries();
      std::vector<std::string> keysMppcName;
      // fill a vector with the leaves names
      std::string crystal_prefix("Crystal ");
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
      }
      for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
      {
        //  std::cout << CrystalFolders[iCry] << std::endl;
        gDirectory->cd(CrystalFolders[iCry].c_str());
        TList *listCry = gDirectory->GetListOfKeys();
        int nKeysCry = listCry->GetEntries();
        std::vector<std::string> keysCryName;
        if(nKeysCry) //if directory not empty
        {
          for(int i = 0 ; i < nKeysCry ; i++){
            keysCryName.push_back(listCry->At(i)->GetName());
          }
          // std::string CalibName;
          // std::string CutName;
          // std::vector<std::string> cutgNames;
          std::string spectrum_prefix("Charge Spectrum");
          std::string spectrumCorr_prefix("Charge Spectrum Corrected");
          std::string basicCTR_prefix("Basic CTR");
          std::string bigspectra_prefix("bigCanvasPosition");

          int bigCanvasPosition = 0;

          for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
          {
            if(!keysCryName[i].compare(0,bigspectra_prefix.size(),bigspectra_prefix)) // find detector channel
            {
              //  std::cout << keysCryName[i] << std::endl;
              std::stringstream snameCh;
              snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
              //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             //  istringstream()
              bigCanvasPosition = atoi(snameCh.str().c_str());
              //  std::cout <<temp_crystal.detectorChannel << std::endl;
             //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                       //  << temp_crystal.detectorChannel << std::endl;
            }
          }
          for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
          {
            if(!keysCryName[i].compare(0,spectrum_prefix.size(),spectrum_prefix)) //find calibration graph
            {
              //  std::cout << keysCryName[i] << std::endl;
              TCanvas* Canvas = NULL;
              TH1F *spectrum = NULL;
              Canvas = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
              if(Canvas)
              {
                spectrum = (TH1F*) Canvas->GetPrimitive(keysCryName[i].c_str());
              }
              if(spectrum)
              {
                spectra_t temp_spect;
                temp_spect.spectrum = spectrum;
                temp_spect.pos = bigCanvasPosition;
                // counter++;
                chargeSpectra.push_back(temp_spect);

                // std::cout << spectrum->GetName() << std::endl;
              }
            }
            if(!keysCryName[i].compare(0,spectrumCorr_prefix.size(),spectrumCorr_prefix)) //find calibration graph
            {
              //  std::cout << keysCryName[i] << std::endl;
              TCanvas* Canvas = NULL;
              TH1F *spectrum = NULL;
              Canvas = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
              if(Canvas)
              {
                spectrum = (TH1F*) Canvas->GetPrimitive(keysCryName[i].c_str());
              }
              if(spectrum)
              {
                spectra_t temp_spect;
                temp_spect.spectrum = spectrum;
                temp_spect.pos = bigCanvasPosition;
                // counterCorr++;
                chargeSpectraCorrected.push_back(temp_spect);
              }
            }
            if(!keysCryName[i].compare(0,basicCTR_prefix.size(),basicCTR_prefix)) //find calibration graph
            {
              //  std::cout << keysCryName[i] << std::endl;
              TCanvas* Canvas = NULL;
              TH1F *spectrum = NULL;
              Canvas = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
              if(Canvas)
              {
                spectrum = (TH1F*) Canvas->GetPrimitive(keysCryName[i].c_str());
              }
              if(spectrum)
              {
                spectra_t temp_spect;
                temp_spect.spectrum = spectrum;
                temp_spect.pos = bigCanvasPosition;
                // counterCTR++;
                CTRSpectra.push_back(temp_spect);

                // std::cout << spectrum->GetName() << std::endl;
              }
            }
          }
        }
      }
    }
  }






  // TH2F *spectrum2dModule = (TH2F*)h1->Clone("h");
  spectrum2dModule->Reset();
  spectrum2dModule->Merge(list2D);
  TCanvas *canvas2D = new TCanvas("Flood Histogram 2D","Flood Histogram 2D",800,800);
  canvas2D->cd();
  spectrum2dModule->Draw("COLZ");

  spectrum3dModule->Reset();
  spectrum3dModule->Merge(list3D);
  TCanvas *canvas3D = new TCanvas("Flood Histogram 3D","Flood Histogram 3D",800,800);
  canvas3D->cd();
  spectrum3dModule->Draw();

  if(chargeSpectraCorrected.size() > 0)
  {
    for(unsigned int i = 0 ; i < chargeSpectraCorrected.size(); i++)
    {
      bigSpectra->cd(chargeSpectraCorrected[i].pos);
      chargeSpectraCorrected[i].spectrum->Draw();
    }
  }
  else
  {
    for(unsigned int i = 0 ; i < chargeSpectra.size(); i++)
    {
      bigSpectra->cd(chargeSpectra[i].pos);
      chargeSpectra[i].spectrum->Draw();
    }
  }


  for(unsigned int i = 0 ; i < CTRSpectra.size(); i++)
  {
    BigCTRbasic->cd(CTRSpectra[i].pos);
    CTRSpectra[i].spectrum->Draw();
  }



  TFile *output = new TFile("merge_output.root","RECREATE");
  output->cd();
  canvas2D->Write();
  canvas3D->Write();
  bigSpectra->Write();
  BigCTRbasic->Write();
  output->Close();

  // //prepare output text file
  // std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  // textFileName += ".txt";
  // // std::cout << textFileName << std::endl;
  //
  // std::ofstream textfile;
  // textfile.open (textFileName.c_str(),std::ofstream::out);



  return 0;
}
