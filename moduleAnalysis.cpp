// Program to analyze the performance of a single module

// compile with
// g++ -o ../build/moduleAnalysis moduleAnalysis.cpp `root-config --cflags --glibs` 

//HOWTO
//This program analyzes the output of a moduleCalibration run
//The generic syntax is

//moduleAnalysis moduleCalibration.root [calibration_params.txt] [pointsFromDoi]

// moduleCalibration.root   // name of the root file to analyze
// calibration_params.txt   // OPTIONAL - name of the file with results of DOI tagging bench scan. The scan can have as many points as you want, for as many crystals as you like. In any case
                            // only the non-edge modules will be considered and this data is used only to calcolate the sigma w in order to get a value for DOI res - if no value is given, DOI res is not analyzed
// pointsFromDoi            // OPTIONAL - the number of points in vertical DOI scan. By default is 1

// So in order to perform a complete analysis, the DOI scan has to be already performed and analyzed following the procedure described in the doiAnalysis repostory, file procedure.txt

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <getopt.h>
#include <iostream>
#include <set>
#include <assert.h>
#include <string.h>
#include <vector>
#include <algorithm>
#include <sstream>
#include <fstream>
#include <string>
#include <iomanip>      // std::setprecision
#include <TROOT.h>
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TPad.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


// class of input points from doi tag bench
class inputDoi_t
{
public:
  int i;
  int j;
  double m;
  double q;
  double doires;
  double avgs;
  std::vector<double> w;
  std::vector<double> sw;
  std::vector<double> sqrt_nentries;
  std::vector<double> z;
  std::vector<double> sz;  
  int pointsFromDoi;
  inputDoi_t(int a){ pointsFromDoi = a;};
  void clear()
  {
    w.clear();
    sw.clear();
    sqrt_nentries.clear();
    z.clear();
    sz.clear();
  };
  friend std::istream& operator>>(std::istream& input, inputDoi_t& s)
  {
    input >> s.i; 
    input >> s.j;           
    input >> s.m; 
    input >> s.q;    
    input >> s.doires; 
    input >> s.avgs;
    for(int p = 0; p < s.pointsFromDoi; p++)
    {
      double wValue,swValue,sqrtValue;
      input >> wValue >> swValue >> sqrtValue; 
      s.w.push_back(wValue);
      s.sw.push_back(swValue);
      s.sqrt_nentries.push_back(sqrtValue);
    }
    return input;
  }
};                   


//----------------//
//  MAIN PROGRAM  //
//----------------//
int main(int argc, char **argv)
{
  if(argc < 2) // check input from command line
  {
    std::cout << "USAGE:\t\t moduleAnalysis moduleCalibration.root [calibration_params.txt] [pointsFromDoi]" << std::endl;
    std::cout << std::endl;
    return 1;
  }
  
  
  //play with strings to extract the name 
  std::string rootFileName = argv[1];
  std::string rootFileNameNoExtension = rootFileName.substr(0,rootFileName.length() -5 );
  
  //file with spectra from module calibration
  TFile *f = new TFile(rootFileName.c_str());
  f->cd("Module 0.0");
  
  
  //file with doi tag m, q and points per crystal
  //if it is given, set a flag to calculate the doi precision
  bool doDoiPrecision = false;
  std::ifstream fDoiTag;
  if(argc > 2)
  {
    doDoiPrecision = true;
    fDoiTag.open(argv[2],std::ios::in);
  }
  int pointsFromDoi = 1;
  if(argc > 3)
  {
    pointsFromDoi = atoi(argv[3]);
  }
  
  
  TString outputFileName = "output_" + rootFileNameNoExtension + ".root";
//   if(argc > 3)
//     outputFileName = argv[3];
  TFile *fOut = new TFile(outputFileName,"RECREATE");
  
  //useful variables...
  int nmodulex = 1;
  int nmoduley = 1;
  int nmppcx = 4;
  int nmppcy = 4;
  int ncrystalsx = 2;
  int ncrystalsy = 2;
  std::string letter[4] = {"A","B","C","D"};
  std::string number[4] = {"1","2","3","4"};
  
  std::vector<inputDoi_t> inputDoi;
  inputDoi_t tempInputDoi(pointsFromDoi);
  while(fDoiTag >> tempInputDoi)
  {
    inputDoi.push_back(tempInputDoi);
    tempInputDoi.clear();
  }
  
  
  
// DEBUG
//   for(int i = 0 ; i < inputDoi.size() ; i++)
//   {
//     std::cout << inputDoi[i].i << " " << inputDoi[i].j << " " << inputDoi[i].m << " " << inputDoi[i].q << " ";
//     for(int k = 0 ; k < pointsFromDoi ; k++)
//     {
//       std::cout << inputDoi[i].w[k] << " " << inputDoi[i].sw[k] << " ";
// //       std::cout << inputDoi[i].z[k] << " " << inputDoi[i].sz[k] << " ";
//     }
//     
//     std::cout << std::endl;
//   }
  double minSigma = 0.01;
  double maxSigma = 0.03;
  TH1F* sigmaWdoi = new TH1F("sigmaWdoi","Distribution of measured sigma w",50,minSigma,maxSigma);
  TH1F* sigmaWdoiCentral = new TH1F("sigmaWdoiCentral","Central Crystals - Distribution of measured sigma w",50,minSigma,maxSigma);
  double averageSigma = 0;
  double averageSigmaError = 0;
  double averageDoiResFWHM = 0;
  double averageDoiResFWHMerr = 0;
  double averageEnResFWHM = 0;
//   double averageLO = 0;
  
  
  //-----------------------------//
  //         Photopeaks
  //-----------------------------//
  
  TCanvas *c;/* = (TCanvas*) gDirectory->Get("Energy res FWHM vs. i,j");*/
  TCanvas *cc;
  TH2F *spectrum2d;
  TH1F* all;
  TH1F* central;
  THStack* hs;
  TLegend* legend;
  double sumEN_corr = 0;
  double stdEN_corr = 0;
  double sumLO = 0;
  double stdLO = 0;
  f->cd("Module 0.0");
  
  TCanvas *c2 = (TCanvas*) gDirectory->Get("Photopeak positions vs. i,j");
  spectrum2d = (TH2F*)c2->GetPrimitive("Photopeak positions vs. i,j");
  all = new TH1F("all","",100,0,12000);
  all->GetXaxis()->SetTitle("Photopeak Position [ADC Channels]");
  //all->GetXaxis()->SetLabelSize(0.5);
  all->SetName("511 KeV Photopeak Position");
  //   all->GetYaxis()->SetTitle("N");
  all->SetStats(0);
  central = new TH1F("central","",100,0,12000);
  //   all->GetYaxis()->SetTitle("N");
  central->SetStats(0);
  central->SetFillStyle(3001);
  central->SetFillColor(kRed);
  
  int validEntry = 0;
  for(int i = 0 ; i < nmodulex*nmppcx*ncrystalsx ; i++)
  {
    for(int j = 0 ; j < nmoduley*nmppcy*ncrystalsy ; j++)
    {
      if(i > (ncrystalsx - 1) && i < nmodulex*nmppcx*ncrystalsx - (ncrystalsx) && j > (ncrystalsy-1) && j < nmoduley*nmppcy*ncrystalsy - (ncrystalsy)) //only crystals not from frame channels
      {
	if(spectrum2d->GetBinContent(i+1,j+1))
        {
	  central->Fill(spectrum2d->GetBinContent(i+1,j+1));
	  sumLO += spectrum2d->GetBinContent(i+1,j+1);
	  validEntry++;
	}
      }
      all->Fill(spectrum2d->GetBinContent(i+1,j+1));
    } 
  }
  hs = new THStack("hs","");
  double averageLO = sumLO/validEntry ;
  
  for(int i = 0 ; i < nmodulex*nmppcx*ncrystalsx ; i++)
  {
    for(int j = 0 ; j < nmoduley*nmppcy*ncrystalsy ; j++)
    {
      if(i > (ncrystalsx - 1) && i < nmodulex*nmppcx*ncrystalsx - (ncrystalsx) && j > (ncrystalsy-1) && j < nmoduley*nmppcy*ncrystalsy - (ncrystalsy)) //only crystals not from frame channels
      {
	if(spectrum2d->GetBinContent(i+1,j+1))
        {
	  stdLO += (1.0/validEntry)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageLO,2));
	}
      }
    }
  }
  
  cc = new TCanvas("cc","",1200,800);
  all->SetFillStyle(3001);
  all->SetFillColor(kBlue);
  cc->cd();
  hs->Add(all);
  hs->Add(central);
  //   all->Draw();
  //   central->Draw("same");
  legend = new TLegend(0.15,0.62,0.54,0.89,"");
  legend->SetFillStyle(0);
  legend->AddEntry(all,"All channels","f");
  legend->AddEntry(central,"Central channels","f");
  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("ADC Channels");
  hs->GetXaxis()->SetTitleOffset(1);
  hs->GetXaxis()->SetTitleSize(0.045);
  hs->GetXaxis()->SetLabelSize(0.045);
  hs->GetYaxis()->SetLabelSize(0.045);
  //hs->GetXaxis()->SetTickSize(0.08);
  legend->Draw();
  
  spectrum2d->GetYaxis()->SetTitleOffset(1.8);
  spectrum2d->GetXaxis()->SetTitleOffset(1.5);
  spectrum2d->GetZaxis()->SetTitleOffset(2.1);
  spectrum2d->GetXaxis()->SetLabelSize(0.045);
  spectrum2d->GetYaxis()->SetLabelSize(0.045);
  spectrum2d->GetZaxis()->SetLabelSize(0.045);
  spectrum2d->GetYaxis()->SetTitleSize(0.045);
  spectrum2d->GetXaxis()->SetTitleSize(0.045);
  spectrum2d->GetZaxis()->SetTitleSize(0.045);
  spectrum2d->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
  spectrum2d->GetYaxis()->SetNdivisions(8,2,0, kTRUE);
  
  TCanvas *clightOutput2D = new TCanvas("clightOutput2D","clightOutput2D",800,800);
  clightOutput2D->SetLeftMargin(0.18);
  spectrum2d->Draw("LEGO2");   
  //     clightOutput2D->Print("./plots/lightOutput2D.png");
  //     cc->Print("./plots/lightOutput.png");
  fOut->cd();
  cc->Write();
  delete spectrum2d;
  delete all;
  delete central;
  delete hs;
  delete cc;
  delete c2;
  delete legend;
  std::cout << "Light Output Central Channels [ADC Ch.] = " << averageLO << " +/- " << stdLO << std::endl;
  
  
  //-----------------------------//
  //      En Res Corrected
  //-----------------------------//
  
  f->cd("Module 0.0");
  c = (TCanvas*) gDirectory->Get("Corrected Energy res FWHM vs. i,j");
  
  spectrum2d = new TH2F();
  spectrum2d = (TH2F*)c->GetPrimitive("Corrected Energy res FWHM vs. i,j");
  all = new TH1F("Corrected Energy res FWHM","",250,0,0.5);
  all->GetXaxis()->SetTitle("Corrected Energy Resolution FWHM");
  all->SetName("Corrected Energy Resolution FWHM");
  //   all->GetYaxis()->SetTitle("N");
  all->SetStats(0);
  central = new TH1F("central","",250,0,0.5);
  //   all->GetYaxis()->SetTitle("N");
  central->SetStats(0);
  central->SetFillStyle(3001);
  central->SetFillColor(kRed);
  
  validEntry = 0;
  
  for(int i = 0 ; i < nmodulex*nmppcx*ncrystalsx ; i++)
  {
    for(int j = 0 ; j < nmoduley*nmppcy*ncrystalsy ; j++)
    {
      if(i > (ncrystalsx - 1) && i < nmodulex*nmppcx*ncrystalsx - (ncrystalsx) && j > (ncrystalsy-1) && j < nmoduley*nmppcy*ncrystalsy - (ncrystalsy)) //only crystals not from frame channels
      {
	if(spectrum2d->GetBinContent(i+1,j+1))
        {
	  central->Fill(spectrum2d->GetBinContent(i+1,j+1));
	  sumEN_corr += spectrum2d->GetBinContent(i+1,j+1);
	  validEntry++;
	}
      }
      all->Fill(spectrum2d->GetBinContent(i+1,j+1));
    } 
  }
  hs = new THStack("hs","");
  double averageEN_corr = sumEN_corr/validEntry ;
  
  for(int i = 0 ; i < nmodulex*nmppcx*ncrystalsx ; i++)
  {
    for(int j = 0 ; j < nmoduley*nmppcy*ncrystalsy ; j++)
    {
      if(i > (ncrystalsx - 1) && i < nmodulex*nmppcx*ncrystalsx - (ncrystalsx) && j > (ncrystalsy-1) && j < nmoduley*nmppcy*ncrystalsy - (ncrystalsy)) //only crystals not from frame channels
      {
	if(spectrum2d->GetBinContent(i+1,j+1))
        {
	  stdEN_corr += (1.0/validEntry)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageEN_corr,2));
	}
      }
    }
  }
  cc = new TCanvas("cc","",1200,800);
  all->SetFillStyle(3001);
  all->SetFillColor(kBlue);
  cc->cd();
  hs->Add(all);
  hs->Add(central);
  //   all->Draw();
  //   central->Draw("same");
  legend = new TLegend(0.5,0.62,0.893,0.89,"");
  legend->SetFillStyle(0);
  legend->AddEntry(all,"All channels","f");
  legend->AddEntry(central,"Central channels","f");
  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("Energy Resolution FWHM");
  hs->GetXaxis()->SetTitleOffset(1);
  hs->GetXaxis()->SetTitleSize(0.045);
  hs->GetXaxis()->SetLabelSize(0.045);
  hs->GetYaxis()->SetLabelSize(0.045);
  legend->Draw();
  
  spectrum2d->GetYaxis()->SetTitleOffset(1.8);
  spectrum2d->GetXaxis()->SetTitleOffset(1.5);
  spectrum2d->GetZaxis()->SetTitleOffset(1.7);
  spectrum2d->GetXaxis()->SetLabelSize(0.045);
  spectrum2d->GetYaxis()->SetLabelSize(0.045);
  spectrum2d->GetZaxis()->SetLabelSize(0.045);
  spectrum2d->GetYaxis()->SetTitleSize(0.045);
  spectrum2d->GetXaxis()->SetTitleSize(0.045);
  spectrum2d->GetZaxis()->SetTitleSize(0.045);
  spectrum2d->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
  spectrum2d->GetYaxis()->SetNdivisions(8,2,0, kTRUE);
  
  TCanvas *cenergyResolution2DCorrected = new TCanvas("cenergyResolution2DCorrected","cenergyResolution2DCorrected",800,800);
  cenergyResolution2DCorrected->SetLeftMargin(0.16);
  spectrum2d->Draw("LEGO2");
  //     cenergyResolution2DCorrected->Print("./plots/energyResolution2DCorrected.png");
  //     cc->Print("./plots/energyResolutionCorrected.png");
  fOut->cd();
  cc->Write();
  delete spectrum2d;
  delete all;
  delete central;
  delete hs;
  delete cc;
  delete legend;
  
    std::cout << "En. Res Corrected Central Channels [FWHM] = " << averageEN_corr << " +/- " << stdEN_corr << std::endl;
  
  //-----------------------------//
  //         DOI
  //-----------------------------//
    
  //DOI precision. Do this only if a doi file was given
  if(doDoiPrecision)
  {
    //simple histogram of measured sigma w
    
    for(int i = 0; i < inputDoi.size();i++)
    {
      for(int k =0 ; k < pointsFromDoi ; k++)
      {
	//       std::cout << inputDoi[i].sw[k] * inputDoi[i].sqrt_nentries[k] << std::endl;
	//       if(i > 1 && i < 6 && j > 1 && j < 6)
	//       {
	sigmaWdoi->Fill(inputDoi[i].sw[k] * inputDoi[i].sqrt_nentries[k]);
	//       }
      }
    }
    fOut->cd();
    sigmaWdoi->Write();
    
    //simple histogram of measured sigma w - central channels 
    for(int i = 0 ; i < nmodulex*nmppcx*ncrystalsx ; i++)
    {
      for(int j = 0 ; j < nmoduley*nmppcy*ncrystalsy ; j++)
      {
	if(i > (ncrystalsx - 1) && i < nmodulex*nmppcx*ncrystalsx - (ncrystalsx) && j > (ncrystalsy-1) && j < nmoduley*nmppcy*ncrystalsy - (ncrystalsy)) //only crystals not from frame channels
	{
	  for(int k = 0; k < inputDoi.size(); k++)
	  {
	    int ik = inputDoi[k].i;
	    int jk = inputDoi[k].j;
	    if(i == ik && j == jk)
	    {
	      for(int ik =0 ; ik < pointsFromDoi ; ik++)
	      {
		sigmaWdoiCentral->Fill(inputDoi[k].sw[ik] * inputDoi[k].sqrt_nentries[ik]);
	      }
	    }
	  }
	}
      }
    }
    fOut->cd();
    TF1 *gaussF = new TF1("gaussF","[0]*exp(-0.5*((x-[1])/[2])**2)",minSigma,maxSigma);
    gaussF->SetParameter(1,sigmaWdoiCentral->GetMean());
    gaussF->SetParameter(2,sigmaWdoiCentral->GetRMS());
    sigmaWdoiCentral->Fit("gaussF","RQ");
//     averageSigma = gaussF->GetParameter(1);
//     averageSigmaError = gaussF->GetParError(1);
    
    averageSigma = sigmaWdoiCentral->GetMean();
    averageSigmaError = sigmaWdoiCentral->GetRMS()/TMath::Sqrt(sigmaWdoiCentral->GetEntries());
    
    
    sigmaWdoiCentral->Write();
    
    int calibGraphCounter = 0;
    double averageM = 0;
    double averageMErr = 0;
    
    // run on all the crystals, look for the calibration plot
    for(int iModule = 0; iModule < nmodulex ; iModule++)
    {
      for(int jModule = 0; jModule < nmoduley ; jModule++)
      {
	for(int iMppc = 1; iMppc < nmppcx -1 ; iMppc++) // don't run on frame mppcs
	{
	  for(int jMppc = 1; jMppc < nmppcy -1; jMppc++)// don't run on frame mppcs
	  {
	    for(int iCry = 0; iCry < ncrystalsx ; iCry++)
	    {
	      for(int jCry = 0; jCry < ncrystalsy ; jCry++)
	      {
		int cryI = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
		int cryJ = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
		int crystalNumber = (cryI*ncrystalsx*nmppcx + cryJ);
		std::stringstream MppcDirStream;
		MppcDirStream << "MPPC " << letter[jMppc] << number[iMppc] << " - 0.0-" << iMppc << "." << jMppc;
		std::stringstream CrystalDirStream;
		CrystalDirStream << "Crystal " <<  crystalNumber;
		std::stringstream directory;
		directory << "Module 0.0/" << MppcDirStream.str() << "/" << CrystalDirStream.str();		
// 		std::cout << directory.str() << std::endl;
		f->cd();
		bool dir_exists = f->cd(directory.str().c_str()); //try to open the dir in root file
		
// 		int iCrystal = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
// 		int jCrystal = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
		
		
// 		int CentralCalibGraphCounter = 0;
		
		TCanvas* C_graph; 
		TGraph *calibGraph;
		
		if(dir_exists)
		{
		  std::stringstream stream;
		  stream << "Calibration Plot - Crystal " << crystalNumber;
		  C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
		  
		  if(C_graph) 
		    calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());
		  
		  if(calibGraph) //if calibration plot exist
		  {
		    calibGraphCounter++;
// 		    std::cout << calibGraph->GetName() << std::endl;
		    TCanvas* C_new = new TCanvas(stream.str().c_str(),stream.str().c_str(),1200,800);
		    
		    double crystaLenght = calibGraph->Eval(0); 
		    double yTopLimit    = crystaLenght - crystaLenght*0.10;
		    double yBottomLimit = 0 + crystaLenght*0.10;
		    
		    
		    double minBound = -1;
		    double maxBound = -1;
		    
		    //sample the calib plot and find the xlimits
		    int samples = 1000;
		    for(int nSample = 0 ; nSample < samples; nSample++)
		    {
		      double y = calibGraph->Eval((1.0/samples)*nSample);
		      if((minBound == -1) && (y < yTopLimit))
			minBound = (1.0/samples)*nSample;
		      if((maxBound == -1) && (y < yBottomLimit))
			maxBound = (1.0/samples)*nSample;
		    }
		    
// 		    std::cout << minBound << " " << maxBound << std::endl;
		    TF1 *lineFitSigma = new TF1("lineFitSigma","[0]*x + [1]",minBound,maxBound);
		    lineFitSigma->SetLineColor(kRed);
		    
		    calibGraph->Draw("AL");
		    calibGraph->Fit(lineFitSigma,"RQ");
		    
		    averageM += TMath::Abs(lineFitSigma->GetParameter(0));
		    averageMErr += TMath::Power(TMath::Abs(lineFitSigma->GetParError(0)),2);
		    
		    fOut->cd();
		    C_new->Write();
		  }
		}
	      }
	    }
	  }
	}
      }
    }
    averageMErr = TMath::Sqrt(averageMErr);
    averageDoiResFWHM = 2.355 * (averageM * averageSigma) / calibGraphCounter;
    averageDoiResFWHMerr = averageDoiResFWHM * TMath::Sqrt( TMath::Power(averageMErr/averageM,2) + TMath::Power(averageSigmaError/averageSigma,2) );
    
    std::cout << "DOI Resolution FWHM Central Channels [mm] = " << averageDoiResFWHM << " +/- " << averageDoiResFWHMerr << std::endl;
    
  }
  
  // doi bench specific part
  std::ofstream resultsFile;
  std::string resultsFileName = "results_" + rootFileNameNoExtension + ".txt";
  resultsFile.open(resultsFileName.c_str(), std::ofstream::out);
  
  resultsFile << "LO\t\tsLO\t\tEnRes\t\tsEnRes\t\tDoiRes\t\tsDoiRes" << std::endl;
  resultsFile << "[ADC Ch.]\t[ADC Ch.]\t[%]\t\t[%]\t\t[mm]\t\t[mm]" << std::endl;
  resultsFile << averageLO << "\t" <<"\t" << stdLO << "\t" << "\t"<<  averageEN_corr << "\t"<< stdEN_corr << "\t" << averageDoiResFWHM << "\t"<< "\t" << averageDoiResFWHMerr << std::endl;
  
//   std::cout << "Light Output Central Channels [ADC Ch.] = " << averageLO << " +/- " << stdLO << std::endl;
//   std::cout << "En. Res Corrected Central Channels FWHM [%] = " << averageEN_corr << " +/- " << stdEN_corr << std::endl;
//   std::cout << "DOI Resolution FWHM Central Channels [mm] = " << averageDoiResFWHM << " +/- " << averageDoiResFWHMerr << std::endl;
  
  
  fOut->Close();
  return 0;

}