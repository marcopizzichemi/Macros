// g++ -o ../build/modulePlot modulePlot.cpp `root-config --cflags --glibs` 
// -Wl,--no-as-needed -lHist -lCore -lMathCore

#include <stdio.h>
#include <stdlib.h>
#include <math.h>
#include <sys/types.h>
#include <getopt.h>
// #include <boost/lexical_cast.hpp>
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
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <TGraph.h>
#include <THStack.h>
#include <TLegend.h>
#include <TStyle.h>
#include <TMath.h>
#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>

const int pointsFromDoi = 4;

struct inputDoi_t
{
  int i;
  int j;
  double m;
  double q;
  double doires;
  double avgs;
  double w[pointsFromDoi];
  double sw[pointsFromDoi];
  double sqrt_nentries[pointsFromDoi];
  double z[pointsFromDoi];
  double sz[pointsFromDoi];
  
};

struct inputZ_t
{
  int i;
  int j;
  double z[pointsFromDoi];
};


double z[pointsFromDoi] = {10.8,8,5.2,2.4};
double sz[pointsFromDoi] = {0.5,0.5,0.5,0.5};

// double z[pointsFromDoi] = {13.6,10.8,8,5.2,2.4};
// double sz[pointsFromDoi] = {0.5,0.5,0.5,0.5,0.5};

int main(int argc, char * argv[])
{
  if(argc < 4)
  {
    std::cout << "USAGE:\t\t modulePlot moduleCalibration.root doiTagPoints.txt zPositions.txt [isBackgroudRun]" << std::endl;
    std::cout << std::endl;
    return 1;
  }
  
  //file with spectra from module calibration
  TFile *f = new TFile(argv[1]);
  f->cd("Module 0.0");
  
  TFile *fOut = new TFile("output.root","RECREATE");
  //file with doi tag m, q and points per crystal
  std::ifstream fDoiTag;
  fDoiTag.open(argv[2],std::ios::in);
  
  std::ifstream fZPos;
  fZPos.open(argv[3],std::ios::in);
  
  
//   std::ifstream fDoiResFile;
//   fDoiResFile.open(argv[3],std::ios::in);
  
  //see if this is a backgroudRun or not. If it is, no en res and light output calculations
  bool isBackgroudRun = 0;
  if(argc > 4)
    isBackgroudRun = atoi(argv[4]);
  
  std::vector<inputDoi_t> inputDoi;
  std::vector<inputZ_t> inputZ;
  while(!fZPos.eof())
  {
    int a,b;
    inputZ_t tempInput;
    
    fZPos >> tempInput.i >> tempInput.j;
    for(int k = 0 ; k < pointsFromDoi ; k++)
    {
      fZPos >> tempInput.z[k];
    }
    
    if(!fZPos.eof())
    {
      inputZ.push_back(tempInput);
    }
  }
  
  //DEBUG  
//   for(int i = 0 ; i < inputZ.size() ; i++)
//   {
//     std::cout << inputZ[i].i << " " << inputZ[i].j;
//     for(int k = 0 ; k < pointsFromDoi ; k++)
//     {
//       std::cout << " " << inputZ[i].z[k];
//     }
//     
//     std::cout << std::endl;
//   }
// 
//   std::cout << std::endl;
//   std::cout << std::endl;

//   double TagOffset = 2; // offset of tag derived by analisys of tagging bench. NOT USED ANYMORE
  //until we properly characterize the alignment of the tagging setup and modify the z positions at source, we use this temporary evaluation of the offset between the aligment of the x-y-z stage scale that we thought was correct and the one that we derive from the fine analisys of the tagging setup
  
  while(!fDoiTag.eof())
  {
    int a,b;
    inputDoi_t tempInput;
    
    
    fDoiTag >> tempInput.i >> tempInput.j >>tempInput.m >> tempInput.q >> tempInput.doires >> tempInput.avgs;
    
    tempInput.q = tempInput.q /*+ TagOffset*/; //until we properly characterize the alignment of the tagging setup and modify the z positions at source, we use this temporary evaluation of the offset between the aligment of the x-y-z stage scale that we thought was correct and the one that we derive from the fine analisys of the tagging setup
    
    for(int k = 0 ; k < pointsFromDoi ; k++)
    {
      fDoiTag >> tempInput.w[k] >> tempInput.sw[k] >> tempInput.sqrt_nentries[k];  
      //tempInput.z[k] = z[k] + TagOffset; //until we properly characterize the alignment of the tagging setup and modify the z positions at source, we use this temporary evaluation of the offset between the aligment of the x-y-z stage scale that we thought was correct and the one that we derive from the fine analisys of the tagging setup
      //tempInput.sz[k] = sz[k];
    }
    
    
    for(int k = 0; k < inputZ.size(); k++)
    {
      if(tempInput.i == inputZ[k].i && tempInput.j == inputZ[k].j)
      {
	for(int h = 0 ; h < pointsFromDoi ; h++)
        {
	  tempInput.z[h]  = inputZ[k].z[h];
	  tempInput.sz[h] = sz[h];
	}
      }
    }
    
    
    if(!fDoiTag.eof())
    {
      inputDoi.push_back(tempInput);
    }
  }
  fDoiTag.close();
  
// DEBUG
//   for(int i = 0 ; i < inputDoi.size() ; i++)
//   {
//     std::cout << inputDoi[i].i << " " << inputDoi[i].j << " " << inputDoi[i].m << " " << inputDoi[i].q << " ";
//     for(int k = 0 ; k < pointsFromDoi ; k++)
//     {
//       std::cout << inputDoi[i].w[k] << " " << inputDoi[i].sw[k] << " ";
//       std::cout << inputDoi[i].z[k] << " " << inputDoi[i].sz[k] << " ";
//     }
//     
//     std::cout << std::endl;
//   }

  
  
  
  
  struct stat st = {0};

  if (stat("./plots/", &st) == -1) {
    mkdir("./plots/", 0700);
  }
  
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	for(int k = 0; k < inputDoi.size(); k++)
// 	{
// 	  int ik = inputDoi[k].i;
// 	  int jk = inputDoi[k].j;
// 	  if(i == ik && j == jk)
// 	  {
  
  //simple histogram of measured sigma w
  TH1F* sigmaWdoi = new TH1F("sigmaWdoi","Distribution of measured sigma w",50,0.01,0.03);
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
  
  //simple histogram of measured sigma w - central 
  TH1F* sigmaWdoiCentral = new TH1F("sigmaWdoiCentral","Central Crystals - Distribution of measured sigma w",50,0.01,0.03);
  for(int i = 0 ; i < 8 ; i++)
  {
    for(int j = 0 ; j < 8 ; j++)
    {
      if(i > 1 && i < 6 && j > 1 && j < 6)
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
  sigmaWdoiCentral->Write();
  
  
  //-----------------------------//
  //      En Res Corrected
  //-----------------------------//
  
//   isBackgroudRun = true;//FIXME now hardcoded to true since this program is meant to analyze just the doi calibration
  if(!isBackgroudRun)
  {
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
    
    for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 8 ; j++)
      {
	if(i > 1 && i < 6 && j > 1 && j < 6)
	{
	  central->Fill(spectrum2d->GetBinContent(i+1,j+1));
	  sumEN_corr += spectrum2d->GetBinContent(i+1,j+1);
	}
	all->Fill(spectrum2d->GetBinContent(i+1,j+1));
      } 
    }
    hs = new THStack("hs","");
    double averageEN_corr = sumEN_corr/16 ;
    
    for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 8 ; j++)
      {
	if(i > 1 && i < 6 && j > 1 && j < 6)
	{
	  stdEN_corr += (1.0/16.0)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageEN_corr,2));
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
    cenergyResolution2DCorrected->Print("./plots/energyResolution2DCorrected.png");
    cc->Print("./plots/energyResolutionCorrected.png");
    fOut->cd();
    cc->Write();
    delete spectrum2d;
    delete all;
    delete central;
    delete hs;
    delete cc;
    delete legend;
    
    
    //-----------------------------//
    //         Photopeaks
    //-----------------------------//
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
    
    for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 8 ; j++)
      {
	if(i > 1 && i < 6 && j > 1 && j < 6)
	{
	  central->Fill(spectrum2d->GetBinContent(i+1,j+1));
	  sumLO += spectrum2d->GetBinContent(i+1,j+1);
	}
	all->Fill(spectrum2d->GetBinContent(i+1,j+1));
      } 
    }
    hs = new THStack("hs","");
    double averageLO = sumLO/16 ;
    
    for(int i = 0 ; i < 8 ; i++)
    {
      for(int j = 0 ; j < 8 ; j++)
      {
	if(i > 1 && i < 6 && j > 1 && j < 6)
	{
	  stdLO += (1.0/16.0)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageLO,2));
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
    clightOutput2D->Print("./plots/lightOutput2D.png");
    cc->Print("./plots/lightOutput.png");
    fOut->cd();
    cc->Write();
    delete spectrum2d;
    delete all;
    delete central;
    delete hs;
    delete cc;
    delete c2;
    delete legend;
    
    
    std::cout << "En. Res Corrected Central Channels [FWHM] = " << sumEN_corr/16 << " +/- " << stdEN_corr << std::endl;
    std::cout << "Light Output Central Channels [ADC Ch.] = " << sumLO/16 << " +/- " << stdLO << std::endl;
  }
  
  
  //DOI plots
  
  
  
  int nmodulex = 1;
  int nmoduley = 1;
  int nmppcx = 4;
  int nmppcy = 4;
  int ncrystalsx = 2;
  int ncrystalsy = 2;
  
//   f->cd("Module 0.0");
  
  std::string letter[4] = {"A","B","C","D"};
  std::string number[4] = {"1","2","3","4"};
  
  TH1F *histDoiPrecision = new TH1F("histDoiPrecision","histDoiPrecision",40,-5,5);
  TH1F *histLinePrecision = new TH1F("histLinePrecision","histLinePrecision",100,-5,5);
  
  std::ofstream outTagFile;
  outTagFile.open("test.dat",std::ios::out);
  
  
//   int points = 0;
  TH2F* diff2d = new TH2F("diff2d","diff2d",8,0,8,8,0,8);
  std::vector<double> xcorr,ycorr,xcorrCentral,ycorrCentral;
  
  for(int iModule = 0; iModule < nmodulex ; iModule++)
  {
    for(int jModule = 0; jModule < nmoduley ; jModule++)
    {
      
      for(int iMppc = 0; iMppc < nmppcx ; iMppc++)
      {
	for(int jMppc = 0; jMppc < nmppcy ; jMppc++)
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
// 	      std::cout << MppcDirStream.str() << std::endl;
	      
	      
	      std::stringstream CrystalDirStream;
	      CrystalDirStream << "Crystal " <<  crystalNumber;
// 	      std::cout << CrystalDirStream.str() << std::endl;
	      
	      std::stringstream directory;
	      directory << "Module 0.0/" << MppcDirStream.str() << "/" << CrystalDirStream.str();
	      
// 	      f->cd("Module 0.0");
	      std::cout << directory.str() << std::endl;
	      f->cd();
	      f->cd(directory.str().c_str());
	      
// 	      gDirectory->cd(MppcDirStream.str().c_str());
// 	      gDirectory->cd(CrystalDirStream.str().c_str());
	      
	      int iCrystal = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
	      int jCrystal = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
	      std::stringstream stream;
	      stream << "Calibration Plot - Crystal " << crystalNumber;
	      TCanvas* C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
	      TGraph *calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());
	      
	      
	      
	      for(int k = 0; k < inputDoi.size(); k++)
	      {
		int i = inputDoi[k].i;
		int j = inputDoi[k].j;
		if(iCrystal == i && jCrystal == j)
		{
		  TGraphErrors *points = new TGraphErrors(pointsFromDoi,inputDoi[k].w,inputDoi[k].z,inputDoi[k].sw,inputDoi[k].sz);
		  TF1 *lineTag = new TF1("lineTag","[0]*x + [1]",0,1);
		  lineTag->SetParameter(0,inputDoi[k].m);
		  lineTag->SetParameter(1,inputDoi[k].q);
		  points->Fit(lineTag,"Q");
		  
		  lineTag->SetLineColor(2);
		  double totalDelta = 0;
		  for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
		  {
		    double xPoint;
		    double yPoint;
		    points->GetPoint(pointNum,xPoint,yPoint);
		    histDoiPrecision->Fill(calibGraph->Eval(xPoint)-yPoint);
		    totalDelta += calibGraph->Eval(xPoint) - yPoint;
		    histLinePrecision->Fill(lineTag->Eval(xPoint)-yPoint);
		  }
		  diff2d->Fill(i,j,totalDelta/pointsFromDoi);
		  TCanvas* C_new = new TCanvas(stream.str().c_str(),stream.str().c_str(),1200,800);
		  
		  double crystaLenght = 15.0; //FIXME hardcoded!
		  double minBound = -1;
		  double maxBound = -1;
// 		  for(int step = 0 ; step < 10000 ; step ++)
// 		  {
// 		    if(minBound == -1)
// 		    {
// 		      if(calibGraph->Eval((crystaLenght/10000)*step) < 14)
// 		      {  
// 			minBound = (crystaLenght/10000)*step;
// 			
// 		      }
// 		        
// 		    }
// 		    if(maxBound == -1)
// 		    {
// 		      if(calibGraph->Eval((crystaLenght/10000)*step) < 1)
// 			 maxBound = (crystaLenght/10000)*step;
// 		    }
// 		  }
		  minBound = inputDoi[k].w[0];
		  maxBound = inputDoi[k].w[pointsFromDoi-1];
// 		  std::cout << minBound << " " << maxBound << std::endl;
		  TF1 *lineFitSigma = new TF1("lineFitSigma","[0]*x + [1]",minBound,maxBound);
		  lineFitSigma->SetLineColor(4);
		  
		  calibGraph->Draw("AL");
		  calibGraph->Fit(lineFitSigma,"RQ");
		  points->SetMarkerColor(4);
                  points->SetMarkerSize(1.2);
                  points->SetMarkerStyle(21);
		  points->Draw("P same");
		  fOut->cd();
                  C_new->Write();
		  lineTag->Draw("same");
		  lineFitSigma->Draw("same");
// 		  legend = new TLegend(0.5,0.62,0.893,0.89,"");
// 		  legend->SetFillStyle(0);
// 		  legend->AddEntry(points,"Measured point","f");
// 		  legend->AddEntry(lineTag,"Points Regression","f");
// 		  legend->AddEntry(calibGraph,"Cumulative Calibration","f");
// 		  legend->Draw("same");
		  std::stringstream fileName;
		  fileName << "./plots/CalibPlot" << crystalNumber << ".png"; 
		  outTagFile << i << "\t" << j << "\t" << lineTag->GetParameter(0) << "\t" << lineFitSigma->GetParameter(0) << std::endl;
		  xcorr.push_back(lineTag->GetParameter(0));
		  ycorr.push_back(lineFitSigma->GetParameter(0));
		  if(i > 1 && i < 6 && j > 1 && j < 6)
		  {
		    xcorrCentral.push_back(lineTag->GetParameter(0));
		    ycorrCentral.push_back(lineFitSigma->GetParameter(0));
		  }
		  
// 		  correlationPlot->SetPoint(points,lineTag->GetParameter(0),lineFitSigma->GetParameter(0));
// 		  points++;
// 		  outTagFile << i << "\t" << j << "\t" << inputDoi[k].m << "\t" << lineFitSigma->GetParameter(0) << std::endl;
		  C_new->Print(fileName.str().c_str());
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
  outTagFile.close();
  
  f->cd();
  f->cd("Module 0.0");
  TCanvas* C_Average = (TCanvas*) gDirectory->Get("Average DOI res FWHM vs. i,j");
  TH2F* AverageDOI = (TH2F*) C_Average->GetPrimitive("Average DOI res FWHM vs. i,j");
  
  TH1F* histoDoiTag = new TH1F("histoDoiTag","histoDoiTag",50,0,10);
  TH1F* histoDoiFromCalibration = new TH1F("histoDoiFromCalibration","histoDoiFromCalibration",50,0,10);
  TH1F* histoDoiTag_central = new TH1F("histoDoiTag_central","histoDoiTag_central",50,0,10);
  TH1F* histoDoiFromCalibration_central = new TH1F("histoDoiFromCalibration_central","histoDoiFromCalibration_central",50,0,10);
  
  TH2F *DoiResolutionVsIJ = new TH2F("DOI res FWHM vs. i,j","",nmppcx*ncrystalsx,0,nmppcx*ncrystalsx,nmppcy*ncrystalsy,0,nmppcy*ncrystalsy);
  DoiResolutionVsIJ->GetXaxis()->SetTitle("i (U axis)");
  DoiResolutionVsIJ->GetYaxis()->SetTitle("j (V axis)");
  DoiResolutionVsIJ->GetZaxis()->SetTitle("DOI Resolution FWHM [mm]");
  DoiResolutionVsIJ->GetXaxis()->SetTitleOffset(1.8);
  DoiResolutionVsIJ->GetYaxis()->SetTitleOffset(1.8);
  DoiResolutionVsIJ->GetZaxis()->SetTitleOffset(2.2);
  DoiResolutionVsIJ->GetZaxis()->SetRangeUser(0,10);
  
  std::vector<double> doiFromTag;
  std::vector<double> doiFromCalibration;
  std::vector<double> doiFromTag_central;
  std::vector<double> doiFromCalibration_central;
  
  
  for(int i = 0 ; i < 8 ; i++)
  {
    for(int j = 0 ; j < 8 ; j++)
    {
      if(i > 1 && i < 6 && j > 1 && j < 6)
      {
	for(int k = 0; k < inputDoi.size(); k++)
	{
	  int ik = inputDoi[k].i;
	  int jk = inputDoi[k].j;
	  if(i == ik && j == jk)
	  {
	    histoDoiTag_central->Fill(inputDoi[k].doires);
	    doiFromTag_central.push_back(inputDoi[k].doires);
	  }
	}
	histoDoiFromCalibration_central->Fill(AverageDOI->GetBinContent(i+1,j+1));
	doiFromCalibration_central.push_back(AverageDOI->GetBinContent(i+1,j+1));
      }
      
      for(int k = 0; k < inputDoi.size(); k++)
      {
	int ik = inputDoi[k].i;
	int jk = inputDoi[k].j;
	if(i == ik && j == jk)
	{
	  histoDoiTag->Fill(inputDoi[k].doires);
	  doiFromTag.push_back(inputDoi[k].doires);
	  DoiResolutionVsIJ->Fill(i,j,inputDoi[k].doires);
	}
      }
      histoDoiFromCalibration->Fill(AverageDOI->GetBinContent(i+1,j+1));
      doiFromCalibration.push_back(AverageDOI->GetBinContent(i+1,j+1));
      
    } 
  }
  
  fOut->cd();
  histoDoiTag->Write();
  histoDoiFromCalibration->Write();
  histoDoiTag_central->Write();
  histoDoiFromCalibration_central->Write();
  
  
  TCanvas* C_DoiResolutionVsIJ = new TCanvas("C_DoiResolutionVsIJ","C_DoiResolutionVsIJ",800,800);
  C_DoiResolutionVsIJ->SetName(DoiResolutionVsIJ->GetName());
  C_DoiResolutionVsIJ->cd();
  DoiResolutionVsIJ->Draw("LEGO2");
  C_DoiResolutionVsIJ->SetLeftMargin(0.15);
  C_DoiResolutionVsIJ->Write();
  
  TCanvas* C_correlationDoiRes = new TCanvas("C_correlationDoiRes","C_correlationDoiRes",1200,800);
  TGraph *correlationDoiRes = new TGraph(doiFromCalibration.size(),&doiFromTag[0],&doiFromCalibration[0]);
  correlationDoiRes->SetTitle("Correlation Plot DOI res");
  correlationDoiRes->GetXaxis()->SetTitle("DOI res FWHM from Tagging bench [mm]");
  correlationDoiRes->GetYaxis()->SetTitle("DOI res FWHM from Calibration [mm]");
  correlationDoiRes->Draw("A*");  
  fOut->cd();
  C_correlationDoiRes->Write();

  TCanvas* C_correlationDoiRes_central = new TCanvas("C_correlationDoiRes_central","C_correlationDoiRes_central",1200,800);  
  TGraph *correlationDoiRes_central = new TGraph(doiFromCalibration_central.size(),&doiFromTag_central[0],&doiFromCalibration_central[0]);
  correlationDoiRes_central->SetTitle("Correlation Plot DOI res - Central");
  correlationDoiRes_central->GetXaxis()->SetTitle("DOI res FWHM from Tagging bench [mm]");
  correlationDoiRes_central->GetYaxis()->SetTitle("DOI res FWHM from Calibration [mm]");
  correlationDoiRes_central->Draw("A*");
  fOut->cd();
  C_correlationDoiRes_central->Write();
  
  TCanvas* C_histDoiPrecision = new TCanvas("C_histDoiPrecision","C_histDoiPrecision",1200,800);
  C_histDoiPrecision->cd();
  histDoiPrecision->GetXaxis()->SetTitle("Delta (Calibration - Tagging bench) [mm]");
  histDoiPrecision->Draw();
//   gStyle->SetOptFit(1);
//   histDoiPrecision->Fit("gaus");
  C_histDoiPrecision->Print("./plots/HistoDOIprecision.png");
  fOut->cd();
  C_histDoiPrecision->Write();
  
  TCanvas* C_histLinePrecision = new TCanvas("C_histLinePrecision","C_histLinePrecision",1200,800);
  C_histLinePrecision->cd();
  histLinePrecision->GetXaxis()->SetTitle("Delta (Regression Tagging - Tagging bench) [mm]");
  histLinePrecision->Draw();
//   gStyle->SetOptFit(1);
//   histLinePrecision->Fit("gaus");
  C_histLinePrecision->Print("./plots/HistoLINEprecision.png");
  fOut->cd();
  C_histLinePrecision->Write();
  
  TGraph *correlationPlot = new TGraph(xcorr.size(),&xcorr[0],&ycorr[0]);
  TCanvas* C_correlationPlot = new TCanvas("C_correlationPlot","C_correlationPlot",1200,800);
  C_correlationPlot->cd();
  correlationPlot->SetTitle("Correlation Plot");
  correlationPlot->GetXaxis()->SetTitle("Tagging Bench m coeff.");
  correlationPlot->GetYaxis()->SetTitle("Cumulative plot m coeff.");
  correlationPlot->Draw("A*");
//   TF1* corr1line = new TF1("corr1line","x");
//   corr1line->Draw("same");
  std::cout << "Correlation = " << correlationPlot->GetCorrelationFactor() << std::endl;
  C_correlationPlot->Print("./plots/correlation.png");
  fOut->cd();
  C_correlationPlot->Write();
  
  TGraph *correlationPlotCentral = new TGraph(xcorrCentral.size(),&xcorrCentral[0],&ycorrCentral[0]);
  TCanvas* C_correlationPlotCentral = new TCanvas("C_correlationPlotCentral","C_correlationPlotCentral",1200,800);
  C_correlationPlotCentral->cd();
  correlationPlotCentral->SetTitle("Central Crystals - Correlation Plot");
  correlationPlotCentral->GetXaxis()->SetTitle("Tagging Bench m coeff.");
  correlationPlotCentral->GetYaxis()->SetTitle("Cumulative plot m coeff.");
  correlationPlotCentral->Draw("A*");
//   TF1* corr1lineCentral = new TF1("corr1lineCentral","x");
//   corr1lineCentral->Draw("same");
  std::cout << "Correlation = " << correlationPlotCentral->GetCorrelationFactor() << std::endl;
  C_correlationPlotCentral->Print("./plots/correlationCentral.png");
  fOut->cd();
  C_correlationPlotCentral->Write();
  
  
  TH1F* doiResFixedSigma = new TH1F("doiResFixedSigma","doiResFixedSigma",50,0,10);
  TH1F* mHisto = new TH1F("mHisto","mHisto",30,30,120);
  for(int i = 0 ; i < ycorr.size() ; i++)
  {
    doiResFixedSigma->Fill(std::abs(ycorr[i])*sigmaWdoi->GetMean()*2.355);
    mHisto->Fill(std::abs(ycorr[i]));
  }
  fOut->cd();
  mHisto->Write();
  doiResFixedSigma->Write();
  
  
  TH1F* doiResFixedSigma_central = new TH1F("doiResFixedSigma_central","doiResFixedSigma_central",50,0,10);
  TH1F* mHisto_central = new TH1F("mHisto_central","mHisto_central",30,30,120);
  for(int i = 0 ; i < ycorrCentral.size() ; i++)
  {
    doiResFixedSigma_central->Fill(std::abs(ycorrCentral[i])*sigmaWdoiCentral->GetMean()*2.355);
    mHisto_central->Fill(std::abs(ycorrCentral[i]));
  }
  fOut->cd();
  mHisto_central->Write();
  doiResFixedSigma_central->Write();
  
  
  std::cout << "Average m(ALL) = " << mHisto->GetMean() << " +/- " << mHisto->GetMeanError() << std::endl; 
  std::cout << "Average m(CENTRAL) = " << mHisto_central->GetMean() << " +/- " << mHisto_central->GetMeanError()  << std::endl; 
  
  std::cout << "Average sigmaW(ALL) = " << sigmaWdoi->GetMean() << " +/- " << sigmaWdoi->GetMeanError() <<  std::endl; 
  std::cout << "Average sigmaW(CENTRAL) = " << sigmaWdoiCentral->GetMean() << " +/- " << sigmaWdoiCentral->GetMeanError() << std::endl; 
  
  std::cout << "Average DOI Res FWHM [mm] - from DOI TAG, ALL = " << histoDoiTag->GetMean()<< " +/- " << histoDoiTag->GetMeanError()  << std::endl; 
  std::cout << "Average DOI Res FWHM [mm] - from DOI TAG, CENTRAL = "<< histoDoiTag_central->GetMean()<< " +/- " << histoDoiTag_central->GetMeanError()  << std::endl; 
  
  std::cout << "Average DOI Res FWHM [mm] - from CALIBRATION, ALL = " << mHisto->GetMean() * sigmaWdoi->GetMean() * 2.355 << " +/- " <<  mHisto->GetMean() * sigmaWdoi->GetMean() * 2.355 *TMath::Sqrt(TMath::Power((mHisto->GetMeanError()/mHisto->GetMean()),2) + TMath::Power((sigmaWdoi->GetMeanError()/sigmaWdoi->GetMean()),2) ) << std::endl; 
  std::cout << "Average DOI Res FWHM [mm] - from CALIBRATION, CENTRAL = "<< mHisto_central->GetMean() * sigmaWdoiCentral->GetMean() * 2.355<< " +/- " <<  mHisto_central->GetMean() * sigmaWdoiCentral->GetMean() * 2.355 *TMath::Sqrt(TMath::Power((mHisto_central->GetMeanError()/mHisto_central->GetMean()),2) + TMath::Power((sigmaWdoiCentral->GetMeanError()/sigmaWdoiCentral->GetMean()),2) )   << std::endl; 
  
  
  
  TCanvas* C_diff2d = new TCanvas("C_diff2d","C_diff2d",800,800);
  C_diff2d->cd();
  diff2d->Draw("COLZ");
  fOut->cd();
  C_diff2d->Write();
  
  fOut->Close();
  return 0;

}