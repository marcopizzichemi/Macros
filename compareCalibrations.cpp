// g++ -o ../build/compareCalibrations compareCalibrations.cpp `root-config --cflags --glibs` 
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
#include <TMultiGraph.h>

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


struct data_t
{
  std::string file;
  std::string title;
  int color;
  data_t (std::string a, std::string b, int c):file(a),title(b),color(c){}
};

int main(int argc, char * argv[])
{
  std::vector<data_t> dataset;
  data_t a("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/background_good/Run_2016_13_04_11_37_39/nearChannels.root","background",2);
  data_t b("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/lateral/Run_2016_05_04_08_48_47/nearChannels.root","lateral",3);
  data_t c("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/top/Run_2016_05_04_11_59_20/nearChannels.root","top",4);
  
  dataset.push_back(a);
  dataset.push_back(b);
  dataset.push_back(c);
  
  std::vector<TFile*> file;
  for(int i = 0 ; i < dataset.size(); i++)
  {
    TFile *temp = new TFile(dataset[i].file.c_str());
    file.push_back(temp);
  }
  
  //doi tag file
  std::vector<inputDoi_t> inputDoi;
  std::ifstream fDoiTag;
  fDoiTag.open("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/calibration_params.txt",std::ios::in);
  
  std::ifstream fZPos;
  fZPos.open("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/z_positions.txt",std::ios::in);
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
  fZPos.close();
  
  double TagOffset = 1.4;
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
  
  //output file
  std::stringstream fOutName;
  fOutName << "comparison";
  for(int i = 0 ; i < dataset.size(); i++)
  {
    fOutName << "-" << dataset[i].title;
  }
  fOutName << ".root";
  TFile *fOut = new TFile(fOutName.str().c_str(),"RECREATE");
  
  
  const int nmodulex = 1;
  const int nmoduley = 1;
  const int nmppcx = 4;
  const int nmppcy = 4;
  const int ncrystalsx = 2;
  const int ncrystalsy = 2;
  std::string letter[4] = {"A","B","C","D"};
  std::string number[4] = {"1","2","3","4"};
  
  
  //all 64 curves - prepare canvases
  int nCurves = nmodulex*nmoduley*nmppcx*nmppcy*ncrystalsx*ncrystalsy;
  std::vector<TCanvas *> C_curves;
  std::vector<TMultiGraph *> curves;
  for(int i = 0 ; i < nCurves ; i++)
  {
    std::stringstream title;
    title << "Curves - Crystal " << i;
    TCanvas *temp = new TCanvas(title.str().c_str(),title.str().c_str(),1200,800);
    C_curves.push_back(temp);
    TMultiGraph *mg = new TMultiGraph();
    curves.push_back(mg);
  }
  
  //doi precision
  std::vector<TH1F*> histDoiPrecision;
  std::vector<TH1F*> histDoiPrecision_central;
  for(int i = 0 ; i < dataset.size(); i++)
  {
    std::stringstream title;
    title << "Calibration Precision - " << dataset[i].title;
    TH1F* temp = new TH1F(title.str().c_str(),title.str().c_str(),60,-5,5);
    temp->SetLineWidth(3);
    temp->SetLineColor(dataset[i].color);
    histDoiPrecision.push_back(temp);
    
    
    title.str("");
    title << "Calibration Precision, Central Crystals - " << dataset[i].title;
    temp = new TH1F(title.str().c_str(),title.str().c_str(),60,-5,5);
    temp->SetLineWidth(3);
    temp->SetLineColor(dataset[i].color);
    histDoiPrecision_central.push_back(temp);
    
  }
  
  
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
	      //index and stuff
	      int cryI = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
	      int cryJ = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
	      int crystalNumber = (cryI*ncrystalsx*nmppcx + cryJ);
	      int iCrystal = (iModule*nmppcx*ncrystalsx)+(iMppc*ncrystalsx)+(iCry);
	      int jCrystal = (jModule*nmppcy*ncrystalsy)+(jMppc*ncrystalsy)+(jCry);
	      //compose dir string
	      std::stringstream MppcDirStream;
	      MppcDirStream << "MPPC " << letter[jMppc] << number[iMppc] << " - 0.0-" << iMppc << "." << jMppc;
	      std::stringstream CrystalDirStream;
	      CrystalDirStream << "Crystal " <<  crystalNumber;
	      std::stringstream directory;
	      directory << "Module 0.0/" << MppcDirStream.str() << "/" << CrystalDirStream.str();
	      
	      TLegend* legend;
	      legend = new TLegend(0.65,0.7,0.89,0.89,"");
	      legend->SetFillStyle(0);
	      
	      
	      TGraphErrors *points;
	      //draw the points
	      for(int k = 0; k < inputDoi.size(); k++)
	      {
		int i = inputDoi[k].i;
		int j = inputDoi[k].j;
		if(iCrystal == i && jCrystal == j)
		{
		  points = new TGraphErrors(pointsFromDoi,inputDoi[k].w,inputDoi[k].z,inputDoi[k].sw,inputDoi[k].sz);
		  points->SetLineColor(kBlack);
		  points->SetLineWidth(2);
		  points->SetMarkerStyle(21);
		  points->SetMarkerSize(1);
		  points->SetMarkerColor(kBlack);
		  curves[crystalNumber]->Add(points,"P");
		  legend->AddEntry(points,"tagging bench","lep");
		  // 		    TF1 *lineTag = new TF1("lineTag","[0]*x + [1]",0,1);
		  // 		    lineTag->SetParameter(0,inputDoi[k].m);
		  // 		    lineTag->SetParameter(1,inputDoi[k].q);
		  // 		    points->Fit(lineTag);
		  // 		    
		  // 		    lineTag->SetLineColor(2);
		  // 		    double totalDelta = 0;
		  // 		    for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
		  // 		    {
		  // 		      double xPoint;
		  // 		      double yPoint;
		  // 		      points->GetPoint(pointNum,xPoint,yPoint);
		  // 		      histDoiPrecision->Fill(calibGraph->Eval(xPoint)-yPoint);
		  // 		      totalDelta += calibGraph->Eval(xPoint) - yPoint;
		  // 		      histLinePrecision->Fill(lineTag->Eval(xPoint)-yPoint);
		  // 		    }
		  // 		    diff2d->Fill(i,j,totalDelta/pointsFromDoi);
		  // 		    TCanvas* C_new = new TCanvas(stream.str().c_str(),stream.str().c_str(),1200,800);
		  // 		    
		  // 		    double crystaLenght = 15.0;
		  // 		    double minBound = -1;
		  // 		    double maxBound = -1;
		  // 		    for(int step = 0 ; step < 10000 ; step ++)
		  // 		    {
		  // 		      if(minBound == -1)
		  // 		      {
		  // 			if(calibGraph->Eval((crystaLenght/10000)*step) < 14)
		  // 			{  
		  // 			  minBound = (crystaLenght/10000)*step;
		  // 			  
		  // 			}
		  // 			
		  // 		      }
		  // 		      if(maxBound == -1)
		  // 		      {
		  // 			if(calibGraph->Eval((crystaLenght/10000)*step) < 1)
		  // 			  maxBound = (crystaLenght/10000)*step;
		  // 		      }
		  // 		    }
		  // 		    std::cout << minBound << " " << maxBound << std::endl;
		  // 		    TF1 *lineFitSigma = new TF1("lineFitSigma","[0]*x + [1]",minBound,maxBound);
		  // 		    lineFitSigma->SetLineColor(4);
		  // 		    
		  // 		    calibGraph->Draw("AL");
		  // 		    calibGraph->Fit(lineFitSigma,"RQ");
		  // 		    points->Draw("P same");
		  // 		    lineTag->Draw("same");
		  // 		    lineFitSigma->Draw("same");
		  // 		    // 		  legend = new TLegend(0.5,0.62,0.893,0.89,"");
		  // 		    // 		  legend->SetFillStyle(0);
		  // 		    // 		  legend->AddEntry(points,"Measured point","f");
		  // 		    // 		  legend->AddEntry(lineTag,"Points Regression","f");
		  // 		    // 		  legend->AddEntry(calibGraph,"Cumulative Calibration","f");
		  // 		    // 		  legend->Draw("same");
		  // 		    std::stringstream fileName;
		  // 		    fileName << "./plots/CalibPlot" << crystalNumber << ".png"; 
		  // 		    outTagFile << i << "\t" << j << "\t" << lineTag->GetParameter(0) << "\t" << lineFitSigma->GetParameter(0) << std::endl;
		  // 		    xcorr.push_back(lineTag->GetParameter(0));
		  // 		    ycorr.push_back(lineFitSigma->GetParameter(0));
		  // 		    
		  // 		    // 		  correlationPlot->SetPoint(points,lineTag->GetParameter(0),lineFitSigma->GetParameter(0));
		  // 		    // 		  points++;
		  // 		    // 		  outTagFile << i << "\t" << j << "\t" << inputDoi[k].m << "\t" << lineFitSigma->GetParameter(0) << std::endl;
		  // 		    C_new->Print(fileName.str().c_str());
		  // 		    fOut->cd();
		  // 		    C_new->Write();
		}
	      }
	      
	      
	      //draw the curves
	      for(int i = 0 ; i < dataset.size(); i++)
	      {
		
		file[i]->cd();
		//go to the crystal
		file[i]->cd(directory.str().c_str());
		
		//take the graph
		std::stringstream stream;
		stream << "Calibration Plot - Crystal " << crystalNumber;
		TCanvas* C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
		TGraph *calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());
		
		
		for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
		{
		  double xPoint;
		  double yPoint;
		  points->GetPoint(pointNum,xPoint,yPoint);
		  histDoiPrecision[i]->Fill(calibGraph->Eval(xPoint)-yPoint);
		  if(iCrystal > 1 && iCrystal < 6 && jCrystal  > 1 && jCrystal < 6)
		  {
		    histDoiPrecision_central[i]->Fill(calibGraph->Eval(xPoint)-yPoint);
		  }
		  // 		    totalDelta += calibGraph->Eval(xPoint) - yPoint;
		  // 		    histLinePrecision->Fill(lineTag->Eval(xPoint)-yPoint);
		}
		
		
		//inclinations and offsets in case of "line"
		double crystaLenght = 15.0;
		double minBound = -1;
		double maxBound = -1;
		for(int step = 0 ; step < 10000 ; step ++)
		{
		  if(minBound == -1)
		  {
		    if(calibGraph->Eval((crystaLenght/10000)*step) < 14)
		    {  
		      minBound = (crystaLenght/10000)*step; 
		    }
		  }
		  if(maxBound == -1)
		  {
		    if(calibGraph->Eval((crystaLenght/10000)*step) < 1)
		      maxBound = (crystaLenght/10000)*step;
		  }
		}
		std::cout << minBound << " " << maxBound << std::endl;
		TF1 *lineFitSigma = new TF1("lineFitSigma","[0]*x + [1]",minBound,maxBound);
		// 		  lineFitSigma->SetLineColor(4);
		// 		  calibGraph->Draw("AL");
		calibGraph->Fit(lineFitSigma,"RNQ");
		
		
		//add the graph to the multigraph
		calibGraph->SetLineColor(dataset[i].color);
		calibGraph->SetLineWidth(2);
		curves[crystalNumber]->Add(calibGraph,"L");
		legend->AddEntry(calibGraph,dataset[i].title.c_str(),"l");
		
		
	      }
	      fOut->cd();
	      C_curves[crystalNumber]->cd();
	      // 		std::cout << directory.str() << std::endl;
	      // 		gPad->Update();
	      
	      
	      curves[crystalNumber]->Draw("a");
	      std::stringstream stream;
	      stream << "Calibrations comparison - Crystal " << crystalNumber;
	      
	      
	      curves[crystalNumber]->SetTitle(stream.str().c_str());
	      curves[crystalNumber]->GetXaxis()->SetTitle("W");
	      curves[crystalNumber]->GetYaxis()->SetTitle("Z [mm]");
	      legend->Draw();
	      
	      
	      C_curves[crystalNumber]->Write(); 
	      
	      
	      //comparisons part
	      for(int i = 0 ; i < dataset.size() - 1; i++)
	      {  
		for(int j = i+1 ; j < dataset.size(); j++)
		{
		  std::cout << dataset[i].title << " - " << dataset[j].title << std::endl;
		}
	      }
	      
	      
	      
	    }
	  }
	}
      }
    }
  }
  
  
  
  //   for(int i = 0 ; i < nCurves ; i++)
  //   {
  //     curves[i]->Write();
  //   }  
  //     
  
  

  fOut->cd();
  for(int i = 0 ; i < dataset.size() ; i++)
  {

    histDoiPrecision[i]->Write();  
    histDoiPrecision_central[i]->Write();  
  }
  

  TCanvas *C_allHisto = new TCanvas("Precision Comparison","Precision Comparison",1200,800);
  TLegend* legend;
  legend = new TLegend(0.65,0.7,0.89,0.89,"");
  
  legend->SetFillStyle(0);
  
  THStack *hs = new THStack("hs","");
  for(int i = 0 ; i < dataset.size() ; i++)
  {
    hs->Add(histDoiPrecision[i]);
    legend->AddEntry(histDoiPrecision[i],dataset[i].title.c_str(),"l");
  }
  hs->SetTitle("Precision Comparison");
  hs->Draw("nostack");
  hs->GetXaxis()->SetTitle("delta [mm]");
  legend->Draw();
  fOut->cd();
  C_allHisto->Write();
  
  
  TCanvas *C_allHisto_central = new TCanvas("Precision Comparison_central","Precision Comparison_central",1200,800);
  TLegend* legend_central;
  legend_central = new TLegend(0.65,0.7,0.89,0.89,"");
  legend_central->SetFillStyle(0);
  THStack *hs_central = new THStack("hs_central","");
  for(int i = 0 ; i < dataset.size() ; i++)
  {
    hs_central->Add(histDoiPrecision_central[i]);
    legend_central->AddEntry(histDoiPrecision_central[i],dataset[i].title.c_str(),"l");
  }
  
  hs_central->SetTitle("Precision Comparison_central");
  hs_central->Draw("nostack");
  hs_central->GetXaxis()->SetTitle("delta [mm]");
  legend_central->Draw();
  fOut->cd();
  C_allHisto_central->Write();
  
  
  
  fOut->Close();
  return 0;
}