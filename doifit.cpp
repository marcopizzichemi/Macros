//simple analysis

// compile with 
// g++ -o doifit ./Macros/doifit.cpp `root-config --cflags --glibs`


#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TH2F.h"
#include "TH1F.h"
#include "TChain.h"
#include "TF1.h"
#include "TCanvas.h"
#include "TGraphErrors.h"
#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include "TMath.h"

struct Crystal_t
{
  std::vector<double> x;
  std::vector<double> y;
  std::vector<double> w;
  std::vector<double> s;
  TCanvas *Cplot;
  TGraphErrors *doiR;
  double m,avgs,doires;
  
  
} Crystal[8];


int main (int argc, char** argv)
{
  std::vector<std::string> file;
//   file.push_back("./z0/doiData.txt");
//   file.push_back("./z1/doiData.txt");
  file.push_back("./z2/doiData.txt");
  file.push_back("./z3/doiData.txt");
  file.push_back("./z4/doiData.txt");
  
  
  std::vector<double> z;
//   z.push_back(18.87);
//   z.push_back(16.07);
  z.push_back(13.23);
  z.push_back(10.47);
  z.push_back(7.67);
  
  std::ifstream data[5];
  double x,y,w,s;
  for(int i = 0 ; i < file.size() ; i++)
  {
    data[i].open(file[i].c_str());
    int counter = 0;
    while(!data[i].eof())
    {
      data[i] >> x >> y >> w >> s;
      Crystal[counter].x.push_back(x);
      Crystal[counter].y.push_back(y);
      Crystal[counter].w.push_back(w);
      Crystal[counter].s.push_back(s);
      counter++;
    }
  }
  
  
//   for(int i = 0 ; i < file.size() ; i++)
//   {
//     int counter = 0;
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       std::cout <<  Crystal[counter].x[i] << " "<<  Crystal[counter].y[i] << " " << Crystal[counter].w[i] << " " << Crystal[counter].s[i] << std::endl;
//       counter++;
//     }
//     std::cout << std::endl;
//   }
  
  for(int i = 0 ; i < 8 ; i++)
  {
    double partial = 0;
    for(int j = 0 ; j < Crystal[i].s.size(); j++)
    {
      partial += Crystal[i].s[j];
    }
    Crystal[i].avgs = partial/Crystal[i].s.size();
  }
  
  
  std::stringstream tempStream;
  
  for(int j = 0 ; j < 8 ; j++)
  {
    std::cout << "Correlation Plot" << std::endl;
    tempStream.str("");
    tempStream << "Canvas_DOIvsR_" << j;
    
    Crystal[j].Cplot = new TCanvas(tempStream.str().c_str(),"",800,600);
    tempStream << ".png";
    std::string plotFile = tempStream.str().c_str();
    
    std::vector<double> tempEy,tempEx;
    for(int tt = 0 ; tt < 5  ;tt++)
    {
      tempEy.push_back(0.5);
      tempEx.push_back(0);
    }
    
    //Crystal[j].doiR = new TGraphErrors(5,&Crystal[j].x[0],&Crystal[j].y[0],&Crystal[j].ex[0],&Crystal[j].ey[0]);
    Crystal[j].doiR = new TGraphErrors(file.size(),&Crystal[j].w[0],&z[0],&tempEx[0],&tempEy[0]);
    Crystal[j].doiR->GetXaxis()->SetTitle("W");
    Crystal[j].doiR->GetYaxis()->SetTitle("DOI [mm]");
    tempStream.str("");
    tempStream << "DOI vs. R - Crystal " << j ;
    Crystal[j].doiR->SetTitle("");
    Crystal[j].Cplot->cd();
    Crystal[j].doiR->Draw("A*");
    
    TF1 *fit = new TF1("fit", "[0]+x*[1]",0,1);
    fit->SetParameter(0,80);
    fit->SetParameter(1,-100);
    Crystal[j].doiR->Fit("fit","QR");
    fit->SetLineColor(2);
    Crystal[j].Cplot->cd();
    fit->Draw("same");
    
    Crystal[j].Cplot->Print(plotFile.c_str());
//     Crystal[j].fitPar0 = fit->GetParameter(0);
//     Crystal[j].fitPar1 = fit->GetParameter(1);
    Crystal[j].m = fit->GetParameter(1);
    Crystal[j].doires = fit->GetParameter(1) * Crystal[j].avgs;
    std::cout << "Crystal " << j << " = " << fit->GetParameter(1) << " , " << Crystal[j].doires<< std::endl;
    
    
    
  }
  
  std::ofstream doi;
  doi.open("doi.txt");
  
  for(int j = 7 ; j > -1 ; j--)
  {
    doi << Crystal[j].x[0] << " " << Crystal[j].y[0] << " " << TMath::Abs(Crystal[j].doires) << std::endl;
  }
  doi.close();
  
  
  
  
  
  return 0;
}