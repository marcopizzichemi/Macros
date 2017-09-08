// g++ -o ../build/saturationAnalysis saturationAnalysis.cpp `root-config --cflags --glibs`
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
#include <TROOT.h>
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

// class of input from saturation file
struct inputSaturation_t
{
// public:
  std::string label;
  int mppcI;
  int mppcJ;
  int cryID;
  int cryI;
  int cryJ;
  float energy;
  float q;
  float sq;
  // inputSaturation_t(){};
  // friend std::istream& operator>>(std::istream& input, inputSaturation_t& s)
  // {
  //   input >> s.label;
  //   input >> s.mppcI;
  //   input >> s.mppcJ;
  //   input >> s.cryID;
  //   input >> s.cryI;
  //   input >> s.cryJ;
  //   input >> s.energy;
  //   input >> s.q;
  //   input >> s.sq;
  //
  //   return input;
  // }
};

struct channel_t
{
  std::string label;
  std::vector<float> energy;
  std::vector<float> q;
  std::vector<float> sq;
  std::vector<float> se;
};


int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cout << "Usage: saturationAnalysis inputFile [outputROOTfile] [outputTEXTfile]" << std::endl;
    std::cout << "       inputFile           = name of the input test file"<< std::endl;
    std::cout << "       outputROOTfile      = name of the output ROOT file - OPTIONAL, default to saturation_fits.root"<< std::endl;
    std::cout << "       outputTEXTfile      = name of the output TEXT file - OPTIONAL, dafault to saturation_parameters.txt"<< std::endl;
    return 1;
  }

  std::string inputFile = argv[1];
  std::string outputROOTfile = "saturation_fits.root";
  std::string outputTEXTfile = "saturation_parameters.txt";
  if(argc > 2)
  {
    outputROOTfile = argv[2];
  }
  if(argc > 3)
  {
    outputROOTfile = argv[3];
  }


  std::cout << "inputFile           = " << inputFile      << std::endl;
  std::cout << "outputROOTfile      = " << outputROOTfile << std::endl;
  std::cout << "outputTEXTfile      = " << outputTEXTfile << std::endl;
  TFile *fOut = new TFile(outputROOTfile.c_str(),"RECREATE");

  std::string label[16] = {"D1","C1","B1","A1","D2","C2","B2","A2","D3","C3","B3","A3","D4","C4","B4","A4"};
  channel_t mppc[16];
  for(int i = 0 ; i < 16 ; i++)
  {
    mppc[i].label = label[i];
  }

  std::ifstream fSaturation;
  fSaturation.open(inputFile.c_str(),std::ios::in);

  while(!fSaturation.eof())
  {
    std::string label;
    int mppcI;
    int mppcJ;
    int cryID;
    int cryI;
    int cryJ;
    float energy;
    float q;
    float sq;
    fSaturation >> label;
    fSaturation >> mppcI;
    fSaturation >> mppcJ;
    fSaturation >> cryID;
    fSaturation >> cryI;
    fSaturation >> cryJ;
    fSaturation >> energy;
    fSaturation >> q;
    fSaturation >> sq;

    if(!fSaturation.eof())
    {
      std::cout <<   label << "\t";
      std::cout <<   mppcI << "\t";
      std::cout <<   mppcJ << "\t";
      std::cout <<   cryID << "\t";
      std::cout <<   cryI << "\t";
      std::cout <<   cryJ << "\t";
      std::cout <<   energy << "\t";
      std::cout <<   q << "\t";
      std::cout <<   sq << std::endl;

      for(int i = 0 ; i < 16 ; i++)
      {
        if(mppc[i].label.compare(label) == 0)
        {
          mppc[i].energy.push_back(energy);
          mppc[i].q.push_back(q);
          mppc[i].sq.push_back(sq);
          mppc[i].se.push_back(0.0);
        }
      }
    }
  }

  std::vector<TGraphErrors*> graphs;
  std::vector<TCanvas*> C_graphs;

  for(int i = 0 ; i < 16 ; i++)
  {
    std::cout << mppc[i].label << std::endl;
    for(int j = 0 ; j < mppc[i].energy.size(); j++)
    {
      std::cout << mppc[i].energy[j] << "\t"
                << mppc[i].q[j] << "\t"
                << mppc[i].se[j] << "\t"
                << mppc[i].sq[j] << "\t"
                << std::endl;
    }
  }

  std::ofstream saturation_parameters;
  saturation_parameters.open(outputTEXTfile.c_str(),std::ios::out);

  for(int i = 0 ; i < 16 ; i++)
  {
    std::stringstream sname;
    sname << "Saturation " << mppc[i].label;
    // TCanvas* temp_Canvas = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    TGraphErrors* temp_graph = new TGraphErrors(mppc[i].energy.size(),&mppc[i].energy[0],&mppc[i].q[0],&mppc[i].se[0],&mppc[i].sq[0]);

    temp_graph->SetName(sname.str().c_str());
    //get max e and q
    float maxQ = 0.0;
    float maxE = 0.0;
    for(int j = 0 ; j < mppc[i].energy.size(); j++)
    {
      if(mppc[i].energy[j]>maxE)
      {
        maxE = mppc[i].energy[j];
      }
      if(mppc[i].q[j]>maxQ)
      {
        maxQ = mppc[i].q[j];
      }
    }

    TF1* fit = new TF1("fit","[0]*(1-exp(-x*[1]/[0]))",0,maxE);
    fit->SetParameter(0,maxQ);
    fit->SetParameter(1,maxQ/maxE); //to start as if there was no saturation...
    temp_graph->Fit("fit","R");

    saturation_parameters << mppc[i].label << "\t" << fit->GetParameter(0) << std::endl;
    graphs.push_back(temp_graph);
    // temp_Canvas->cd();
    // temp_graph->Draw("A*");
    // C_graphs.push_back(temp_Canvas);
  }

  fOut->cd();
  for(int i = 0 ; i < graphs.size() ; i++)
  {
    TCanvas * C_spectrum;
    C_spectrum = new TCanvas("C_spectrum","C_spectrum",1200,800);
    C_spectrum->SetName(graphs[i]->GetName());
    C_spectrum->cd();
    graphs[i]->Draw("A*");
    C_spectrum->Write();
    delete C_spectrum;
  }

  fOut->Close();
  saturation_parameters.close();
  fSaturation.close();



  return 0;
}
