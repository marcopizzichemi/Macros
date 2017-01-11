// compile with
// g++ -o ../build/convertToZ convertToZ.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/

#include "TROOT.h"
#include "TTree.h"
#include "TFile.h"
#include "TChain.h"
#include "TH2F.h"
#include "TObjArray.h"
#include "TObject.h"
#include "TH1F.h"
#include "TCanvas.h"
#include "TF1.h"
#include "TMath.h"
#include "TGraph.h"
#include "TGraph2DErrors.h"
#include "TGraph2D.h"
#include "TGraphDelaunay.h"
#include "TF2.h"
#include "TH3I.h"


#include "../code/struct.hh"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <algorithm>    // std::sort
#include <numeric>      // std::accumulate
#include <iomanip>      // std::setprecision
#include <vector>


struct CrystalData
{
  int id;
  double x;
  double y;
  int mppci;
  int mppcj;
  // TGraph* wz; //w(z) graph for this crystal
  TGraph* calibrationGraph; //z(w) graph for this crystal
  // TGraphDelaunay*** gd; //pointers to the pi_(w,E) for this crystal
  Float_t correction;
};


//MAIN
int main (int argc, char** argv)
{
  if(argc < 2) //check input, provide usage explainations
  {
    std::cout << "USAGE:\t\t\t convertToZ out0.root ... outN.root" << std::endl;
    std::cout << "a file calibration.root if assumed" << std::endl;
    std::cout << std::endl;
    return 1;
  }


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        ROOT STUFF                                    //
  //                                                                      //
  //----------------------------------------------------------------------//
  gROOT->ProcessLine("#include <vector>"); //needed by ROOT to deal with standard vectors
  //HACK to use the dictionary easily
  // Code taken from: http://www.gamedev.net/community/forums/topic.asp?topic_id=459511
  std::string fullFileName = "";
  std::string path = "";
  pid_t pid = getpid();
  char buf[20] = {0};
  sprintf(buf,"%d",pid);
  std::string _link = "/proc/";
  _link.append( buf );
  _link.append( "/exe");
  char proc[512];
  int ch = readlink(_link.c_str(),proc,512);
  if (ch != -1) {
    proc[ch] = 0;
    path = proc;
    std::string::size_type t = path.find_last_of("/");
    path = path.substr(0,t);
  }
  fullFileName = path + std::string("/");
  std::string command = ".L " + fullFileName + "structDictionary.C+";
  gROOT->ProcessLine(command.c_str());

  int numOfCry = 64;//this needs to be input somewhere...
  int numOfCh = 16;//this needs to be input somewhere...

  //input ttrees
  TChain *tree =  new TChain("tree"); // read input files
  for (int i = 1 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }

  //calibration file part
  CrystalData*** crystal;
  std::string mppcLabel [16] = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"};
  // std::string calibrationFileName = argv[1];
  std::string calibrationFileName = "calibration.root";
  TFile *calibrationFile = TFile::Open(calibrationFileName.c_str()); // hardcoded calibration file, output of ModuleCalibration on the out* files..
  int cryLateralNum = (int) sqrt(numOfCry);
  int cryLatPerMppc = ((int) sqrt(numOfCry)) / ((int) sqrt(numOfCh));
  double crystalPitch = 1.6;
  double matrixShiftX = (crystalPitch * cryLateralNum)/2.0 -  crystalPitch/2.0;
  double matrixShiftY = (crystalPitch * cryLateralNum)/2.0 -  crystalPitch/2.0;
  crystal = new CrystalData**[cryLateralNum];
  for(int i = 0 ; i < cryLateralNum ; i ++) crystal[i] = new CrystalData* [cryLateralNum];
  //fill them with positions etc
  for(int iCry = 0 ; iCry < cryLateralNum; iCry++)
  {
    for(int jCry = 0 ; jCry < cryLateralNum; jCry++)
    {

      crystal[iCry][jCry] = new CrystalData();
      crystal[iCry][jCry]->id = iCry*cryLateralNum + jCry;
      crystal[iCry][jCry]->x = iCry*crystalPitch - matrixShiftX;
      crystal[iCry][jCry]->y = jCry*crystalPitch - matrixShiftX;
      //find mppc i and j
      crystal[iCry][jCry]->mppci = iCry / cryLatPerMppc;
      crystal[iCry][jCry]->mppcj = jCry / cryLatPerMppc;
      if( (iCry > (cryLatPerMppc-1)) && (iCry < (cryLateralNum - cryLatPerMppc)) && (jCry > (cryLatPerMppc-1)) && (jCry < (cryLateralNum - cryLatPerMppc)) ) // avoid frame channels
      {
        //fetch the wz TGraph;
        // std::cout << iCry << "," << jCry << std::endl;
        std::stringstream  sscrystal;
        sscrystal << "Crystal " << crystal[iCry][jCry]->id;
        std::stringstream ssmppc;
        ssmppc << "MPPC " << mppcLabel[crystal[iCry][jCry]->mppcj] << crystal[iCry][jCry]->mppci+1 << " - 0.0-" << crystal[iCry][jCry]->mppci << "."<< crystal[iCry][jCry]->mppcj;
        std::stringstream sscryfolder;
        sscryfolder << "Module 0.0/" << ssmppc.str() << "/" << sscrystal.str();
        std::stringstream sswPlot;
        sswPlot << "Calibration Plot - " << sscrystal.str();
        // std::cout << sscryfolder.str().c_str() << std::endl;
        calibrationFile->cd(sscryfolder.str().c_str());
        TCanvas *canvas = (TCanvas*) gDirectory->Get(sswPlot.str().c_str());
        crystal[iCry][jCry]->calibrationGraph = (TGraph*) canvas->GetPrimitive(sswPlot.str().c_str());
        // TGraphDelaunay ***gd;
        // crystal[iCry][jCry]->gd = new TGraphDelaunay**[(int) sqrt(numOfCh)];
        // for(int igd = 0; igd < sqrt(numOfCh) ; igd++) crystal[iCry][jCry]->gd[igd] = new TGraphDelaunay* [(int) sqrt(numOfCh)];
        // TGraph2D ***camp;
        // camp = new TGraph2D**[(int) sqrt(numOfCh)];
        // for(int igd = 0; igd < sqrt(numOfCh) ; igd++) camp[igd] = new TGraph2D* [(int) sqrt(numOfCh)];

        //retrieve the calibrationGraph


        //calculate the TGraphDelaunays (after sampling)
        for(int iMPPC = 0; iMPPC < sqrt(numOfCh); iMPPC++)
        {
          for(int jMPPC = 0; jMPPC < sqrt(numOfCh) ; jMPPC++)
          {
            std::stringstream sstream;
            sstream << "Graph_Pi(E,w)[" << iMPPC <<  "][" << jMPPC <<  "]_" << crystal[iCry][jCry]->id;
            TCanvas *canv = (TCanvas*) gDirectory->Get(sstream.str().c_str());
            TGraph2D* graph2d = (TGraph2D*) canv-> GetPrimitive(sstream.str().c_str()); //get the 3d graph

            // int Npoints = graph2d->GetN();

            Double_t *x = graph2d->GetX();
            Double_t *y = graph2d->GetY();
            Double_t *z = graph2d->GetZ();
            std::ofstream myfile;
            std::stringstream fileStream;
            fileStream << "Gr_[" << iMPPC <<  "][" << jMPPC <<  "]_" << crystal[iCry][jCry]->id << ".dat";
            myfile.open (fileStream.str().c_str(),std::ios::out);
            for(int i = 0 ; i < graph2d->GetN() ; i++){
              myfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
            }
            myfile.close();
            sstream.str("");
             // std::cout <<  << std::endl;

            // sstream << "_camp";
            // camp[iMPPC][jMPPC] = new TGraph2D();
            // camp[iMPPC][jMPPC]->SetName(sstream.str().c_str());
            // int campPoint = 0;
            // for(int iBin = 1; iBin < histo3d->GetXaxis()->GetNbins()-1 ; iBin++)
            // {
            //   for(int jBin = 1; jBin < histo3d->GetYaxis()->GetNbins()-1 ; jBin++)
            //   {
            //     double sum = 0.0;
            //     double part = 0.0;
            //     for(int kBin = 1; kBin < histo3d->GetZaxis()->GetNbins()-1 ; kBin++)
            //     {
            //       if(histo3d->GetBinContent(iBin,jBin,kBin) > 0)
            //       {
            //         sum += histo3d->GetBinContent(iBin,jBin,kBin);
            //         part += histo3d->GetBinContent(iBin,jBin,kBin) * histo3d->GetZaxis()->GetBinCenter(kBin);
            //       }
            //     }
            //     if(sum > 0){
            //       part = part/sum;
            //       camp[iMPPC][jMPPC]->SetPoint(campPoint,histo3d->GetXaxis()->GetBinCenter(iBin),histo3d->GetYaxis()->GetBinCenter(jBin),part);
            //       campPoint++;
            //     }
            //   }
            // }
            // std::cout << crystal[iCry][jCry]->id << " " << iMPPC << " " << jMPPC << std::endl;
            // crystal[iCry][jCry]->gd[iMPPC][jMPPC] = new TGraphDelaunay(camp[iMPPC][jMPPC]);
            // crystal[iCry][jCry]->gd[iMPPC][jMPPC]->SetMaxIter(100000);
            // crystal[iCry][jCry]->gd[iMPPC][jMPPC]->SetMarginBinsContent(0);
            // crystal[iCry][jCry]->gd[iMPPC][jMPPC]->ComputeZ(0,0);
            // crystal[iCry][jCry]->gd[iMPPC][jMPPC]->FindAllTriangles();
            // sstream.str("");
          }
        }

        // std::stringstream sstream1;
        // sstream1 << "Charge Spectrum - Crystal " << crystal[iCry][jCry]->id << " - MPPC " <<  mppcLabel[crystal[iCry][jCry]->mppcj] << crystal[iCry][jCry]->mppci+1;
        // // std::cout << sstream1.str() << std::endl;
        // TCanvas *c_spectrum0 = (TCanvas*) gDirectory->Get(sstream1.str().c_str());
        // TH1F* spectrum0 = (TH1F*) c_spectrum0->GetPrimitive(sstream1.str().c_str());
        // sstream1.str("");
        // sstream1 << "gaussCharge - Crystal " << crystal[iCry][jCry]->id << " - MPPC " <<  mppcLabel[crystal[iCry][jCry]->mppcj] << crystal[iCry][jCry]->mppci+1;
        // TF1* gaussCharge0 = (TF1*) c_spectrum0->GetPrimitive(sstream1.str().c_str());
        // Float_t peak0 = gaussCharge0->GetParameter(1);
        // crystal[iCry][jCry]->correction = peak0 / 0.511;
      }
    }
  }

  return 0;
}
