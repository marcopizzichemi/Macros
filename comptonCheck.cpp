// compile with
// g++ -o ../build/comptonCheck comptonCheck.cpp `root-config --cflags --glibs` && cp structDictionary.C ../build/
// syntax
// comptonCheck `ls out*`

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



//Forward declarations
Float_t distance3D(Float_t ax, Float_t ay, Float_t az, Float_t bx, Float_t by, Float_t bz);
Float_t angle3D(Float_t a[3], Float_t b[3], Float_t c[3]);
Float_t scatteredGammaEnergy(Float_t energy,Float_t a[3],Float_t b[3],Float_t c[3]);
int generateZ(Float_t* zValue);
Float_t computeW(TGraph* wz,Float_t z);
// Float_t computeU(TGraphDelaunay ***gd,Double_t *xmppc,Float_t e,Float_t w,int mppcI0, int mppcJ0, int numOfCh);
Float_t computeU(TGraphDelaunay ***gd,Double_t *xmppc,Float_t charge,Float_t w,std::vector<int>& relevantMppcs,int numOfCh);
// Float_t computeV(TGraphDelaunay ***gd,Double_t *ymppc,Float_t e,Float_t w,int mppcI0, int mppcJ0, int numOfCh);
Float_t computeV(TGraphDelaunay ***gd,Double_t *ymppc,Float_t charge,Float_t w,std::vector<int>& relevantMppcs,int numOfCh);
Float_t averageW(Float_t e0,Float_t e1,TGraph* w0,Float_t z0,TGraph* w1,Float_t z1);
Float_t averageU(Float_t e0,Float_t e1,Float_t u0,Float_t u1);
Float_t averageV(Float_t e0,Float_t e1,Float_t v0,Float_t v1);
Float_t computeCrystalContribution(TF2* f2[4][4], Double_t position[16],Float_t charge,Float_t w);
Float_t simTravel1(Float_t pathLength, Float_t lambda511);
Float_t simCompton(Float_t comptonAngle);
Float_t simTravel2(Float_t pathLength, Float_t lambdaE);
Float_t simPhotoelectric(Float_t csPE);





//struct of averaged crystal events
struct avgCryEnergyDep
{
  int id;
  float x;
  float y;
  float z;
  float sx;
  float sy;
  float sz;
  float energy;
  float time;
};

//function to compare deposition event struct vectors using the field time
bool compareByTime(const enDep &a,const enDep  &b)
{
  return a.DepositionTime < b.DepositionTime;
}


//MAIN
int main (int argc, char** argv)
{
  if(argc < 2) //check input, provide usage explainations
  {
    std::cout << "USAGE:\t\t\t comptonCheck performAnalysis allAnalysis cEvent [inputFile1 inputFile2 ... inputFileN]" << std::endl;
    std::cout << "performAnalysis \t 1 = look for z0 and z1, 0 = don't" << std::endl;
    std::cout << "allAnalysis \t\t 1 = analyze all entries from input, 0 = analyze only 1 entry, specified by next value" << std::endl;
    std::cout << "cEvent \t\t\t entry to analyze if previous field is 1 (ignored if previous field is true)" << std::endl;
    std::cout << "inputFile1 \t\t list of .root files (output of simulation)" << std::endl;
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


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        CROSS SECTIONS ETC.                           //
  //                                                                      //
  //----------------------------------------------------------------------//
  Float_t vEnergyComptonPhotoelLYSO[103] = {
    5, 10, 15, 20, 25, 30, 35, 40, 45, 50, 55, 60, 65, 70, 75, 80, 85, 90, 95, 100,
    105, 110, 115, 120, 125, 130, 135, 140, 145, 150, 155, 160, 165, 170, 175, 180, 185, 190, 195, 200,
    205, 210, 215, 220, 225, 230, 235, 240, 245, 250, 255, 260, 265, 270, 275, 280, 285, 290, 295, 300,
    305, 310, 315, 320, 325, 330, 335, 340, 345, 350, 355, 360, 365, 370, 375, 380, 385, 390, 395, 400,
    405, 410, 415, 420, 425, 430, 435, 440, 445, 450, 455, 460, 465, 470, 475, 480, 485, 490, 495, 500,
    505, 510, 515

  };
  //change units to MeV in energy for photoelectric and compton cross section
  for(int i=0; i<103; i++)
  {
    vEnergyComptonPhotoelLYSO[i]=vEnergyComptonPhotoelLYSO[i]/1000.0;
  }

  Float_t vCSComptonLYSO[103] = {
    3.9897,6.5516,8.0894,9.1666,9.9455,10.505,10.902,11.18,11.37,11.494,11.568,11.606,11.615,11.603,
    11.574,11.534,11.484,11.427,11.365,11.299,11.23,11.16,11.088,11.015,10.942,10.868,10.796,10.723,
    10.651,10.58,10.509,10.44,10.371,10.303,10.236,10.17,10.105,10.041,9.9781,9.9161,9.855,9.7949,9.7357,
    9.6774,9.62,9.5635,9.5079,9.4532,9.3992,9.3461,9.2938,9.2423,9.1915,9.1415,9.0922,9.0436,8.9957,
    8.9485,8.9019,8.856,8.8107,8.7661,8.722,8.6786,8.6357,8.5934,8.5516,8.5104,8.4697,8.4296,8.3899,
    8.3507,8.312,8.2738,8.2361,8.1988,8.162,8.1255,8.0896,8.054,8.0189,7.9841,7.9498,7.9158,7.8822,7.849,
    7.8162,7.7837,7.7515,7.7197,7.6883,7.6572,7.6264,7.5959,7.5658,7.5359,7.5064,7.477,7.4482,7.4195,
    7.3912,7.363,7.3353
  };
  Float_t vCSPhotoelectricLYSO[103] = {
    417.63,169.13,93.706,43.83,23.927,14.502,9.4737,6.5458,4.724,3.5298,2.7133,2.1355,
    8.9356,7.3898,6.1779,5.216,4.4434,3.8162,3.302,2.8766,2.5217,2.2233,1.9706,1.7552,1.5704,1.411,
    1.2728,1.1524,1.0469,0.95419,0.87231,0.79973,0.73516,0.67752,0.6259,0.57953,0.53776,0.50002,
    0.46584,0.4348,0.40655,0.38077,0.35721,0.33562,0.31581,0.29759,0.28079,0.2653,0.25097,0.2377,
    0.2254,0.21397,0.20334,0.19344,0.1842,0.17558,0.16752,0.15997,0.15289,0.14625,0.14002,0.13415,
    0.12863,0.12343,0.11853,0.1139,0.10952,0.10538,0.10147,0.097761,0.094244,0.090909,0.087741,0.084732,
    0.08187,0.079147,0.076554,0.074083,0.071728,0.06948,0.067335,0.065285,0.063327,0.061453,0.059661,
    0.057945,0.056301,0.054725,0.053214,0.051765,0.050374,0.049037,0.047754,0.046519,0.045333,0.044191,
    0.043092,0.042034,0.041014,0.040025,0.039078,0.038164,0.037282
  };
  //change units to mm2/g in photoelectric cross section
  for(int i=0; i<103; i++)
  {
    vCSPhotoelectricLYSO[i]=vCSPhotoelectricLYSO[i]*100;
  }
  Float_t vEnergyLambdaLYSO[23] = {
    0.001,0.0015,0.002,0.003,0.004,0.005,0.006,0.008,0.01,0.015,0.02,
    0.03,0.04,0.05,0.06,0.08,0.1,0.15,0.2,0.3,0.4,0.5,0.6

  };
  Float_t vLambdaLYSO[23] = {
    4.0546e-05,9.7126e-05,4.6136e-05,9.6224e-05,0.00019359,0.00033754,
    0.00053432,0.0011114,0.00078548,0.0014043,0.0029794,0.008657,0.018449,
    0.032981,0.052543,0.024507,0.043318,0.11986,0.23536,0.53397,0.84433,
    1.1232,1.6516
  };
  //change units to mm in lambda
  for(int i=0; i<23; i++)
  {
    vLambdaLYSO[i] = vLambdaLYSO[i]*10;
  }
  //create canvas
  TCanvas* Canvas1 = new TCanvas("Canvas1", "Canvas1", 1200, 800);
  TCanvas* Canvas2 = new TCanvas("Canvas2", "Canvas2", 1200, 800);
  TCanvas* Canvas3 = new TCanvas("Canvas3", "Canvas3", 1200, 800);
  TCanvas* Canvas4 = new TCanvas("Canvas4", "Canvas4", 1200, 800);
  TGraph* comptonCrossSectionLYSO = new TGraph(103, vEnergyComptonPhotoelLYSO, vCSComptonLYSO);
  comptonCrossSectionLYSO->SetNameTitle ("Gamma energy vs Compton cross section in LYSO","Gamma energy vs Compton cross section in LYSO");
  TGraph* photoelectricCrossSectionLYSO = new TGraph(103, vEnergyComptonPhotoelLYSO, vCSPhotoelectricLYSO);
  photoelectricCrossSectionLYSO->SetNameTitle ("Gamma energy vs photoelectric cross section in LYSO","Gamma energy vs photoelectric cross section in LYSO");
  TGraph* lambdaLYSO = new TGraph(23, vEnergyLambdaLYSO, vLambdaLYSO);
  lambdaLYSO->SetNameTitle ("Gamma energy vs lambda in LYSO","Gamma energy vs lambda in LYSO");
  Canvas2->cd();
  comptonCrossSectionLYSO->Draw("A*");
  comptonCrossSectionLYSO->GetXaxis()->SetTitle("gamma energy [MeV]");
  comptonCrossSectionLYSO->GetYaxis()->SetTitle("cross section [mm2/g]");
  Canvas3->cd();
  photoelectricCrossSectionLYSO->Draw("A*");
  photoelectricCrossSectionLYSO->GetXaxis()->SetTitle("gamma energy [MeV]");
  photoelectricCrossSectionLYSO->GetYaxis()->SetTitle("cross section [mm2/g]");
  Canvas4->cd();
  lambdaLYSO->Draw("A*");
  lambdaLYSO->GetXaxis()->SetTitle("gamma energy [MeV]");
  lambdaLYSO->GetYaxis()->SetTitle("cross lambda [mm]");


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        INPUT FILES                                   //
  //                                                                      //
  //----------------------------------------------------------------------//
  TChain *tree =  new TChain("tree"); // read input files
  for (int i = 4 ; i < argc ; i++)
  {
    std::cout << "Adding file " << argv[i] << std::endl;
    tree->Add(argv[i]);
  }
  // find the number of channels directly from the tchain file
  // before creating the variables
  // first, get the list of leaves
  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++){
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  std::string det_prefix("detector");
  for(int i = 0 ; i < nLeaves ; i++){
    if (!leavesName[i].compare(0, det_prefix.size(), det_prefix)) numOfCh++;
  }
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        ANALYSIS PARAMETERS                           //
  //                                                                      //
  //----------------------------------------------------------------------//
  //load the TGraphs from the calibration run
  // choose the crystals
  std::stringstream ss1(argv[1]);
  bool performAnalysis;
  ss1 >> performAnalysis;
  std::stringstream ss2(argv[2]);
  bool allAnalysis;
  ss2 >> allAnalysis;
  std::stringstream ss3(argv[3]);
  int cEvent;
  ss3 >> cEvent;

  bool isNewDataset = true;
  int mppcI0 = 1;
  int mppcJ0 = 2;
  int mppcI1 = 1;
  int mppcJ1 = 1;
  // int cry0N = 28;
  // int cry1N = 29;
  // double cry0x = -0.8;
  // double cry0y = 0.8;
  // double cry1x = -0.8;
  // double cry1y = 2.4;
  int cry0N = 28;
  int cry1N = 26;

  double cry0centerX = -0.8;
  double cry0centerY = +0.8;
  double cry1centerX = -0.8;
  double cry1centerY = -2.4;
  double cry0x = cry0centerX;
  double cry0y = cry0centerY;
  double cry1x = cry1centerX;
  double cry1y = cry1centerY;


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        INPUT TTREES                                  //
  //                                                                      //
  //----------------------------------------------------------------------//
  //create the branches in the input ttree and connect to the variables
  // ----Global variables
  // these are 1 number per TTree entry - so 1 number per gamma shot
  Long64_t Seed;                      // seed of the simulation (read every time, but always the same)
  int Run;                            // run id (usually just 1)(read every time, but always the same)
  int Event;                          // event id
  float totalEnergyDeposited;         // total energy deposited in this event, in all the matrix
  Float_t sourcex;
  Float_t sourcey;
  Float_t sourcez;
  Float_t sourceMomentumX;
  Float_t sourceMomentumY;
  Float_t sourceMomentumZ;
  int NumOptPhotons;                  // number of optical photons generated in this event, in the entire matrix
  int NumCherenkovPhotons;            // number of Cherenkov photons generated in this event, in the entire matrix
  // energy deposition, each gamma 511 event has a std::vector of struct (type enDep) with all the data of each energy deposition
  std::vector<enDep> *energyDeposition = 0;
  // Total number of photons detected in this event
  // for each TTree entry, a simple number saying how many optical photons entered that
  // specific detector, passed the PDE check and where "detected" (i.e. saved)
  Short_t  *detector;
  detector = new Short_t [numOfCh];
  // optical photons. for each gamma 511 event, every optical photon detected is a struct of type optPhot. a std::vector<optPhot> is saved for each gamma 511
  std::vector<optPhot> *photons = 0;
  // ----Brach addresses
  tree->SetBranchAddress("Seed",&Seed);
  tree->SetBranchAddress("Run",&Run);
  tree->SetBranchAddress("Event",&Event);
  tree->SetBranchAddress("totalEnergyDeposited",&totalEnergyDeposited);
  if(isNewDataset){
    tree->SetBranchAddress("SourceX",&sourcex);
    tree->SetBranchAddress("SourceY",&sourcey);
    tree->SetBranchAddress("SourceZ",&sourcez);
    tree->SetBranchAddress("SourceMomentumX",&sourceMomentumX);
    tree->SetBranchAddress("SourceMomentumY",&sourceMomentumY);
    tree->SetBranchAddress("SourceMomentumZ",&sourceMomentumZ);
  }
  tree->SetBranchAddress("NumOptPhotons",&NumOptPhotons);
  tree->SetBranchAddress("NumCherenkovPhotons",&NumCherenkovPhotons);
  tree->SetBranchAddress("optical",&photons);
  tree->SetBranchAddress("energyDeposition",&energyDeposition);
  for (int i = 0 ; i < numOfCh ; i++)
  {
    std::stringstream snames;
    snames << "detector" << i;
    tree->SetBranchAddress(snames.str().c_str(),&detector[i]);
  }
  // MPPC geometry
  // Double_t xmppc[16]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  // Double_t ymppc[16]={-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8,-4.8,-1.6,1.6,4.8};
  // MPPC geometry for 64 channels and 64x4 crystals
  Double_t *xmppc,*ymppc;
  xmppc = new Double_t [numOfCh];
  ymppc = new Double_t [numOfCh];
  Double_t mppcPitch = 3.2;
  int detCounter;
  for(int i = 0 ; i < sqrt(numOfCh) ; i++)
  {
    for(int j = 0 ; j < sqrt(numOfCh); j++)
    {
      detCounter = i*sqrt(numOfCh) + j;
      xmppc[detCounter] = (i+0.5*(1.0-sqrt(numOfCh)))*mppcPitch;
      ymppc[detCounter] = (j+0.5*(1.0-sqrt(numOfCh)))*mppcPitch;
    }
  }


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        FETCH CALIBRATION PLOTS                       //
  //                                                                      //
  //----------------------------------------------------------------------//
  std::string mppcLabel [16] = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N","O","P"};
  std::stringstream  crystal0;
  crystal0 << "Crystal " << cry0N;
  std::stringstream crystal1;
  crystal1 << "Crystal " << cry1N;
  std::stringstream mppc0;
  mppc0 << "MPPC " << mppcLabel[mppcJ0] << mppcI0+1 << " - 0.0-" << mppcI0 << "."<< mppcJ0;
  std::stringstream mppc1;
  mppc1 << "MPPC " << mppcLabel[mppcJ1] << mppcI1+1 << " - 0.0-" << mppcI1 << "."<< mppcJ1;
  std::stringstream cry0folder;
  std::stringstream cry1folder;
  cry0folder << "Module 0.0/" << mppc0.str() << "/" << crystal0.str();
  cry1folder << "Module 0.0/" << mppc1.str() << "/" << crystal1.str();
  std::stringstream w0Plot;
  std::stringstream w1Plot;
  w0Plot << "w(z) Plot - " << crystal0.str();
  w1Plot << "w(z) Plot - " << crystal1.str();
  TFile *calibrationFile = TFile::Open("calibration.root"); // hardcoded calibration file, output of ModuleCalibration on the out* files..
  std::cout << cry0folder.str().c_str() << std::endl;
  calibrationFile->cd(cry0folder.str().c_str());
  TCanvas *canvas = (TCanvas*) gDirectory->Get(w0Plot.str().c_str());
  TGraph *wzgraph0,*wzgraph1;
  wzgraph0 = (TGraph*) canvas->GetPrimitive(w0Plot.str().c_str());
  calibrationFile->cd(cry1folder.str().c_str());// cry1
  TCanvas* canvas1 = (TCanvas*) gDirectory->Get(w1Plot.str().c_str());
  wzgraph1 = (TGraph*) canvas1->GetPrimitive(w1Plot.str().c_str());
  //now load the TGraph2Ds from the calibration run
  //and calculate the "surfaces"
  TGraphDelaunay ***gd0;
  gd0 = new TGraphDelaunay**[(int) sqrt(numOfCh)];
  for(int i = 0; i < sqrt(numOfCh) ; i++) gd0[i] = new TGraphDelaunay* [(int) sqrt(numOfCh)];
  TGraphDelaunay ***gd1;
  gd1 = new TGraphDelaunay**[(int) sqrt(numOfCh)];
  for(int i = 0; i < sqrt(numOfCh) ; i++) gd1[i] = new TGraphDelaunay* [(int) sqrt(numOfCh)];
  TGraph2D ***camp0;
  TGraph2D ***camp1;
  camp0 = new TGraph2D**[(int) sqrt(numOfCh)];
  for(int i = 0; i < sqrt(numOfCh) ; i++) camp0[i] = new TGraph2D* [(int) sqrt(numOfCh)];
  camp1 = new TGraph2D**[(int) sqrt(numOfCh)];
  for(int i = 0; i < sqrt(numOfCh) ; i++) camp1[i] = new TGraph2D* [(int) sqrt(numOfCh)];


  bool dumpPoint = false; //dump points to a text file, in case you want to compare with a matlab analysis for example
  if(dumpPoint)
  {
    calibrationFile->cd(cry0folder.str().c_str());// cry0
    for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
    {
      for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
      {
        std::stringstream sstream;
        sstream << "Graph_Pi(E,w)[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry0N;
        TCanvas *canv = (TCanvas*) gDirectory->Get(sstream.str().c_str());
        TGraph2D *graph = (TGraph2D*) canv->GetPrimitive(sstream.str().c_str());
        Double_t *x = graph->GetX();
        Double_t *y = graph->GetY();
        Double_t *z = graph->GetZ();
        std::ofstream myfile;
        std::stringstream fileStream;
        fileStream << "Gr_[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry0N << ".dat";
        myfile.open (fileStream.str().c_str(),std::ios::out);
        for(int i = 0 ; i < graph->GetN() ; i++){
          myfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
        }
        myfile.close();
        sstream.str("");
      }
    }
    calibrationFile->cd(cry1folder.str().c_str());// cry1
    for(int iMPPC = mppcI1-1; iMPPC < mppcI1+2; iMPPC++)
    {
      for(int jMPPC = mppcJ1-1; jMPPC < mppcJ1+2 ; jMPPC++)
      {
        std::stringstream sstream;
        sstream << "Graph_Pi(E,w)[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry1N;
        TCanvas *canv = (TCanvas*) gDirectory->Get(sstream.str().c_str());
        TGraph2D *graph = (TGraph2D*) canv->GetPrimitive(sstream.str().c_str());
        Double_t *x = graph->GetX();
        Double_t *y = graph->GetY();
        Double_t *z = graph->GetZ();
        std::ofstream myfile;
        std::stringstream fileStream;
        fileStream << "Gr_[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry1N << ".dat";
        myfile.open (fileStream.str().c_str(),std::ios::out);
        for(int i = 0 ; i < graph->GetN() ; i++){
          myfile << x[i] << " " << y[i] << " " << z[i] << std::endl;
        }
        myfile.close();
        sstream.str("");
      }
    }
    std::cout << "Points dumped." << std::endl;
  }


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        TGraph2D INTERPOLATION                        //
  //                                                                      //
  //----------------------------------------------------------------------//
  int interpolationWay = 2;  // 0 = TGraphDelaunay, 1 = planes, 2 = TGraphDelaunay on sampled data
  //take and interpolate all TGraph2D, regardless of how many will be used later
  if(interpolationWay == 2)
  {
    calibrationFile->cd(cry0folder.str().c_str());// cry0
    // for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
    // {
    //   for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
    //   {
    for(int iMPPC = 0; iMPPC < sqrt(numOfCh); iMPPC++)
    {
      for(int jMPPC = 0; jMPPC < sqrt(numOfCh) ; jMPPC++)
      {
        std::stringstream sstream;
        sstream << "Pi_[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry0N;
        TCanvas *canv = (TCanvas*) gDirectory->Get(sstream.str().c_str());
        TH3I* histo3d = (TH3I*) canv-> GetPrimitive(sstream.str().c_str()); //get the 3d histo
        sstream << "_camp";
        camp0[iMPPC][jMPPC] = new TGraph2D();
        camp0[iMPPC][jMPPC]->SetName(sstream.str().c_str());
        int campPoint = 0;
        for(int iBin = 1; iBin < histo3d->GetXaxis()->GetNbins()-1 ; iBin++)
        {
          for(int jBin = 1; jBin < histo3d->GetYaxis()->GetNbins()-1 ; jBin++)
          {
            double sum = 0.0;
            double part = 0.0;
            for(int kBin = 1; kBin < histo3d->GetZaxis()->GetNbins()-1 ; kBin++)
            {
              if(histo3d->GetBinContent(iBin,jBin,kBin) > 0)
              {
                sum += histo3d->GetBinContent(iBin,jBin,kBin);
                part += histo3d->GetBinContent(iBin,jBin,kBin) * histo3d->GetZaxis()->GetBinCenter(kBin);
              }
            }
            if(sum > 0){
              part = part/sum;
              camp0[iMPPC][jMPPC]->SetPoint(campPoint,histo3d->GetXaxis()->GetBinCenter(iBin),histo3d->GetYaxis()->GetBinCenter(jBin),part);
              campPoint++;
            }
          }
        }
        std::cout << cry0N << " " << iMPPC << " " << jMPPC << std::endl;
        gd0[iMPPC][jMPPC] = new TGraphDelaunay(camp0[iMPPC][jMPPC]);
        gd0[iMPPC][jMPPC]->SetMaxIter(100000);
        gd0[iMPPC][jMPPC]->SetMarginBinsContent(0);
        gd0[iMPPC][jMPPC]->ComputeZ(0,0);
        gd0[iMPPC][jMPPC]->FindAllTriangles();
        sstream.str("");
      }
    }
    calibrationFile->cd(cry1folder.str().c_str());// cry1
    // for(int iMPPC = mppcI1-1; iMPPC < mppcI1+2; iMPPC++)
    // {
    //   for(int jMPPC = mppcJ1-1; jMPPC < mppcJ1+2 ; jMPPC++)
    //   {
    for(int iMPPC = 0; iMPPC < sqrt(numOfCh); iMPPC++)
    {
      for(int jMPPC = 0; jMPPC < sqrt(numOfCh) ; jMPPC++)
      {
        std::stringstream sstream;
        sstream << "Pi_[" << iMPPC <<  "][" << jMPPC <<  "]_" << cry1N;
        TCanvas *canv = (TCanvas*) gDirectory->Get(sstream.str().c_str());
        TH3I* histo3d = (TH3I*) canv-> GetPrimitive(sstream.str().c_str()); //get the 3d histo
        sstream << "_camp";
        camp1[iMPPC][jMPPC] = new TGraph2D();
        camp1[iMPPC][jMPPC]->SetName(sstream.str().c_str());
        int campPoint = 0;
        for(int iBin = 1; iBin < histo3d->GetXaxis()->GetNbins()-1 ; iBin++)
        {
          for(int jBin = 1; jBin < histo3d->GetYaxis()->GetNbins()-1 ; jBin++)
          {
            double sum = 0.0;
            double part = 0.0;
            for(int kBin = 1; kBin < histo3d->GetZaxis()->GetNbins()-1 ; kBin++)
            {
              if(histo3d->GetBinContent(iBin,jBin,kBin) > 0)
              {
                sum += histo3d->GetBinContent(iBin,jBin,kBin);
                part += histo3d->GetBinContent(iBin,jBin,kBin) * histo3d->GetZaxis()->GetBinCenter(kBin);
              }
            }
            if(sum > 0){
              part = part/sum;
              camp1[iMPPC][jMPPC]->SetPoint(campPoint,histo3d->GetXaxis()->GetBinCenter(iBin),histo3d->GetYaxis()->GetBinCenter(jBin),part);
              campPoint++;
            }
          }
        }
        std::cout << cry1N << " " << iMPPC << " " << jMPPC << std::endl;
        gd1[iMPPC][jMPPC] = new TGraphDelaunay(camp1[iMPPC][jMPPC]);
        gd1[iMPPC][jMPPC]->SetMaxIter(100000);
        gd1[iMPPC][jMPPC]->SetMarginBinsContent(0);
        gd1[iMPPC][jMPPC]->ComputeZ(0,0);
        gd1[iMPPC][jMPPC]->FindAllTriangles();
        sstream.str("");
      }
    }
  }


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        CHARGE CALIBRATION                            //
  //                                                                      //
  //----------------------------------------------------------------------//
  //take the charge spectra for the two channels 28 and 29
  calibrationFile->cd(cry0folder.str().c_str());// cry0
  std::stringstream sstream1;
  sstream1 << "Charge Spectrum - Crystal " << cry0N << " - MPPC " <<  mppcLabel[mppcJ0] << mppcI0+1;
  // std::cout << sstream1.str() << std::endl;
  TCanvas *c_spectrum0 = (TCanvas*) gDirectory->Get(sstream1.str().c_str());
  TH1F* spectrum0 = (TH1F*) c_spectrum0->GetPrimitive(sstream1.str().c_str());
  sstream1.str("");
  sstream1 << "gaussCharge - Crystal " << cry0N << " - MPPC " <<  mppcLabel[mppcJ0] << mppcI0+1;
  TF1* gaussCharge0 = (TF1*) c_spectrum0->GetPrimitive(sstream1.str().c_str());
  Float_t peak0 = gaussCharge0->GetParameter(1);
  sstream1.str("");
  calibrationFile->cd(cry1folder.str().c_str());// cry1
  sstream1 << "Charge Spectrum - Crystal " << cry1N << " - MPPC " <<  mppcLabel[mppcJ1] << mppcI1+1;
  TCanvas *c_spectrum1 = (TCanvas*) gDirectory->Get(sstream1.str().c_str());
  TH1F* spectrum1 = (TH1F*) c_spectrum1->GetPrimitive(sstream1.str().c_str());
  sstream1.str("");
  sstream1 << "gaussCharge - Crystal " << cry1N << " - MPPC " <<  mppcLabel[mppcJ1] << mppcI1+1;
  TF1* gaussCharge1 = (TF1*) c_spectrum1->GetPrimitive(sstream1.str().c_str());
  Float_t peak1 = gaussCharge1->GetParameter(1);
  //assuming linear relation charge - energy deposited...
  Float_t correction0 = peak0 / 0.511; //charge bin / MeV
  Float_t correction1 = peak1 / 0.511; //charge bin / MeV


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        OUTPUT HISTOGRAMS                             //
  //                                                                      //
  //----------------------------------------------------------------------//
  double columsum = 0;
  double rowsum = 0;
  double total = 0;
  double frontSum = 0;
  double floodx = 0;
  double floody = 0;
  double floodz = 0;
  long int foundCandidate = 0;
  long int singleCounter = 0;
  long int doubleCounter = 0;
  long int tripleCounter = 0;
  long int multipleCounter = 0;
  long int globalReturned = 0;
  TH2F *flood = new TH2F("FloodHisto","FloodHisto",1000,-7,7,1000,-7,7);
  flood->GetXaxis()->SetTitle("u");
  flood->GetYaxis()->SetTitle("v");
  flood->GetZaxis()->SetTitle("N");

  TH1F *histoProb = new TH1F("histoProb","histoProb",1000,0,30);

  TH2F *RealEnergy_RealCoord = new TH2F("RealEnergy_RealCoord","RealEnergy_RealCoord",100,0,1,100,0,1);
  RealEnergy_RealCoord->GetXaxis()->SetTitle("RealCoord");
  RealEnergy_RealCoord->GetYaxis()->SetTitle("RealEnergy");

  TH2F *RealEnergy_DetectorCoord = new TH2F("RealEnergy_DetectorCoord","RealEnergy_DetectorCoord",100,0,1,100,0,1);
  RealEnergy_DetectorCoord->GetXaxis()->SetTitle("DetectorCoord");
  RealEnergy_DetectorCoord->GetYaxis()->SetTitle("RealEnergy");

  int divisions = 50;
  TH2F *w_z0z1 = new TH2F("w_z0z1","w_z0z1",divisions,-7.5,7.5,divisions,-7.5,7.5);
  w_z0z1->GetXaxis()->SetTitle("z0");
  w_z0z1->GetYaxis()->SetTitle("z1");

  TH2F *u_z0z1 = new TH2F("u_z0z1","u_z0z1",divisions,-7.5,7.5,divisions,-7.5,7.5);
  u_z0z1->GetXaxis()->SetTitle("z0");
  u_z0z1->GetYaxis()->SetTitle("z1");

  TH2F *delta_z01 = new TH2F("delta_z01","delta_z01",1000,-15,15,1000,-15,15);
  delta_z01->SetTitle("Delta z if 0->1 (z_reco - z_real)");
  delta_z01->GetXaxis()->SetTitle("delta z0 [mm]");
  delta_z01->GetYaxis()->SetTitle("delta z1 [mm]");

  TH2F *delta_z10 = new TH2F("delta_z10","delta_z10",1000,-15,15,1000,-15,15);
  delta_z10->SetTitle("Delta z if 1->0 (z_reco - z_real)");
  delta_z10->GetXaxis()->SetTitle("delta z0 [mm]");
  delta_z10->GetYaxis()->SetTitle("delta z1 [mm]");


  std::vector<float> z0_v;
  std::vector<float> z1_v;
  std::vector<float> wDiff;
  std::vector<float> uDiff;
  std::vector<float> vDiff;
  std::vector<float> wuDiff;
  std::vector<float> wvDiff;
  std::vector<float> wuvDiff;

  std::vector<float> z0_v_10;
  std::vector<float> z1_v_10;
  std::vector<float> wDiff_10;
  std::vector<float> uDiff_10;
  std::vector<float> vDiff_10;
  std::vector<float> wuDiff_10;
  std::vector<float> wvDiff_10;
  std::vector<float> wuvDiff_10;


  // float minPointZ0 = 0.0;
  // float minPointZ1 = 0.0;
  // float minPointV = 0.0;


  //----------------------------------------------------------------------//
  //                                                                      //
  //                        LOOP ON EVENTS                                //
  //                                                                      //
  //----------------------------------------------------------------------//
  long int counter = 0;
  int nEntries = tree->GetEntries();
  std::cout << "nEntries = " << nEntries << std::endl;
  // int cEvent = 5197;
  // int cEvent = 400;
  // int cEvent = 398;

  // int cEvent = 3620;
  // start and stop events
  int startEvent = 0;
  int stopEvent = nEntries;
  if(!allAnalysis)
  {
    startEvent = cEvent;
    stopEvent = cEvent+1;
  }

  //find front mppcs and frame mppcs
  std::vector<int> frontMppcs;
  std::vector<int> frameMppcs;
  std::vector<int> relevantMppcs;
  //front mppcs
  int trigger0 = (int) (mppcI0 * sqrt(numOfCh) + mppcJ0);
  int trigger1 = (int) (mppcI1 * sqrt(numOfCh) + mppcJ1);
  frontMppcs.push_back(trigger0); //push first crystal mppc to the frontMppcs
  if(trigger0 != trigger1){ // if the second crystal is in front of another mppc...
    frontMppcs.push_back(trigger1); //... push also that mppc id into the frontMppcs vector
  }
  //frameMppcs (actually frame plus front)
  //frame of mppc0
  for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
  {
    if( (iMPPC >= 0) && (iMPPC <= sqrt(numOfCh)) )
    {
      for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
      {
        if( (jMPPC >= 0) && (jMPPC <= sqrt(numOfCh)) )
        {
          int thisChannel = iMPPC * sqrt(numOfCh) + jMPPC;
          // std::cout << "thisChannel " << thisChannel << std::endl;
          bool isThere = false;
          for(int a = 0; a < frameMppcs.size(); a++) //check if id is already there
          {
            // std::cout << "frameMppcs[a] " << frameMppcs[a] << " thisChannel " << thisChannel << std::endl;
            if(frameMppcs[a] == thisChannel)
              isThere = true;
          }
          frameMppcs.push_back(thisChannel);
        }
      }
    }
  }
  //frame of mppc1
  for(int iMPPC = mppcI1 -1; iMPPC < mppcI1+2; iMPPC++)
  {
    if( (iMPPC >= 0) && (iMPPC <= sqrt(numOfCh)) )
    {
      for(int jMPPC = mppcJ1-1; jMPPC < mppcJ1+2 ; jMPPC++)
      {
        if( (jMPPC >= 0) && (jMPPC <= sqrt(numOfCh)) )
        {
          int thisChannel = iMPPC * sqrt(numOfCh) + jMPPC;
          // std::cout << "thisChannel " << thisChannel << std::endl;
          bool isThere = false;
          for(int a = 0; a < frameMppcs.size(); a++) //check if id is already there
          {
            // std::cout << "frameMppcs[a] " << frameMppcs[a] << " thisChannel " << thisChannel << std::endl;
            if(frameMppcs[a] == thisChannel)
              isThere = true;
          }
          frameMppcs.push_back(thisChannel);
        }
      }
    }
  }
  //now the relevantMppcs. Intersection of the two frame groups and front group
  for(int iFront = 0; iFront < frontMppcs.size(); iFront++) //first the two front
  {
    relevantMppcs.push_back(frontMppcs[iFront]);
  }
  for(int iFrame = 0; iFrame < frameMppcs.size(); iFrame++) //run on all frames and add if it's not already there
  {
    bool isThere = false;
    for(int a = 0; a < relevantMppcs.size(); a++) //check if id is already there
    {
      if(relevantMppcs[a] == frameMppcs[iFrame]){
        isThere = true;
      }
    }
    if(!isThere){
      relevantMppcs.push_back(frameMppcs[iFrame]);
    }
  }

  //TESTING - use all channels. CAREFUL: then also u,v,w need to be calculated on all channels here (which is automatically the case
  // since they are calculated on the relevantMppcs) and in the calibraton file (at least for the w(z), that's sure).
  bool hardcodeMppcs = true; //TEMPORARY
  if(hardcodeMppcs)
  {
    relevantMppcs.clear();
    frontMppcs.clear();
    relevantMppcs.push_back(0);
    relevantMppcs.push_back(1);
    relevantMppcs.push_back(2);
    relevantMppcs.push_back(3);
    relevantMppcs.push_back(4);
    relevantMppcs.push_back(5);
    relevantMppcs.push_back(6);
    relevantMppcs.push_back(7);
    relevantMppcs.push_back(8);
    relevantMppcs.push_back(9);
    relevantMppcs.push_back(10);
    relevantMppcs.push_back(11);
    relevantMppcs.push_back(12);
    relevantMppcs.push_back(13);
    relevantMppcs.push_back(14);
    relevantMppcs.push_back(15);
    frontMppcs.push_back(5);
    frontMppcs.push_back(6);
  }

  for(int i = 0; i < relevantMppcs.size(); i++){
    std::cout << relevantMppcs[i] << " " ;
  }
  std::cout << std::endl;

  for(int i = 0; i < frontMppcs.size(); i++){
    std::cout << frontMppcs[i] << " " ;
  }
  std::cout << std::endl;

  for(int iEvent = startEvent; iEvent < stopEvent ; iEvent++)
  // for(int iEvent = 0; iEvent < nEntries ; iEvent++)
  // for(int iEvent = cEvent; iEvent < cEvent + 1 ; iEvent++)
  {
    tree->GetEvent(iEvent);

    // calculate the u,v,w that we measure (we already know about the totalEnergyDeposited. it is worth to mention here that in reality
    // the check on the energy deposited is necessarily less accurate, given the energy resolution of the detector, but we will deal with that later)
    // u,v,w calculated with only near channels

    // to implement the different mppc modality, we need to define the
    // front mppc and framemppc. Front are the 2 (or more) directly coupled to the crystals, frame the
    // frame around each front

    columsum = 0.0;
    rowsum = 0.0;
    total = 0.0;
    frontSum = 0.0;
    floodx = 0.0;
    floody = 0.0;
    floodz = 0.0;
    // std::cout << "triggerChannel " << triggerChannel << std::endl;
    //now scan the channels and accept them only if they are into the relevantMppcs. Also, sum the contribution of frontMppcs
    for(int i = 0; i < numOfCh ; i++)
    {
      // check if it's rlevant and front
      bool isRelevant = false;
      bool isFront = false;
      for(int a = 0; a < frontMppcs.size(); a++)
      {
        if(frontMppcs[a] == i)
          isFront = true;
      }
      for(int a = 0; a < relevantMppcs.size(); a++)
      {
        if(relevantMppcs[a] == i)
          isRelevant = true;
      }

      if(isRelevant)
      {
        columsum += detector[i]*ymppc[i];
        rowsum += detector[i]*xmppc[i];
        total += detector[i];
        // std::cout << i << " " << detector[i] << std::endl;;
      }
      if(isFront)
      {
        frontSum += detector[i];
      }
      // std::cout  << " ";
    }

    // std::cout << std::endl;
    floodx = rowsum/total;
    floody = columsum/total;
    // std::cout << "total " << total << std::endl;
    if(total != 0.0) floodz = ((double) (frontSum/total));
    // std::cout << "<-------------------------------------------------------------------"<<std::endl;
    flood->Fill(floodx,floody);
    // std::cout << "flood x,y,z " << floodx << "," << floody << "," << floodz << std::endl;

    // from this 2d plot we can easily select the limits in u,v for what we will consider interscatter events, then select energy cuts etc. for now, we want to check
    // "the model" for the selection will be on the simulation data (inaccessible in reality, of course)

    // "average" what happened into each crystal
    // stackingaction in g4 is last-in-first-out so we need to sort energyDeposition before we check what crystal was hit first
    // check crystal sequence, check for returning gammas, calculate an average xyz per crystal
    std::sort(energyDeposition->begin(), energyDeposition->end(), compareByTime);
    std::vector<std::vector < enDep > > separatedEnDep;

    int CurrentID = -1;
    float RealX = 0.0;
    float RealY = 0.0;
    float RealZ = 0.0;
    std::vector < enDep > CrystalEnDepCollection;
    for(int eEvent = 0; eEvent < energyDeposition->size(); eEvent++) //run on energy depositions and find in how many crystals energy was deposited
    {
      //this for cycles on all the endep events. it stores them (each one is a struct) in a std::vector of structs, until the gamma changes crystal
      //when we enter a new crystal, the std::vector of that crystal is pushed_back into (guess what?) a std::vector of std::vector and another std::vector is created
      // so to sumarize:
      // energy deposition event -> struct
      // all energy deposition events in a given crystal -> std::vector<struct> -> i call this a "collection" below
      // all crystals touched by the gamma -> std::vector< std:: vector <struct> > this would be the std::vector of "collections"
      //read the crystal where energy was deposited
      int cry = energyDeposition->at(eEvent).CrystalID;
      if(eEvent == 0) CurrentID = cry; //needed for the first step
      //create temp enDep variable and copy this eEvent into it
      enDep tempCrystalEnDep;
      tempCrystalEnDep = energyDeposition->at(eEvent);
      if(cry != CurrentID) // if this endep event is happening in a new crystal wrt the one before
      {
        separatedEnDep.push_back(CrystalEnDepCollection); // save the collection of this crystal into the std::vector of collections
        CrystalEnDepCollection.clear(); //clear this collection
        CurrentID = cry; //change the current id
      }

      CrystalEnDepCollection.push_back(tempCrystalEnDep); // save this enDep event into the collection of this crystal
      if(eEvent == energyDeposition->size() -1)
      {
        separatedEnDep.push_back(CrystalEnDepCollection);
      }
    }

    // now for each crystal average what happened inside
    // what happened in a crystal cannot be seen in split details by the detector so we need some useful
    // variables to gain the chance to filter the simulation dataset, if we want to understand
    // what's going on
    std::vector<avgCryEnergyDep> averageDepEvents;
    // now the en dep events are collected by crystal
    // run on each collection and find average
    for(int iColl = 0 ; iColl < separatedEnDep.size(); iColl++)
    {
      avgCryEnergyDep tempAvgEnDep;

      //initialize
      tempAvgEnDep.id = separatedEnDep.at(iColl).at(0).CrystalID;
      tempAvgEnDep.x = 0;
      tempAvgEnDep.y = 0;
      tempAvgEnDep.z = 0;
      tempAvgEnDep.time = 0;
      tempAvgEnDep.energy = 0;
      tempAvgEnDep.sx = 0;
      tempAvgEnDep.sy = 0;
      tempAvgEnDep.sz = 0;
      //now run on the energy deposition events in this crystal and calculate the averages
      for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
      {
        tempAvgEnDep.energy += separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.x += separatedEnDep.at(iColl).at(iEndep).DepositionX * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.y += separatedEnDep.at(iColl).at(iEndep).DepositionY * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.z += separatedEnDep.at(iColl).at(iEndep).DepositionZ * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
        tempAvgEnDep.time += separatedEnDep.at(iColl).at(iEndep).DepositionTime * separatedEnDep.at(iColl).at(iEndep).EnergyDeposited;
      }
      tempAvgEnDep.x = tempAvgEnDep.x / tempAvgEnDep.energy;
      tempAvgEnDep.y = tempAvgEnDep.y / tempAvgEnDep.energy;
      tempAvgEnDep.z = tempAvgEnDep.z / tempAvgEnDep.energy;
      tempAvgEnDep.time = tempAvgEnDep.time / tempAvgEnDep.energy;
      float varx = 0.0; // variance (needed for stdev afterwards)
      float vary = 0.0;
      float varz = 0.0;
      for(int iEndep = 0; iEndep < separatedEnDep.at(iColl).size(); iEndep++)
      {
        varx += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionX  - tempAvgEnDep.x,2)) / tempAvgEnDep.energy;
        vary += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionY  - tempAvgEnDep.y,2)) / tempAvgEnDep.energy;
        varz += (separatedEnDep.at(iColl).at(iEndep).EnergyDeposited * pow(separatedEnDep.at(iColl).at(iEndep).DepositionZ  - tempAvgEnDep.z,2)) / tempAvgEnDep.energy;
      }
      tempAvgEnDep.sx = sqrt(varx);
      tempAvgEnDep.sy = sqrt(vary);
      tempAvgEnDep.sz = sqrt(varz);

      //save into the std::vector of averages
      averageDepEvents.push_back(tempAvgEnDep);
    }

    //two crystals or more hit, but a crystal is hit more than once (like, 28 -> 36 -> 28)
    std::vector <int> checkReturning; // this std::vector has same size of averageDepEvents if there is no returning, smaller if there is
    int returned = 0;
    bool isReturned = false;
    for (int iAvg = 0 ; iAvg < averageDepEvents.size() ; iAvg++)
    {
      int cry = averageDepEvents.at(iAvg).id;
      bool sameID = false;
      for(int i = 0 ; i < checkReturning.size(); i++)
      {
        if(cry == checkReturning[i])
        {
          sameID = true;
          isReturned = true;
          returned++;
        }
      }
      if(!sameID) checkReturning.push_back(cry);
    }
    globalReturned += returned;

    // now let's start filtering the dataset.
    // we start by default
    bool isCandidate = true;
    if(totalEnergyDeposited < 0.510) //filter out events with non total depostion
    {
      isCandidate = false; // cut on totalEnergyDeposited in this event
      // and count
      if(averageDepEvents.size() == 1) singleCounter++;
      if(averageDepEvents.size() == 2) doubleCounter++;
      if(averageDepEvents.size() == 3) tripleCounter++;
      if(averageDepEvents.size() > 3)  multipleCounter++;
    }

    if(averageDepEvents.size() != 2) //filter out events where energy deposition is in 1 or more than 2 crystals
    {
      isCandidate = false;
    }
    else // restric only to two specific crystals
    {
      if(averageDepEvents[0].id != cry0N) isCandidate = false;
      if(averageDepEvents[1].id != cry1N) isCandidate = false;
    }

    if(isCandidate) //let the fun begin
    {
      //check
      // std::cout << totalEnergyDeposited << " " << averageDepEvents.size() << " " << averageDepEvents[0].id << " " <<  averageDepEvents[1].id << std::endl;
      foundCandidate++;
      //take the measurables
      double u = floodx;
      double v = floody;
      double w = floodz;

      std::cout << "u-v-w = " <<  u << " " << v << " " << w << std::endl;

      //first approximation, beta = 0. Beta is gonna be tricky in reality. in this simulation is actually basically 0
      double beta = 0.0;
      double d = 1.6; //distance between crystals centers, in mm
      double en = 0.511; // total energy deposited in MeV. here again, in this dataset is true, in reality the energy resolution will make this tricky
      //check
      // std::cout.precision(3);
      Float_t source[3];
      source [0] = -1.6;  //old source position
      source [1] = +1.6;
      source [2] = -208.11;
      if(isNewDataset){
        source [0] = sourcex;
        source [1] = sourcey;
        source [2] = sourcez;
      }
      std::cout << "Source = " <<  source[0] << "," << source[1] << "," << source[2] << std::endl;
      Float_t realP0[3] = {averageDepEvents[0].x,averageDepEvents[0].y,averageDepEvents[0].z};
      Float_t realP1[3] = {averageDepEvents[1].x,averageDepEvents[1].y,averageDepEvents[1].z};
      Float_t P0[3] = {cry0x,cry0y,averageDepEvents[0].z};
      Float_t P1[3] = {cry1x,cry1y,averageDepEvents[1].z};
      Float_t P0_a[3];
      Float_t P1_a[3];
      Float_t P0_b[3];
      Float_t P1_b[3];
      Float_t E0_a;
      Float_t E0_b;
      Float_t E1_a;
      Float_t E1_b;

      if(fabs(averageDepEvents[0].energy - averageDepEvents[1].energy) > 0.0) // useless like this, just in case we want to investigate...
      {
        // std::cout << angle3D(source,realP0,realP1) << " " << angle3D(source,P0,P1) << " " << std::endl;
        // std::cout << iEvent << " " << 1000*averageDepEvents[1].energy << " " << 1000*scatteredGammaEnergy(en,source,realP0,realP1) << " " << 1000*scatteredGammaEnergy(en,source,P0,P1) << std::endl;
        RealEnergy_RealCoord->Fill(scatteredGammaEnergy(en,source,realP0,realP1),averageDepEvents[1].energy);
        RealEnergy_DetectorCoord->Fill(scatteredGammaEnergy(en,source,P0,P1),averageDepEvents[1].energy);

        //identify the event
        std::cout << "|------------------------------------------------|" << std::endl;
        std::cout << "|            EVENT "<< iEvent <<  "                         |" << std::endl;
        std::cout << "|------------------------------------------------|" << std::endl;
        // std::cout << sqrt(pow(averageDepEvents[0].x - cry0centerX,2) +pow(averageDepEvents[0].y - cry0centerY,2)) << " " << sqrt(pow(averageDepEvents[1].x - cry1centerX,2) +pow(averageDepEvents[1].y - cry1centerY,2)) << std::endl;


        Float_t startZ0;
        Float_t stopZ0;
        Float_t startZ1;
        Float_t stopZ1;
        Float_t stepZ0;
        Float_t stepZ1;
        Float_t range = 6.0;
        // startZ0 = averageDepEvents[0].z - range;
        // stopZ0  = averageDepEvents[0].z + range;
        // startZ1 = averageDepEvents[1].z - range;
        // stopZ1  = averageDepEvents[1].z + range;
        startZ0 = startZ1 = -7.4;  //crystal limits. if you go to +/-7.5 it crashes. why?
        stopZ0 = stopZ1 = 7.4;
        // minPointZ0 = averageDepEvents[0].z;
        // minPointZ1 = averageDepEvents[1].z;
        stepZ0  = (stopZ0 - startZ0)/divisions;
        stepZ1  = (stopZ1 - startZ1)/divisions;
        Float_t zValue[2];
        Float_t eValue[2];
        Float_t incidentEnergy = 0.511;
        Float_t point0[3];
        Float_t point1[3];
        Float_t recoW;
        Float_t u0,u1;
        Float_t v0,v1;
        Float_t recoU;
        Float_t recoV;
        Float_t totalProb_1_0;
        Float_t totalProb_0_1;
        Float_t min01 = INFINITY;
        Float_t min10 = INFINITY;
        Float_t z0Best_01,z1Best_01,z0Best_10,z1Best_10;
        // std::cout << "CASE A - 0-> 1 " << std::endl;

        if(performAnalysis)
        {
          Float_t RecoEn01[2],RecoEn10[2];
          for(int iDiv = 0 ; iDiv < divisions ; iDiv++)
          {
            for(int jDiv = 0 ; jDiv < divisions ; jDiv++)
            {
              // now first, case A
              // try two values, z0 and z1
              zValue[0] = startZ0 + stepZ0*iDiv;
              zValue[1] = startZ1 + stepZ1*jDiv;

              //point 0 is always in crystal 0
              point0[0] = cry0centerX;
              point0[1] = cry0centerY;
              // point0[0] = averageDepEvents[0].x;
              // point0[1] = averageDepEvents[0].y;
              point0[2] = zValue[0];

              point1[0] = cry1centerX;
              point1[1] = cry1centerY;
              // point1[0] = averageDepEvents[1].x;
              // point1[1] = averageDepEvents[1].y;
              point1[2] = zValue[1];

              //calculate first the values of E for these 2 z
              eValue[0] = incidentEnergy - scatteredGammaEnergy(incidentEnergy,source,point0,point1);
              eValue[1] = scatteredGammaEnergy(incidentEnergy,source,point0,point1);

              //calc w from z0 and z1
              recoW = averageW(eValue[0],eValue[1],wzgraph0,zValue[0],wzgraph1,zValue[1]);

              //calc u0 and u1 from e0,e1 and w0,w1
              Float_t charge0 = eValue[0] * correction0;
              Float_t charge1 = eValue[1] * correction1;

              //now the sum to find u0 and u1 have to be modified to use always the relevantMppcs
              if(interpolationWay == 2)
              {
                u0 = computeU(gd0,xmppc,charge0,computeW(wzgraph0,zValue[0]),relevantMppcs,numOfCh);
                u1 = computeU(gd1,xmppc,charge1,computeW(wzgraph1,zValue[1]),relevantMppcs,numOfCh);
              }
              if(interpolationWay == 2)
              {
                v0 = computeV(gd0,ymppc,charge0,computeW(wzgraph0,zValue[0]),relevantMppcs,numOfCh);
                v1 = computeV(gd1,ymppc,charge1,computeW(wzgraph1,zValue[1]),relevantMppcs,numOfCh);
              }

              //TEMPORARY HACK
              //skip the points where one the interpolation failed
              bool skip = false;
              if(u0 == 0 | u1 == 0 | v0 == 0 | v1 == 0)
              {
                skip = true;
              }

              if(!skip)
              {
                recoU = averageU(eValue[0],eValue[1],u0,u1);
                recoV = averageV(eValue[0],eValue[1],v0,v1);
                z0_v.push_back(zValue[0]);
                z1_v.push_back(zValue[1]);
                wDiff.push_back(fabs(w - recoW));
                uDiff.push_back(fabs(u - recoU)/12.0);
                vDiff.push_back(fabs(v - recoV)/12.0);
                wuvDiff.push_back(sqrt(pow((u - recoU)/12.0,2)+pow(w - recoW,2)+pow((v -recoV)/(12.0),2)));
                wuDiff.push_back(sqrt(pow((u - recoU)/12.0,2)+pow(w - recoW,2)));
                wvDiff.push_back(sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)));
                if(min01 > sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)))
                {
                  min01 = sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2));
                  z0Best_01 = zValue[0];
                  z1Best_01 = zValue[1];
                  RecoEn01[0] = eValue[0];
                  RecoEn01[1] = eValue[1];
                }


                //intersection source and plane z = -7.5
                Float_t v_sp0[3] = {source[0]-point0[0],source[1]-point0[1],source[2]-point0[2]};
                Float_t y_int = ((-7.5 - point0[2]) / (v_sp0[2])) * v_sp0[1] + point0[1];
                Float_t x_int = point0[0] + ( (y_int - point0[1])/(v_sp0[1]) )*v_sp0[0];
                Float_t intersection0[3] = {x_int,y_int,-7.5};
                // std::cout <<  "Inters.\t" << intersection0[0] << " " << intersection0[1] << " " << intersection0[2] << std::endl;
                Float_t disTravel1_0 = distance3D(intersection0[0],intersection0[1],intersection0[2],point0[0], point0[1], point0[2]);
                // std::cout << "point0\t" << point0[0] << " " << point0[1] << " " << point0[2] << std::endl;
                // std::cout <<  "disTravel1_0\t" << disTravel1_0 << std::endl;
                Float_t p_travel1 = simTravel1(disTravel1_0, lambdaLYSO->Eval(0.511));
                // std::cout <<  "p_travel1\t" << p_travel1 << std::endl;
                // std::cout <<  "angle\t\t" << angle3D(source,point0,point1) << std::endl;
                Float_t p_compton = simCompton(angle3D(source,point0,point1));
                // std::cout <<  "p_compton\t" << p_compton << std::endl;
                Float_t comptonPhotoelDistance = distance3D(point0[0], point0[1], point0[2],point1[0], point1[1], point1[2]);
                // std::cout <<  "comptonPhotoelDistance\t" << comptonPhotoelDistance << std::endl;
                Float_t p_travel2 = simTravel1(comptonPhotoelDistance, lambdaLYSO->Eval(eValue[1]));
                // std::cout <<  "p_travel2\t" << p_travel2 << std::endl;
                Float_t p_photoelect = simPhotoelectric(photoelectricCrossSectionLYSO->Eval(eValue[1]));
                // std::cout <<  "p_photoelect\t" << p_photoelect << std::endl;

                // std::cout << std::endl;
                if(min01 > sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)))
                  totalProb_0_1  = p_travel1*p_compton*p_travel2*p_photoelect;
              }
            }
          }

          //compute delta z in the hypo that 0->1
          delta_z01->Fill(z0Best_01-averageDepEvents[0].z,z1Best_01-averageDepEvents[1].z);
          //TEMPORARY
          if(fabs(z0Best_01-averageDepEvents[0].z) < 0.1 && fabs(z1Best_01-averageDepEvents[1].z) < 0.1)
          {
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<------ " << iEvent << "---" <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
            std::cout << "<---------------------- " <<std::endl;
          }

          bool doCaseB = true; // was useful before, now not really but let's leave the option here...
          if(doCaseB)
          {
            std::cout << std::endl;
            std::cout << "CASE B - 1->0" << std::endl;

            // startZ0 = averageDepEvents[0].z - range;
            // stopZ0  = averageDepEvents[0].z + range;
            // startZ1 = averageDepEvents[1].z - range;
            // stopZ1  = averageDepEvents[1].z + range;
            startZ0 = startZ1 = -7.4;
            stopZ0 = stopZ1 = 7.4;
            stepZ0  = (stopZ0 - startZ0)/divisions;
            stepZ1  = (stopZ1 - startZ1)/divisions;
            for(int iDiv = 0 ; iDiv < divisions ; iDiv++)
            {
              for(int jDiv = 0 ; jDiv < divisions ; jDiv++)
              {
                // case B
                // try two values, z0 and z1
                zValue[0] = startZ0 + stepZ0*iDiv;
                zValue[1] = startZ1 + stepZ1*jDiv;

                //point 0 is always in crystal 0
                point0[0] = cry0centerX;
                point0[1] = cry0centerY;
                // point0[0] = averageDepEvents[0].x;
                // point0[1] = averageDepEvents[0].y;
                point0[2] = zValue[0];

                point1[0] = cry1centerX;
                point1[1] = cry1centerY;
                // point1[0] = averageDepEvents[1].x;
                // point1[1] = averageDepEvents[1].y;
                point1[2] = zValue[1];

                // std::cout << zValue[0] << " " << zValue[1] << std::endl;
                //calculate first the values of E for these 2 z
                eValue[0] = scatteredGammaEnergy(incidentEnergy,source,point1,point0);
                eValue[1] = incidentEnergy - scatteredGammaEnergy(incidentEnergy,source,point1,point0);

                //calc w from z0 and z1
                recoW = averageW(eValue[0],eValue[1],wzgraph0,zValue[0],wzgraph1,zValue[1]);

                //calc u0 and u1 from e0,e1 and w0,w1
                Float_t charge0 = eValue[0] * correction0;
                Float_t charge1 = eValue[1] * correction1;
                // std::cout << "eValue/charge\t" << eValue[0] << " " << charge0 << " " << eValue[1] <<  " " << charge1 << std::endl;

                if(interpolationWay == 2)
                {
                  u0 = computeU(gd0,xmppc,charge0,computeW(wzgraph0,zValue[0]),relevantMppcs,numOfCh);
                  u1 = computeU(gd1,xmppc,charge1,computeW(wzgraph1,zValue[1]),relevantMppcs,numOfCh);
                }

                // Float_t v0,v1;
                if(interpolationWay == 2)
                {
                  v0 = computeV(gd0,ymppc,charge0,computeW(wzgraph0,zValue[0]),relevantMppcs,numOfCh);
                  v1 = computeV(gd1,ymppc,charge1,computeW(wzgraph1,zValue[1]),relevantMppcs,numOfCh);
                }
                //TEMPORARY HACK
                //skip the points where one the interpolation failed
                bool skip = false;
                if(u0 == 0 | u1 == 0 | v0 == 0 | v1 == 0)
                {
                  skip = true;
                }

                if(!skip)
                {

                  recoU = averageU(eValue[0],eValue[1],u0,u1);
                  recoV = averageV(eValue[0],eValue[1],v0,v1);

                  z0_v_10.push_back(zValue[0]);
                  z1_v_10.push_back(zValue[1]);
                  wDiff_10.push_back(fabs(w - recoW));
                  uDiff_10.push_back(fabs(u - recoU)/12.0);
                  vDiff_10.push_back(fabs(v - recoV)/12.0);
                  wuvDiff_10.push_back(sqrt(pow((u - recoU)/12.0,2)+pow(w - recoW,2)+pow((v -recoV)/(12.0),2)));
                  wuDiff_10.push_back(sqrt(pow((u - recoU)/12.0,2)+pow(w - recoW,2)));
                  wvDiff_10.push_back(sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)));
                  if(min10 > sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)))
                  {
                    min10 = sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2));
                    z0Best_10 = zValue[0];
                    z1Best_10 = zValue[1];
                    RecoEn10[0] = eValue[0];
                    RecoEn10[1] = eValue[1];
                  }

                  Float_t v_sp1[3] = {source[0]-point1[0],source[1]-point1[1],source[2]-point1[2]};
                  Float_t y_int = ((-7.5 - point1[2]) / (v_sp1[2])) * v_sp1[1] + point1[1];
                  Float_t x_int = point1[0] + ( (y_int - point1[1])/(v_sp1[1]) )*v_sp1[0];
                  Float_t intersection1[3] = {x_int,y_int,-7.5};
                  // std::cout <<  "Inters.\t" << intersection1[0] << " " << intersection1[1] << " " << intersection1[2] << std::endl;

                  Float_t disTravel1_0 = distance3D(intersection1[0],intersection1[1],intersection1[2],point1[0], point1[1], point1[2]);
                  // std::cout << "point1\t" << point1[0] << " " << point1[1] << " " << point1[2] << std::endl;
                  // std::cout <<  "disTravel1_0\t" << disTravel1_0 << std::endl;
                  // disTravel1_0 = distance3D(0 + paramT0*(point1[0] - 0),0 + paramT0*(point1[1] - 0),-7.5, point1[0], point1[1], point1[2]);
                  Float_t p_travel1 = simTravel1(disTravel1_0, lambdaLYSO->Eval(0.511));
                  // std::cout <<  "p_travel1\t" << p_travel1 << std::endl;
                  // std::cout <<  "angle\t\t" << angle3D(source,point1,point0) << std::endl;
                  Float_t p_compton = simCompton(angle3D(source,point1,point0));
                  // std::cout <<  "p_compton\t" << p_compton << std::endl;
                  Float_t comptonPhotoelDistance = distance3D(point1[0], point1[1], point1[2],point0[0], point0[1], point0[2]);
                  // std::cout <<  "comptonPhotoelDistance\t" << comptonPhotoelDistance << std::endl;
                  Float_t p_travel2 = simTravel1(comptonPhotoelDistance, lambdaLYSO->Eval(eValue[0]));
                  // std::cout <<  "p_travel2\t" << p_travel2 << std::endl;
                  Float_t p_photoelect = simPhotoelectric(photoelectricCrossSectionLYSO->Eval(eValue[0]));
                  // std::cout <<  "p_photoelect\t" << p_photoelect << std::endl;
                  if(min10 > sqrt(pow((v - recoV)/12.0,2)+pow(w - recoW,2)))
                    totalProb_1_0 = p_travel1*p_compton*p_travel2*p_photoelect;
                }
              }
            }
          }

          delta_z10->Fill(z0Best_10-averageDepEvents[0].z,z1Best_10-averageDepEvents[1].z);
          std::cout << std::endl;
          std::cout << "Real Z " << averageDepEvents[0].z << " " << averageDepEvents[1].z << std::endl;
          std::cout << "CASE 0->1, z " << z0Best_01 << " " << z1Best_01 << std::endl;
          std::cout << "CASE 1->0, z " << z0Best_10 << " " << z1Best_10 << std::endl;
          std::cout << std::endl;
          std::cout << "Real E " << averageDepEvents[0].energy << " " << averageDepEvents[1].energy << std::endl;
          std::cout << "CASE 0->1, E " << RecoEn01[0] << " " << RecoEn01[1] << std::endl;
          std::cout << "CASE 1->0, E " << RecoEn10[0] << " " << RecoEn10[1] << std::endl;
          std::cout << std::endl;
          std::cout <<  "P(CASE_A)\tP(CASE_B)\tP(CASE_A)/P(CASE_B)"<< std::endl;
          std::cout << totalProb_0_1 << "\t\t" << totalProb_1_0 << "\t\t" << totalProb_0_1/totalProb_1_0 << std::endl;
          histoProb->Fill(totalProb_0_1/totalProb_1_0);
          std::cout << "|--------------------------------------|"<< std::endl;
        }
        else // don't perform analysis, just output the event values
        {

          eValue[0] = incidentEnergy - scatteredGammaEnergy(incidentEnergy,source,realP0,realP1);
          eValue[1] = scatteredGammaEnergy(incidentEnergy,source,realP0,realP1);


          std::cout << std::endl;
          std::cout << "Event " << iEvent << std::endl;
          std::cout << "cry = " << averageDepEvents[0].id << " "<< averageDepEvents[1].id << std::endl;
          std::cout << "Z = " <<  averageDepEvents[0].z << " " << averageDepEvents[1].z << std::endl;
          std::cout << "E = "<< averageDepEvents[0].energy << " " << averageDepEvents[1].energy <<  std::endl;
          std::cout << "reco E = " << eValue[0] << " " << eValue[1] << std::endl;
          Float_t axis[3] = {averageDepEvents[0].x,averageDepEvents[0].y,sourcez};
          if(!isNewDataset) axis[2] = -208.11;
          Float_t firstImpact[3] = {averageDepEvents[0].x,averageDepEvents[0].y,averageDepEvents[0].z};
          if(isNewDataset) std::cout << "source = " << sourcex << " " << sourcey << " " << sourcez << std::endl;
          std::cout << "angle = " << (angle3D(source,point0,axis) / TMath::Pi())*180.0  << " deg" << std::endl;
        }
      }
    }

    counter++;

    int perc = ((100*counter)/nEntries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      // std::cout << "\r";
      // std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }

  std::cout << std::endl;
  std::cout << "1 cry [511 KeV deposition] events = "    << singleCounter << std::endl;
  std::cout << "2 cry [511 KeV deposition] events = "    << doubleCounter << std::endl;
  std::cout << "3 cry [511 KeV deposition] events = "    << tripleCounter << std::endl;
  std::cout << "Multi cry [511 KeV deposition] events = "<< multipleCounter << std::endl;
  std::cout << "Returned = " << globalReturned << std::endl;
  std::cout << "Candidates = "<< foundCandidate << std::endl;

  TCanvas *c_uDiff = new TCanvas();
  TGraph2D *g_uDiff = new TGraph2D();
  TCanvas *c_vDiff = new TCanvas();
  TGraph2D *g_vDiff = new TGraph2D();
  TCanvas *c_wDiff = new TCanvas();
  TGraph2D *g_wDiff = new TGraph2D();
  TCanvas *c_wuDiff = new TCanvas();
  TGraph2D *g_wuDiff = new TGraph2D();
  TCanvas *c_wvDiff = new TCanvas();
  TGraph2D *g_wvDiff = new TGraph2D();
  TCanvas *c_wuvDiff = new TCanvas();
  TGraph2D *g_wuvDiff = new TGraph2D();
  TCanvas *c_uDiff_10 = new TCanvas();
  TGraph2D *g_uDiff_10 = new TGraph2D();
  TCanvas *c_vDiff_10 = new TCanvas();
  TGraph2D *g_vDiff_10 = new TGraph2D();
  TCanvas *c_wDiff_10 = new TCanvas();
  TGraph2D *g_wDiff_10 = new TGraph2D();
  TCanvas *c_wuDiff_10 = new TCanvas();
  TGraph2D *g_wuDiff_10 = new TGraph2D();
  TCanvas *c_wuvDiff_10 = new TCanvas();
  TGraph2D *g_wuvDiff_10 = new TGraph2D();
  TCanvas *c_wvDiff_10 = new TCanvas();
  TGraph2D *g_wvDiff_10 = new TGraph2D();

  std::string outFile = "compton.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  comptonCrossSectionLYSO->Write();
  photoelectricCrossSectionLYSO->Write();
  lambdaLYSO->Write();
  histoProb->Write();
  wzgraph0->Write();
  wzgraph1->Write();

  if(performAnalysis && !allAnalysis)
  {
    std::cout << "Writing TGraph2Ds..." << std::endl;

    c_uDiff = new TCanvas("c_uDiff","c_uDiff",1200,800);
    g_uDiff = new TGraph2D(uDiff.size(),&z0_v[0],&z1_v[0],&uDiff[0]);
    c_uDiff->cd();
    g_uDiff->SetName("c_uDiff");
    g_uDiff->SetTitle("c_uDiff");
    g_uDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_uDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_uDiff->GetZaxis()->SetTitle("u diff");
    g_uDiff->Draw("surf1");
    c_uDiff->Write();


    c_vDiff = new TCanvas("c_vDiff","c_vDiff",1200,800);
    g_vDiff = new TGraph2D(vDiff.size(),&z0_v[0],&z1_v[0],&vDiff[0]);
    c_vDiff->cd();
    g_vDiff->SetName("c_vDiff");
    g_vDiff->SetTitle("c_vDiff");
    g_vDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_vDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_vDiff->GetZaxis()->SetTitle("v diff");
    g_vDiff->Draw("surf1");
    c_vDiff->Write();


    c_wDiff = new TCanvas("c_wDiff","c_wDiff",1200,800);
    g_wDiff = new TGraph2D(wDiff.size(),&z0_v[0],&z1_v[0],&wDiff[0]);
    c_wDiff->cd();
    g_wDiff->SetName("c_wDiff");
    g_wDiff->SetTitle("c_wDiff");
    g_wDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_wDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_wDiff->GetZaxis()->SetTitle("w diff");
    g_wDiff->Draw("surf1");
    c_wDiff->Write();


    c_wuDiff = new TCanvas("c_wuDiff","c_wuDiff",1200,800);
    g_wuDiff = new TGraph2D(wuDiff.size(),&z0_v[0],&z1_v[0],&wuDiff[0]);
    c_wuDiff->cd();
    g_wuDiff->SetName("c_wuDiff");
    g_wuDiff->SetTitle("c_wuDiff");
    g_wuDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_wuDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_wuDiff->GetZaxis()->SetTitle("wu diff");
    g_wuDiff->Draw("surf1");
    c_wuDiff->Write();


    c_wvDiff = new TCanvas("c_wvDiff","c_wvDiff",1200,800);
    g_wvDiff = new TGraph2D(wvDiff.size(),&z0_v[0],&z1_v[0],&wvDiff[0]);
    c_wvDiff->cd();
    g_wvDiff->SetName("c_wvDiff");
    g_wvDiff->SetTitle("c_wvDiff");
    g_wvDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_wvDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_wvDiff->GetZaxis()->SetTitle("wv diff");
    g_wvDiff->Draw("surf1");
    c_wvDiff->Write();


    c_wuvDiff = new TCanvas("c_wuvDiff","c_wuvDiff",1200,800);
    g_wuvDiff = new TGraph2D(wuvDiff.size(),&z0_v[0],&z1_v[0],&wuvDiff[0]);
    c_wuvDiff->cd();
    g_wuvDiff->SetName("c_wuvDiff");
    g_wuvDiff->SetTitle("c_wuvDiff");
    g_wuvDiff->GetXaxis()->SetTitle("z0 [mm]");
    g_wuvDiff->GetYaxis()->SetTitle("z1 [mm]");
    g_wuvDiff->GetZaxis()->SetTitle("wvu diff");
    g_wuvDiff->Draw("surf1");
    c_wuvDiff->Write();

    // cry 1 -> 0
    c_uDiff_10 = new TCanvas("c_uDiff_10","c_uDiff_10",1200,800);
    g_uDiff_10 = new TGraph2D(uDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&uDiff_10[0]);
    c_uDiff_10->cd();
    g_uDiff_10->SetName("c_uDiff_10");
    g_uDiff_10->SetTitle("c_uDiff_10");
    g_uDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_uDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_uDiff_10->GetZaxis()->SetTitle("u diff");
    g_uDiff_10->Draw("surf1");
    c_uDiff_10->Write();


    c_vDiff_10 = new TCanvas("c_vDiff_10","c_vDiff_10",1200,800);
    g_vDiff_10 = new TGraph2D(vDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&vDiff_10[0]);
    c_vDiff_10->cd();
    g_vDiff_10->SetName("c_vDiff_10");
    g_vDiff_10->SetTitle("c_vDiff_10");
    g_vDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_vDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_vDiff_10->GetZaxis()->SetTitle("v diff");
    g_vDiff_10->Draw("surf1");
    c_vDiff_10->Write();


    c_wDiff_10 = new TCanvas("c_wDiff_10","c_wDiff_10",1200,800);
    g_wDiff_10 = new TGraph2D(wDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&wDiff_10[0]);
    c_wDiff_10->cd();
    g_wDiff_10->SetName("c_wDiff_10");
    g_wDiff_10->SetTitle("c_wDiff_10");
    g_wDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_wDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_wDiff_10->GetZaxis()->SetTitle("w diff");
    g_wDiff_10->Draw("surf1");
    c_wDiff_10->Write();

    c_wuDiff_10 = new TCanvas("c_wuDiff_10","c_wuDiff_10",1200,800);
    g_wuDiff_10 = new TGraph2D(wuDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&wuDiff_10[0]);
    c_wuDiff_10->cd();
    g_wuDiff_10->SetName("c_wuDiff_10");
    g_wuDiff_10->SetTitle("c_wuDiff_10");
    g_wuDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_wuDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_wuDiff_10->GetZaxis()->SetTitle("wu diff");
    g_wuDiff_10->Draw("surf1");
    c_wuDiff_10->Write();


    c_wvDiff_10 = new TCanvas("c_wvDiff_10","c_wvDiff_10",1200,800);
    g_wvDiff_10 = new TGraph2D(wuvDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&wuvDiff_10[0]);
    c_wvDiff_10->cd();
    g_wvDiff_10->SetName("c_wvDiff_10");
    g_wvDiff_10->SetTitle("c_wvDiff_10");
    g_wvDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_wvDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_wvDiff_10->GetZaxis()->SetTitle("wv diff");
    g_wvDiff_10->Draw("surf1");
    c_wvDiff_10->Write();

    c_wuvDiff_10 = new TCanvas("c_wuvDiff_10","c_wuvDiff_10",1200,800);
    g_wuvDiff_10 = new TGraph2D(wvDiff_10.size(),&z0_v_10[0],&z1_v_10[0],&wvDiff_10[0]);
    c_wuvDiff_10->cd();
    g_wuvDiff_10->SetName("c_wuvDiff_10");
    g_wuvDiff_10->SetTitle("c_wuvDiff_10");
    g_wuvDiff_10->GetXaxis()->SetTitle("z0 [mm]");
    g_wuvDiff_10->GetYaxis()->SetTitle("z1 [mm]");
    g_wuvDiff_10->GetZaxis()->SetTitle("wvu diff");
    g_wuvDiff_10->Draw("surf1");
    c_wuvDiff_10->Write();
  }

  flood->Write();
  if(performAnalysis)
  {
    RealEnergy_RealCoord->Write();
    RealEnergy_DetectorCoord->Write();
  }

  delta_z01->Write();
  delta_z10->Write();

  if(performAnalysis)
  {
    if(interpolationWay == 2)
    {
      TCanvas *C_spectrum;
      for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
      {
        for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
        {
          C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
          C_spectrum->SetName(camp0[iMPPC][jMPPC]->GetName());
          C_spectrum->cd();
          camp0[iMPPC][jMPPC]->Draw("p0");
          // f28[iMPPC][jMPPC]->Draw("surf same");
          C_spectrum->Write();
          delete C_spectrum;


          C_spectrum = new TCanvas("C_spectrum","C_spectrum",800,800);
          C_spectrum->SetName(camp1[iMPPC][jMPPC]->GetName());
          C_spectrum->cd();
          camp1[iMPPC][jMPPC]->Draw("p0");
          // f28[iMPPC][jMPPC]->Draw("surf same");
          C_spectrum->Write();
          delete C_spectrum;
          // g28[iMPPC][jMPPC]->Write();
        }
      }
    }
  }

  fOut->Close();

  return 0;
}




//distance between 2 points in 3d space
Float_t distance3D(Float_t ax, Float_t ay, Float_t az, Float_t bx, Float_t by, Float_t bz)
{
  Float_t v[3] = {bx-ax, by-ay, bz-az};
  Float_t vMod = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
  return vMod;
}

//angle between 3 points in 3d space
Float_t angle3D(Float_t a[3], Float_t b[3], Float_t c[3])
{
  Float_t v1[3] = {b[0]-a[0], b[1]-a[1], b[2]-a[2]};
  Float_t v2[3] = {c[0]-b[0], c[1]-b[1], c[2]-b[2]};
  Float_t v1Mod = sqrt(pow(v1[0],2) + pow(v1[1],2) + pow(v1[2],2));
  Float_t v2Mod = sqrt(pow(v2[0],2) + pow(v2[1],2) + pow(v2[2],2));
  Float_t dotProduct = v1[0] * v2[0] + v1[1] * v2[1] + v1[2] * v2[2];
  Float_t angle = acos(dotProduct/(v1Mod*v2Mod));
  return angle;
}

//compton scattered energy
// in our model, case A is when first crystal hit is crystal 0 (crystal 28 in this dataset)
// case B when first is crystal 1 (crystal 29 in this dataset). So, the args of this function will be
// energy = 0.511
// b = z_0^(A) or z_1^(B)
// c = z_1^(A) or z_0_(B)
Float_t scatteredGammaEnergy(Float_t energy,Float_t a[3],Float_t b[3],Float_t c[3])
{
  return energy/(2.0 - cos(angle3D(a,b,c)));
}

Float_t averageW(Float_t e0,Float_t e1,TGraph* w0,Float_t z0,TGraph* w1,Float_t z1)
{
  return ( (e0*computeW(w0,z0)) + (e1*computeW(w1,z1)) ) / (e0+e1);
}

Float_t averageU(Float_t e0,Float_t e1,Float_t u0,Float_t u1)
{
  return ( (e0*u0) + (e1*u1) ) / (e0+e1);
}

Float_t averageV(Float_t e0,Float_t e1,Float_t v0,Float_t v1)
{
  return ( (e0*v0) + (e1*v1) ) / (e0+e1);
}


Float_t computeW(TGraph* wz,Float_t z) // remember that in the ModuleCalibration z(w) have z goes from 0 to 15, here from -7.5 to 7.5
{
  // std::cout << "------- " << z+7.5 << " " << wz->Eval(z+7.5) << std::endl;
  return wz->Eval(-(z-7.5));
}

Float_t computeCrystalContribution(TF2* f2[4][4], Double_t position[16],Float_t charge,Float_t w)
{
  Float_t a = 0.0;
  Float_t sum = 0.0;
  for(int iMPPC = 0; iMPPC < 3; iMPPC++)
  {
    for(int jMPPC = 1; jMPPC < 4 ; jMPPC++)
    {

      int indexMPPC = iMPPC * 4 + jMPPC;
      // std::cout << "------ " <<  iMPPC << " " << jMPPC << " "<< indexMPPC << std::endl;
      a += f2[iMPPC][jMPPC]->Eval(w,charge) * position[indexMPPC];
      sum += f2[iMPPC][jMPPC]->Eval(w,charge);
    }
  }
  a = a/sum;
  return a;
}


// Float_t computeU(TGraphDelaunay ***gd,Double_t *xmppc,Float_t charge,Float_t w,int mppcI0, int mppcJ0, int numOfCh)
Float_t computeU(TGraphDelaunay ***gd,Double_t *xmppc,Float_t charge,Float_t w,std::vector<int>& relevantMppcs,int numOfCh)
{
  Float_t u = 0.0;
  Float_t sum = 0.0;
  // for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
  // {
  //   for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
  //   {
  for(int i = 0; i < relevantMppcs.size(); i++)
  {
    int iMPPC = ( (int) relevantMppcs[i]) / ( (int) sqrt(numOfCh) );
    int jMPPC = ( (int) relevantMppcs[i]) % ( (int) sqrt(numOfCh) );
    int indexMPPC = relevantMppcs[i];
    // std::cout << "relevantMppcs[i] = " << relevantMppcs[i] << std::endl;
    // std::cout << "iMPPC = " << iMPPC << std::endl;
    // std::cout << "jMPPC = " << jMPPC << std::endl;
    // std::cout << "indexMPPC = " << indexMPPC << std::endl;
    // int indexMPPC = iMPPC*sqrt(numOfCh) + jMPPC;
    u += gd[iMPPC][jMPPC]->ComputeZ(w,charge) * xmppc[indexMPPC];
    // std::cout << w << " "<< charge << " " << gd[iMPPC][jMPPC]->ComputeZ(w,charge) << " " << xmppc[indexMPPC] << std::endl;
    sum += gd[iMPPC][jMPPC]->ComputeZ(w,charge);
    // std::cout << indexMPPC << " " << u << " " << sum << std::endl;
  }
  //   }
  // }
  if(sum >0) u = u/sum;
  else u = 0;
  return u;
}

Float_t computeV(TGraphDelaunay ***gd,Double_t *ymppc,Float_t charge,Float_t w,std::vector<int>& relevantMppcs,int numOfCh)
{
  Float_t v = 0.0;
  Float_t sum = 0.0;
  // for(int iMPPC = mppcI0 -1; iMPPC < mppcI0+2; iMPPC++)
  // {
  //   for(int jMPPC = mppcJ0-1; jMPPC < mppcJ0+2 ; jMPPC++)
  //   {
  // std::cout << "relevantMppcs.size() " << relevantMppcs.size() << std::endl;
  // for(int i = 0; i < relevantMppcs.size(); i++){
  //   std::cout << relevantMppcs[i] << " " ;
  // }
  // std::cout << std::endl;
  for(int i = 0; i < relevantMppcs.size(); i++)
  {

    int iMPPC = ( (int) relevantMppcs[i]) / ( (int) sqrt(numOfCh) );
    int jMPPC = ( (int) relevantMppcs[i]) % ( (int) sqrt(numOfCh) );
    int indexMPPC = relevantMppcs[i];
      // int indexMPPC = iMPPC*sqrt(numOfCh) + jMPPC;

    v += gd[iMPPC][jMPPC]->ComputeZ(w,charge) * ymppc[indexMPPC];
    sum += gd[iMPPC][jMPPC]->ComputeZ(w,charge);
  }
    // }
  // }
  if(sum > 0) v = v/sum;
  else v = 0;
  return v;
}

int generateZ(Float_t* zValue)
{
  //for the moment stupid
  zValue[0] = 3.5;
  zValue[1] = -3.5;
  return 0;
}

Float_t simTravel1(Float_t pathLength, Float_t lambda511)
{
  Float_t probTravel1 = exp(-pathLength/lambda511);
  return probTravel1;
}

//compton effect function
Float_t simCompton(Float_t comptonAngle)
{
  Float_t finalEnergy = 0.511/(2-cos(comptonAngle));
  //std::cout << "energy calculated with angle: " << finalEnergy << std::endl;
  //
  Float_t probCompton = (2.0*3.1415)* sin(comptonAngle) * pow((finalEnergy/0.511),2)*(0.511/finalEnergy + finalEnergy/0.511 - pow(sin(comptonAngle),2));
  return probCompton;
}

//travel function used between compton and photoelectric effects
//energy is 511-energy deposited in first crystal because of compton effect, in keV
//pathLength and lambdaE in cm
Float_t simTravel2(Float_t pathLength, Float_t lambdaE)
{
  Float_t probTravel2 = exp(-pathLength/lambdaE);
  return probTravel2;
}

//photoelectric effect function
//energy is energy deposited in second crystal, in keV
Float_t simPhotoelectric(Float_t csPE)
{
  Float_t probPhotoelectric = csPE;
  return probPhotoelectric;
}
