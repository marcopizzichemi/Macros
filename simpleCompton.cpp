// compile with
// g++ -o ../build/simpleCompton simpleCompton.cpp `root-config --cflags --glibs`
// syntax


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

Float_t transmissionProbability(Float_t pathLength, Float_t lambda)
{
  Float_t probTravel1 = exp(-pathLength/lambda);
  return probTravel1;
};

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
};

//distance between 2 points in 3d space
Float_t distance3D(Float_t ax, Float_t ay, Float_t az, Float_t bx, Float_t by, Float_t bz)
{
  Float_t v[3] = {bx-ax, by-ay, bz-az};
  Float_t vMod = sqrt(pow(v[0],2) + pow(v[1],2) + pow(v[2],2));
  return vMod;
};

Float_t scatteredGammaEnergy(Float_t energy,Float_t a[3],Float_t b[3],Float_t c[3])
{
  return energy/(2.0 - cos(angle3D(a,b,c)));
};

//compton effect function
Float_t simCompton(Float_t comptonAngle)
{
  Float_t finalEnergy = 0.511/(2.0-cos(comptonAngle));
  //std::cout << "energy calculated with angle: " << finalEnergy << std::endl;
  //
  Float_t probCompton = (2.0*3.1415)* sin(comptonAngle) * pow((finalEnergy/0.511),2)*(0.511/finalEnergy + finalEnergy/0.511 - pow(sin(comptonAngle),2));
  return probCompton;
};
//photoelectric effect function
//energy is energy deposited in second crystal, in MeV
Float_t simPhotoelectric(Float_t csPE)
{
  Float_t probPhotoelectric = csPE;
  return probPhotoelectric;
}

double MeasureProb(double Em,double Et,double sigma)
{
  double factor1 = (1.0/(TMath::Sqrt(2.0*TMath::Pi()))*sigma);
  double factor2 = exp(- pow((Em-Et),2) / (2.0*pow(sigma,2)) );
  double measureprob = factor1*factor2;
  // if(measureprob == 0) return 0.0000000000000001;
  return measureprob;
}

int main (int argc, char** argv)
{
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
  lambdaLYSO->GetYaxis()->SetTitle("lambda [mm]");

  Float_t zValue[2];
  int divisions = 50;
  Float_t startZ0 = 0.0;
  Float_t startZ1 = 0.0;
  Float_t stopZ0 = 15.0;
  Float_t stopZ1 = 15.0;
  Float_t stepZ0  = (stopZ0 - startZ0)/divisions;
  Float_t stepZ1  = (stopZ1 - startZ1)/divisions;
  Float_t travel1_true; //true is always 0 to 1. probability for the 511 Kev gamma to reach point in 0
  Float_t travel1_false; //true is always 1 to 0. probability for the 511 Kev gamma to reach point in 1
  Float_t travel2_true; //true is always 0 to 1. probability for the 511 Kev gamma to reach point in 0
  Float_t travel2_false; //true is always 1 to 0. probability for the 511 Kev gamma to reach point in 1
  Float_t compton_true;
  Float_t compton_false;
  Float_t photo_true;
  Float_t photo_false;

  // Float_t travel1;
  std::vector<double> z0,z1,travel1,travel2,compton,total,photo,comptonPhoto,energy;


  int point = 0;
  for(int iDiv = 0 ; iDiv < divisions+1 ; iDiv++) //for all possible z0 and z1.
  {
    for(int jDiv = 0 ; jDiv < divisions+1 ; jDiv++)
    {
      zValue[0] = startZ0 + stepZ0*iDiv;
      zValue[1] = startZ1 + stepZ1*jDiv;
      z0.push_back(zValue[0]);
      z1.push_back(zValue[1]);
      Float_t point0[3];
      Float_t point1[3];
      Float_t source0[3];
      Float_t source1[3];
      point0[0] = 0;
      point0[1] = 0;
      point0[2] = zValue[0];
      point1[0] = 1.6;
      point1[1] = 0;
      point1[2] = zValue[1];
      source0[0] = 0;
      source0[1] = 0;
      source0[2] = -100;
      source1[0] = 1.6;
      source1[1] = 0;
      source1[2] = -100;


      travel1_true = transmissionProbability(zValue[0],lambdaLYSO->Eval(0.511));
      travel1_false = transmissionProbability(zValue[1],lambdaLYSO->Eval(0.511));
      travel1.push_back(travel1_true / travel1_false);
      Float_t eValue_true[2];
      Float_t eValue_false[2];
      eValue_true[0] = 0.511 - scatteredGammaEnergy(0.511,source0,point0,point1); // in crystal 0, energy deposited is the incident energy minus the energy of the outcoming scattered gamma
      eValue_true[1] = scatteredGammaEnergy(0.511,source0,point0,point1); // the scattered gamma will have a photoelectric effect in second crystal
      eValue_false[0] = scatteredGammaEnergy(0.511,source1,point1,point0); // in crystal 0, fot false hypo 1->0, energy deposited is the the scattered gamma energy
      eValue_false[1] = 0.511 - scatteredGammaEnergy(0.511,source1,point1,point0);

      compton_true  = simCompton(angle3D(source0,point0,point1));
      compton_false = simCompton(angle3D(source1,point1,point0));
      compton.push_back(compton_true / compton_false);

      photo_true = simPhotoelectric(photoelectricCrossSectionLYSO->Eval(eValue_true[1]));
      photo_false = simPhotoelectric(photoelectricCrossSectionLYSO->Eval(eValue_false[0]));
      photo.push_back(photo_true / photo_false);

      travel2_true  = transmissionProbability(distance3D(0,0,zValue[0],1.6,0,zValue[1]),lambdaLYSO->Eval(eValue_true[1]));
      travel2_false = transmissionProbability(distance3D(0,0,zValue[0],1.6,0,zValue[1]),lambdaLYSO->Eval(eValue_false[0]));
      travel2.push_back(travel2_true / travel2_false);

      total.push_back(
                     (travel1_true / travel1_false) *
                     (compton_true / compton_false) *
                     (photo_true / photo_false) *
                     (travel2_true / travel2_false)
                     );

      comptonPhoto.push_back(
                             (compton_true / compton_false) *
                             (photo_true / photo_false)
                             );
      // compatibility true and false energy values, given a flat energy resolution of 10% FWHM
      double res = 0.1;
      double prob0_true  = MeasureProb(eValue_true[0] ,eValue_true[0],(eValue_true[0] *res)/2.35);
      double prob1_true  = MeasureProb(eValue_true[1] ,eValue_true[1],(eValue_true[1] *res)/2.35);
      double prob0_false = MeasureProb(eValue_false[0],eValue_true[0],(eValue_false[0]*res)/2.35);
      double prob1_false = MeasureProb(eValue_false[1],eValue_true[1],(eValue_false[1]*res)/2.35);

      if(  (prob0_true*prob1_true) > (prob0_false*prob1_false)  ) energy.push_back(1);
      if(  (prob0_true*prob1_true) < (prob0_false*prob1_false)  ) energy.push_back(0);

      // energy.push_back( (prob0_true*prob1_true) / (prob0_false*prob1_false));

      std::cout << zValue[0] << "\t"
                << zValue[1] << "\t"
                // << distance3D(0,0,zValue[0],1.6,0,zValue[1]) << " "
                // << angle3D(source0,point0,point1) << " "
                // << angle3D(source1,point1,point0) << " "
                // << angle3D(source0,point0,point1) + angle3D(source1,point1,point0)
                << eValue_true[0] << "\t"
                << eValue_true[1] << "\t"
                // << eValue_true[0] + eValue_true[1] << " "
                << eValue_false[0] << "\t"
                << eValue_false[1] << "\t"
                // << eValue_false[0] + eValue_false[1] << " "
                // << photo_true << " "
                // << photo_false << " "
                // << prob0 << " "
                // << prob1 << " "
                << prob0_true << "\t"
                << prob1_true << "\t"
                << prob0_false << "\t"
                << prob1_false << "\t"
                << std::endl;

    }
  }


  std::string outFile = "simpleCompton.root";
  TFile* fOut = new TFile(outFile.c_str(),"recreate");
  comptonCrossSectionLYSO->Write();
  photoelectricCrossSectionLYSO->Write();
  lambdaLYSO->Write();


  TCanvas* c_total = new TCanvas("c_total", "c_total", 1200, 1200);
  TGraph2D* g_total = new TGraph2D(z0.size(),&z0[0],&z1[0],&total[0]);
  c_total->cd();
  g_total->SetName("Ratio Total");
  g_total->SetTitle("Ratio Total");
  g_total->GetXaxis()->SetTitle("z0 [mm]");
  g_total->GetYaxis()->SetTitle("z1 [mm]");
  g_total->GetZaxis()->SetTitle("total");
  g_total->Draw("surf1");
  c_total->Write();

  TCanvas* c_comptonPhoto = new TCanvas("c_comptonPhoto", "c_comptonPhoto", 1200, 1200);
  TGraph2D* g_comptonPhoto = new TGraph2D(z0.size(),&z0[0],&z1[0],&comptonPhoto[0]);
  c_comptonPhoto->cd();
  g_comptonPhoto->SetName("Compton*Photo");
  g_comptonPhoto->SetTitle("Compton*Photo");
  g_comptonPhoto->GetXaxis()->SetTitle("z0 [mm]");
  g_comptonPhoto->GetYaxis()->SetTitle("z1 [mm]");
  g_comptonPhoto->GetZaxis()->SetTitle("Compton*Photo");
  g_comptonPhoto->Draw("surf1");
  c_comptonPhoto->Write();

  TCanvas* c_energy = new TCanvas("c_energy", "c_energy", 1200, 1200);
  TGraph2D* g_energy = new TGraph2D(z0.size(),&z0[0],&z1[0],&energy[0]);
  c_energy->cd();
  g_energy->SetName("energy");
  g_energy->SetTitle("energy");
  g_energy->GetXaxis()->SetTitle("z0 [mm]");
  g_energy->GetYaxis()->SetTitle("z1 [mm]");
  g_energy->GetZaxis()->SetTitle("Energy prob");
  g_energy->Draw("surf1");
  c_energy->Write();

  TCanvas* c_travel1 = new TCanvas("c_travel1", "c_travel1", 1200, 1200);
  TGraph2D* g_travel1 = new TGraph2D(z0.size(),&z0[0],&z1[0],&travel1[0]);
  c_travel1->cd();
  g_travel1->SetName("Ratio Travel 1");
  g_travel1->SetTitle("Ratio Travel 1");
  g_travel1->GetXaxis()->SetTitle("z0 [mm]");
  g_travel1->GetYaxis()->SetTitle("z1 [mm]");
  g_travel1->GetZaxis()->SetTitle("travel1");
  g_travel1->Draw("surf1");
  c_travel1->Write();

  TCanvas* c_travel2 = new TCanvas("c_travel2", "c_travel2", 1200, 1200);
  TGraph2D* g_travel2 = new TGraph2D(z0.size(),&z0[0],&z1[0],&travel2[0]);
  c_travel2->cd();
  g_travel2->SetName("Ratio Travel 2");
  g_travel2->SetTitle("Ratio Travel 2");
  g_travel2->GetXaxis()->SetTitle("z0 [mm]");
  g_travel2->GetYaxis()->SetTitle("z1 [mm]");
  g_travel2->GetZaxis()->SetTitle("travel2");
  g_travel2->Draw("surf1");
  c_travel2->Write();

  TCanvas* c_compton = new TCanvas("c_compton", "c_compton", 1200, 1200);
  TGraph2D* g_compton = new TGraph2D(z0.size(),&z0[0],&z1[0],&compton[0]);
  c_compton->cd();
  g_compton->SetName("Ratio Compton");
  g_compton->SetTitle("Ratio Compton");
  g_compton->GetXaxis()->SetTitle("z0 [mm]");
  g_compton->GetYaxis()->SetTitle("z1 [mm]");
  g_compton->GetZaxis()->SetTitle("compton");
  g_compton->Draw("surf1");
  c_compton->Write();

  TCanvas* c_photo = new TCanvas("c_photo", "c_photo", 1200, 1200);
  TGraph2D* g_photo = new TGraph2D(z0.size(),&z0[0],&z1[0],&photo[0]);
  c_photo->cd();
  g_photo->SetName("Ratio Photoelectric");
  g_photo->SetTitle("Ratio Photoelectric");
  g_photo->GetXaxis()->SetTitle("z0 [mm]");
  g_photo->GetYaxis()->SetTitle("z1 [mm]");
  g_photo->GetZaxis()->SetTitle("photo");
  g_photo->Draw("surf1");
  c_photo->Write();



  fOut->Close();

  return 0;
}
