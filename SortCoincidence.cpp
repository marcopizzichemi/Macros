// compile with 
// g++ -Wl,--no-as-needed -o SortCoincidence SortCoincidence.cpp `root-config --cflags --glibs` -lSpectrum -lMLP -lTreePlayer

#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3F.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TFile.h"
#include "TF2.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TMinuit.h"
#include "TGraph2D.h"

#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <vector>
#include <algorithm>

#include <stdlib.h> 
#include <stdio.h> 
#include <unistd.h>
#include <cmath> 

int main(int argc, char** argv)
{
  gStyle->SetOptStat(0);
  TString filename = argv[1];
  TFile * file = new TFile(filename, "READ");
  TTree *tree = (TTree *)file->Get("adc");
  
  Float_t floodx[2],floody[2],ratio[2];
  int trigger[2];
  bool bad;
  Short_t charge[32];
  
  tree->SetBranchAddress("TriggerChannel_0",&trigger[0]);
  tree->SetBranchAddress("FloodX_0",&floodx[0]);
  tree->SetBranchAddress("FloodY_0",&floody[0]);
  tree->SetBranchAddress("RatioTriggerOnAll_0",&ratio[0]);
  tree->SetBranchAddress("TriggerChannel_1",&trigger[1]);
  tree->SetBranchAddress("FloodX_1",&floodx[1]);
  tree->SetBranchAddress("FloodY_1",&floody[1]);
  tree->SetBranchAddress("RatioTriggerOnAll_1",&ratio[1]);
  tree->SetBranchAddress("BadEvent",&bad);
  
  std::stringstream sname,stype;
  std::string name,type;
   
  for(int i=0; i<32; i++)
  {
    sname.str(std::string());
    stype.str(std::string());
    sname << "ch" << i;
    name = sname.str();
    tree->SetBranchAddress(name.c_str(), &charge[i]);
  }
  
  //TCuts
  TCut Crystals[32]=
  {
    "TriggerChannel_0 == 0"                                                                                                                                                                                                                                                                                                                                                                ,
    "TriggerChannel_0 == 1 && (( TMath::Power(FloodX_0 - -3.668,2) / TMath::Power(0.412041,2))  + ( TMath::Power(FloodY_0 - 1.036,2) / TMath::Power(0.412041,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3786.52 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4915.66"                       ,
    "TriggerChannel_0 == 2 && (( TMath::Power(FloodX_0 - -3.724,2) / TMath::Power(0.385788,2))  + ( TMath::Power(FloodY_0 - -1.708,2) / TMath::Power(0.385788,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3867.03 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4956.38"                      ,
    "TriggerChannel_0 == 3 && (( TMath::Power(FloodX_0 - -3.612,2) / TMath::Power(0.438224,2))  + ( TMath::Power(FloodY_0 - -3.892,2) / TMath::Power(0.438224,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3086.97 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4119.18"                      ,
    "TriggerChannel_0 == 4 && (( TMath::Power(FloodX_0 - -0.98,2) / TMath::Power(0.395637,2))  + ( TMath::Power(FloodY_0 - 3.892,2) / TMath::Power(0.395637,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3629.98 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4637.16"                        ,
    "TriggerChannel_0 == 5 && (( TMath::Power(FloodX_0 - -1.148,2) / TMath::Power(0.326781,2))  + ( TMath::Power(FloodY_0 - 1.092,2) / TMath::Power(0.326781,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 4321.57 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 5574.69"                       ,
    "TriggerChannel_0 == 6 && (( TMath::Power(FloodX_0 - -1.204,2) / TMath::Power(0.318083,2))  + ( TMath::Power(FloodY_0 - -1.54,2) / TMath::Power(0.318083,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 4228.77 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 5387.56"                       ,
    "TriggerChannel_0 == 7 && (( TMath::Power(FloodX_0 - -1.316,2) / TMath::Power(0.38977,2))  + ( TMath::Power(FloodY_0 - -4.228,2) / TMath::Power(0.38977,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3736.52 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4746.55"                        ,
    "TriggerChannel_0 == 8 && (( TMath::Power(FloodX_0 - 1.596,2) / TMath::Power(0.379766,2))  + ( TMath::Power(FloodY_0 - 3.948,2) / TMath::Power(0.379766,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3797.15 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4955.01"                        ,
    "TriggerChannel_0 == 9 && (( TMath::Power(FloodX_0 - 1.428,2) / TMath::Power(0.339249,2))  + ( TMath::Power(FloodY_0 - 1.148,2) / TMath::Power(0.339249,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3851.75 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4933.12"                        ,
    "TriggerChannel_0 == 10 && (( TMath::Power(FloodX_0 - 1.26,2) / TMath::Power(0.31409,2))  + ( TMath::Power(FloodY_0 - -1.54,2) / TMath::Power(0.31409,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 4343.24 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 5622.86"                          ,
    "TriggerChannel_0 == 11 && (( TMath::Power(FloodX_0 - 1.54,2) / TMath::Power(0.439377,2))  + ( TMath::Power(FloodY_0 - -4.172,2) / TMath::Power(0.439377,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3041.24 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 3901.02"                       ,
    "TriggerChannel_0 == 12 && (( TMath::Power(FloodX_0 - 4.004,2) / TMath::Power(0.327012,2))  + ( TMath::Power(FloodY_0 - 4.172,2) / TMath::Power(0.327012,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3441.52 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4418.59"                       ,
    "TriggerChannel_0 == 13 && (( TMath::Power(FloodX_0 - 3.724,2) / TMath::Power(0.389074,2))  + ( TMath::Power(FloodY_0 - 1.036,2) / TMath::Power(0.389074,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3629.07 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4667.87"                       ,
    "TriggerChannel_0 == 14 && (( TMath::Power(FloodX_0 - 3.78,2) / TMath::Power(0.368637,2))  + ( TMath::Power(FloodY_0 - -1.484,2) / TMath::Power(0.368637,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 4076.73 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 5264.54"                       ,
    "TriggerChannel_0 == 15 && (( TMath::Power(FloodX_0 - 3.556,2) / TMath::Power(0.468446,2))  + ( TMath::Power(FloodY_0 - -3.948,2) / TMath::Power(0.468446,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 2822.28 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 3681.71"                      ,
    "TriggerChannel_1 == 16 && (( TMath::Power(FloodX_1 - -3.724,2) / TMath::Power(0.443754,2))  + ( TMath::Power(FloodY_1 - 3.892,2) / TMath::Power(0.443754,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3132.72 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4266.53"  ,
    "TriggerChannel_1 == 17 && (( TMath::Power(FloodX_1 - -3.724,2) / TMath::Power(0.439618,2))  + ( TMath::Power(FloodY_1 - 1.316,2) / TMath::Power(0.439618,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3002.34 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4135.14"  ,
    "TriggerChannel_1 == 18 && (( TMath::Power(FloodX_1 - -3.668,2) / TMath::Power(0.406535,2))  + ( TMath::Power(FloodY_1 - -1.26,2) / TMath::Power(0.406535,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3366.09 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4747.98"  ,
    "TriggerChannel_1 == 19 && (( TMath::Power(FloodX_1 - -3.836,2) / TMath::Power(0.454168,2))  + ( TMath::Power(FloodY_1 - -3.948,2) / TMath::Power(0.454168,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3169.97 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4260.33" ,
    "TriggerChannel_1 == 20 && (( TMath::Power(FloodX_1 - -1.26,2) / TMath::Power(0.435707,2))  + ( TMath::Power(FloodY_1 - 4.06,2) / TMath::Power(0.435707,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3115.91 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4299.74"    ,
    "TriggerChannel_1 == 21 && (( TMath::Power(FloodX_1 - -1.54,2) / TMath::Power(0.375443,2))  + ( TMath::Power(FloodY_1 - 1.428,2) / TMath::Power(0.375443,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3354 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4641.48"      ,
    "TriggerChannel_1 == 22 && (( TMath::Power(FloodX_1 - -1.372,2) / TMath::Power(0.368213,2))  + ( TMath::Power(FloodY_1 - -1.26,2) / TMath::Power(0.368213,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3160.73 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4253.92"  ,
    "TriggerChannel_1 == 23 && (( TMath::Power(FloodX_1 - -1.428,2) / TMath::Power(0.457071,2))  + ( TMath::Power(FloodY_1 - -4.284,2) / TMath::Power(0.457071,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3125.02 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4301.74" ,
    "TriggerChannel_1 == 24 && (( TMath::Power(FloodX_1 - 1.316,2) / TMath::Power(0.437403,2))  + ( TMath::Power(FloodY_1 - 4.172,2) / TMath::Power(0.437403,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 2662.73 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 3593.6"    ,
    "TriggerChannel_1 == 25"                                                                                                                                                                                                                                                                                                                                                                ,
    "TriggerChannel_1 == 26 && (( TMath::Power(FloodX_1 - 1.204,2) / TMath::Power(0.379761,2))  + ( TMath::Power(FloodY_1 - -1.428,2) / TMath::Power(0.379761,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > -510629 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 48433.2"  ,
    "TriggerChannel_1 == 27 && (( TMath::Power(FloodX_1 - 1.372,2) / TMath::Power(0.447345,2))  + ( TMath::Power(FloodY_1 - -4.284,2) / TMath::Power(0.447345,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3241.13 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4519.65"  ,
    "TriggerChannel_1 == 28 && (( TMath::Power(FloodX_1 - 4.228,2) / TMath::Power(0.197975,2))  + ( TMath::Power(FloodY_1 - 4.172,2) / TMath::Power(0.197975,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 2950.12 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 3987.89"   ,
    "TriggerChannel_1 == 29 && (( TMath::Power(FloodX_1 - 3.836,2) / TMath::Power(0.444737,2))  + ( TMath::Power(FloodY_1 - 1.316,2) / TMath::Power(0.444737,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3334.08 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4473.93"   ,
    "TriggerChannel_1 == 30 && (( TMath::Power(FloodX_1 - 3.668,2) / TMath::Power(0.435218,2))  + ( TMath::Power(FloodY_1 - -1.316,2) / TMath::Power(0.435218,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 3244.59 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 4411.88"  ,
    "TriggerChannel_1 == 31 && (( TMath::Power(FloodX_1 - 3.668,2) / TMath::Power(0.478436,2))  + ( TMath::Power(FloodY_1 - -4.172,2) / TMath::Power(0.478436,2))) < 1 && (ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) > 2831.63 &&(ch16+ch17+ch18+ch19+ch20+ch21+ch22+ch23+ch24+ch25+ch26+ch27+ch28+ch29+ch30+ch31) < 3839.56"
  };
  
  
  
  double xPos[32] = 
  {
    0,
    -3.668,
-3.724,
-3.612,
-0.98 ,
-1.148,
-1.204,
-1.316,
1.596 ,
1.428 ,
1.26  ,
3.836, 
4.004, 
3.724, 
3.78, 
3.556 ,
     -3.724,
 -3.724,
 -3.668,
 -3.836,
 -1.26, 
 -1.54, 
 -1.372,
 -1.428,
 1.316,
 0,
 1.204, 
 1.372, 
 3.892, 
 3.836, 
 3.668, 
 3.668
  };
  double yPos[32] =
  {
    0,
    1.036, 
-1.708,
-3.892,
3.892, 
1.092, 
-1.54, 
-4.228,
3.948, 
1.148, 
-1.54, 
3.724, 
4.172, 
1.036, 
-1.484,
-3.948,
 3.892 ,
 1.316 ,
 -1.26 ,
 -3.948,
 4.06  ,
 1.428 ,
 -1.26 ,
 -4.284,
 4.172 ,
 0,
 -1.428,
 -4.284,
 4.06  ,
 1.316 ,
 -1.316,
 -4.172
  };
  
  double xSigma[32]
  {
    0,
    0.412041,
0.385788,
0.438224,
0.395637,
0.326781,
0.318083,
0.38977 ,
0.379766,
0.339249,
0.31409 ,
0.151452,
0.327012,
0.389074,
0.368637,
0.468446,
0.443754,
0.439618,
0.406535,
0.454168,
0.435707,
0.375443,
0.368213,
0.457071,
0.437403,
0,
0.379761,
0.447345,
0.1562,  
0.444737,
0.435218,
0.478436
  };
  
  double ySigma[32]
  {
   0,
    0.412041,
0.385788,
0.438224,
0.395637,
0.326781,
0.318083,
0.38977 ,
0.379766,
0.339249,
0.31409 ,
0.151452,
0.327012,
0.389074,
0.368637,
0.468446,
 0.443754,
0.439618,
0.406535,
0.454168,
0.435707,
0.375443,
0.368213,
0.457071,
0.437403,
0.379761,
0.447345,
0,
0.1562,  
0.444737,
0.435218,
0.478436
  };
  
  
  double eMin[32] =
  {
    0         ,
    3786.52   ,
3867.03       ,
3086.97       ,
3629.98       ,
4321.57       ,
4228.77       ,
3736.52       ,
3797.15       ,
3851.75       ,
4343.24       ,
3663.34       ,
3441.52       ,
3629.07       ,
4076.73       ,
2822.28       ,
3132.72       ,
3002.34       ,
3366.09       ,
3169.97       ,
3115.91       ,
3354          ,
3160.73       ,
3125.02       ,
2662.73       ,
-510629       ,
3241.13       ,
3128.64       ,
3334.08       ,
3244.59       ,
2831.63
    
  };
  
  double eMax[32] =
  {
    0        ,
    4915.66  ,
4956.38      ,
4119.18      ,
4637.16      ,
5574.69      ,
5387.56      ,
4746.55      ,
4955.01      ,
4933.12      ,
5622.86      ,
4609.81      ,
4418.59      ,
4667.87      ,
5264.54      ,
3681.71      ,
4266.53      ,
4135.14      ,
4747.98      ,
4260.33      ,
4299.74      ,
4641.48      ,
4253.92      ,
4301.74      ,
3593.6       ,
48433.2      ,
4519.65      ,
4021.49      ,
4473.93      ,
4411.88      ,
3839.56
  };
  
  
  //increase the limits in TFormula, because the cut string is way too long for the defaults 
  //TFormula::SetMaxima(100000,1000,1000000);
  //create the formula from the cuts
//   TTreeFormula* Formula = new TTreeFormula("Formula",MegaString.c_str(),t1);
  
  //mapping for the new boards
  Double_t xmppc[32]={-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8,-4.8,-4.8,-4.8,-4.8,-1.6,-1.6,-1.6,-1.6,1.6,1.6,1.6,1.6,4.8,4.8,4.8,4.8};
  Double_t ymppc[32]={4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8};
  
  
  TH2F *histo2= new TH2F("ciao2","",4,-6.2,6.2,4,-6.2,6.2);
  Long64_t counter = 0;
  //loop on entries
  Long64_t nentries = tree->GetEntriesFast();   //nentries num eventi
  for (Long64_t jentry=0; jentry<nentries;jentry++)
  {
    tree->GetEntry(jentry);
    
    bool chBoard0 = false;
    bool chBoard1 = false;
    float z0= -3.9;    
    float z1= 3.9;
    float w = z0 - z1;
    float s = (0 - z0)/w;
    float x0,x1,y0,y1;
    
    //calculate sum charge
    
    float sumCharge0 = 0;
    float sumCharge1 = 0;
    
    for(int i = 0 ; i < 16 ; i++)
    {
      sumCharge0 += charge[i];
    }
    for(int i = 16 ; i < 32 ; i++)
    {
      sumCharge1 += charge[i];
    }
    
    
    for(int i = 0 ; i < 32 ; i++)
    {
      //TriggerChannel_0 == 1 &&  (( TMath::Power(FloodX_0 - -3.668,2) / TMath::Power(0.412041,2))  + ( TMath::Power(FloodY_0 - 1.036 ,2 ) / TMath::Power(0.412041,2))) < 1 && (ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) > 3786.52 &&(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) < 4915.66
      if(i < 16) //on board0
      {
	if( (((TMath::Power(floodx[0] - xPos[i],2)/ TMath::Power(xSigma[i],2))  + ( TMath::Power(floody[0] - yPos[i] ,2 ) / TMath::Power(ySigma[i],2))) < 1 )    &&    sumCharge0 > eMin[i] && sumCharge0 < eMax[i]              )
	{
	  if(chBoard0 == true) // 2 511 in the same matrix, get rid of the event
	  {
	    continue;
	  }
	  else
	  {
	    chBoard0 = true;
	    x0 = xmppc[i];
	    y0 = ymppc[i];
	  }
	}
      }
      else //it was on board1
      {
	if( (((TMath::Power(floodx[1] - xPos[i],2)/ TMath::Power(xSigma[i],2))  + ( TMath::Power(floody[1] - yPos[i] ,2 ) / TMath::Power(ySigma[i],2))) < 1 )    &&    sumCharge1 > eMin[i] && sumCharge1 < eMax[i]              )
	{
	  if(chBoard1 == true) // 2 511 in the same matrix, get rid of the event
	  {    
	    continue; 
	  }
	  else    
	  {
	    chBoard1 = true;
	    x1 = xmppc[i];
	    y1 = ymppc[i];
	  }
	}
      }
    }
    
//     for(int i = 0 ; i < 32 ; i++)
//     {
//       TTreeFormula* Formula = new TTreeFormula("Formula",Crystals[i],tree);
//       if( Formula->EvalInstance() )
//       {
// 	if(i < 16) //it was on board0
// 	{
// 	  if(chBoard0 == true) // 2 511 in the same matrix, get rid of the event
// 	  {
// 	    continue;
// 	  }
// 	  else
// 	  {
// 	    chBoard0 = true;
// 	    x0 = xmppc[i];
// 	    y0 = ymppc[i];
// 	  }
// 	}
// 	else //it was on board1
// 	{
// 	  if(chBoard1 == true) // 2 511 in the same matrix, get rid of the event
// 	  {
// 	    continue;
// 	  }
// 	  else
// 	  {
// 	    chBoard1 = true;
// 	    x1 = xmppc[i];
// 	    y1 = ymppc[i];
// 	  }
// 	}
//       }
//     }
    
    if(chBoard0&&chBoard1)
    {
      //std::cout << " coincidence " <<  x0 << " " << y0 << " " << z0 << " " << x1 << " " << y1 << " " << z1 << std::endl; 
      float u = x0-x1;
      float v = y0-y1;
      float mx = x0 + s*u;
      float my = y0 + s*v;
      histo2->Fill(mx,my);
    }
    //std::cout << std::endl;
    /*
    if( Formula->EvalInstance() )
    {
      TimeTag = ExtendedTimeTag;
      CrystalID = 0; //TODO
      t2->Fill();
    }*/
    
    counter++;
    
    int perc = ((100*counter)/nentries); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
    
  }
  
  TCanvas *canvas= new TCanvas("ciao1", "Distribuzione di Poisson");
  histo2->Draw("COLZ");
  TFile *outFile = new TFile ("histo2d.root","recreate");
  //histo->Write();
  histo2->Write();
  canvas->Write();
  outFile->Close();
  
  
  
  return 0;
}