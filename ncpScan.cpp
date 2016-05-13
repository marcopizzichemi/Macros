// g++ -o ../build/ncpScan ncpScan.cpp `root-config --cflags --glibs` 
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
#include <TH3I.h>
#include <TMath.h>
#include <TCut.h>
#include <TCutG.h>
#include <TH1D.h>
#include <TTree.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TLegend.h>

#include <sys/types.h>
#include <sys/stat.h>
#include <unistd.h>


Double_t thetaFunction(Double_t *x, Double_t *par)
{
    Double_t f;
    Float_t xx =x[0];
    if (xx>par[0] && xx<par[1]) 
    {
        f = par[2];    
    }
    else
    {
        f = 0;    
    }
    return f;
}

Double_t sigmoid(Double_t *x, Double_t *par)
{
    Double_t f;
    Float_t xx = x[0];
    f =  par[0] + (par[1] / ( 1 + TMath::Exp(-par[2]*(xx-par[3])))); 
    return f;
}


int main(int argc, char * argv[])
{
  


  TFile *f = new TFile(argv[1]);
  TFile *f2 = new TFile(argv[2]);
  
  TTree *tree =  (TTree*) f->Get("adc");
  
//   TCut CutXYZ = TCut;
//   TCut CutTrigger;
  TCutG *CutZX = (TCutG*) f2->Get("cutg_0_3_B2");
  TCutG *CutZY = (TCutG*) f2->Get("cutg_1_3_B2");
  
  
  TCut broadCut = "ch6 > 1100";
  TCut crystal = "FloodX > -1.5 && FloodX < -1 && FloodY > -1.45 &&FloodY < -0.9";
//   TCut crystal = "FloodX > -4 && FloodX < -3.6 && FloodY > -1.28 &&FloodY < -1.02";
  TCut photopeak = "ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 > 9500 && ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 < 12500";
//   TCut photopeak = "ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 > 8500 && ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 < 11000";
  TCut trigger = "TriggerChannel == 6";
//   TCut trigger = "TriggerChannel == 2";
  std::string w = "ch6/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15)";
//   std::string w = "ch2/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15)";  
  
  TCut wCut = "FloodZ > 0.3 && FloodZ < 0.5";
  
  TH1F *h = new TH1F("h","h",1000,0,35000);
  TH1F *hSum = new TH1F("hSum","hSum",1000,0,35000);
  TH1F *hCorrSig = new TH1F("hCorrSig","hCorrSig",1000,0,35000);
  TH1F *hCorrLinear = new TH1F("hCorrLinear","hCorrLinear",1000,0,35000);
  TH2F *h2 =new TH2F ("h2","h2",1000,-7,7,1000,-7,7);
  TH3I *h3 = new TH3I("h3","h3",100,-3,0,100,-3,0,100,0,1);
  TH1F *hAll = new TH1F("hAll","hAll",500,0,1);
  TH1F *hNear = new TH1F("hNear","hNear",500,0,1);
  TH1F *hFloodZ = new TH1F("hFloodZ","hFloodZ",500,0,1);
  TH2F *scatter = new TH2F("scatter","scatter",500,0,1,500,0,1);
  TH2F *correction = new TH2F("correction","correction",500,0,1,1000,0,35000);
  
  
  
  
  TCanvas *multi = new TCanvas("multi","multi",1800,1200);
  multi->Divide(3,3);
  multi->cd(1);
  
//   h2->GetXaxis()->SetRangeUser(-4.1,-3.5);
//   h2->GetYaxis()->SetRangeUser(-1.3,-1.0);
//   std::cout << "quiiiiiiiiiiiiiiiiii" << std::endl;
  tree->Draw("FloodZ:FloodY:FloodX >> h3",broadCut+trigger+CutZX->GetName()+CutZY->GetName());
//   std::cout << "dopo" << std::endl;
  h2->GetXaxis()->SetRangeUser(-1.5,-1);
  h2->GetYaxis()->SetRangeUser(-1.45,-0.9);

  h3->GetXaxis()->SetRangeUser(-3,0);
  h3->GetYaxis()->SetRangeUser(-3,0);

//   tree->Draw("FloodY:FloodX >> h2","","COLZ");
//   multi->cd(2);
//   tree->Draw("ch1+ch2+ch3+ch5+ch6+ch7+ch9+ch10+ch11 >> h",crystal);
  
  multi->cd(2);
  tree->Draw("ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 >> hSum",broadCut+trigger+CutZX->GetName()+CutZY->GetName());
  
  multi->cd(3);
//   tree->Draw("ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 : ch2/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) >> correction",crystal+photopeak,"COLZ");
  tree->Draw("ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 : FloodZ >> correction",trigger+photopeak+broadCut+CutZX->GetName()+CutZY->GetName(),"COLZ");
  
  multi->cd(4);
  correction->FitSlicesY(0, 0, -1, 0, "RQ");
  TH1D *spectrum2d_1 = (TH1D*)gDirectory->Get("correction_1"); // _1 is the TH1D automatically created by ROOT when FitSlicesX is called, holding the TH1F of the mean values
  spectrum2d_1->GetYaxis()->SetRangeUser(10200,11400);
  TF1 *sig1 = new TF1("sig1", sigmoid,0,1,4);
  sig1->SetParameter( 0, 10500); 
  sig1->SetParameter( 1, 700); 
  sig1->SetParameter( 2, 30);
  sig1->SetParameter( 3, 0.4);
  
  double a = sig1->GetParameter(0);
  double b = sig1->GetParameter(1);
  double c = sig1->GetParameter(2);
  double d = sig1->GetParameter(3);
  
  spectrum2d_1->Fit(sig1,"QR");
  double medianPoint = d - (1/c)*TMath::Log(3);
  medianPoint = 0.3966;
  std::cout << "medianPoint = " << medianPoint << std::endl;
  multi->cd(5);
  
  std::stringstream baseVar,var;
  std::stringstream sSig0,sSigw;
  sSig0 << "(" 
      << a
      << " + (" 
      << b  
      << " * (  0.5 - 1/(1 + TMath::Exp( -"
      << c
      << " * ( "
      << medianPoint
      << " - "
      << d
      << ") ) ) ) ) ) ";
  
  sSigw << "("
      << a
      << " + (" 
      << b  
      << " * (  0.5 - 1/(1 + TMath::Exp( -"
      << c
      << " * ( "
      << "FloodZ"
      << " - "
      << d
      << ") ) ) ) ) )";
  
  baseVar << "((ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) * (" << sSig0.str() << " / " << sSigw.str() << "))";
//   baseVar << "((ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) + ( " 
//       << b  
//       << " * (  0.5 - 1/(1 + TMath::Exp( -"
//       << c
//       << " * ( "
//       << "FloodZ"
//       << " - "
//       << d
//       << ") ) ) ) )) ";
  
  var << baseVar.str()  <<  ">> hCorrSig";
  tree->Draw(var.str().c_str(),broadCut+trigger+CutZX->GetName()+CutZY->GetName());
  TF1 *gauss = new TF1("gauss","[0]*exp(-0.5*((x-[1])/[2])**2)",10000,12000);
  gauss->SetParameter(0,1600);
  gauss->SetParameter(1,10000);
  gauss->SetParameter(2,600);
  hCorrSig->Fit(gauss,"RQ");
  
  std::stringstream cutPeakSig ;
  cutPeakSig << baseVar.str() << ">" <<  gauss->GetParameter(1) - 4.0*gauss->GetParameter(2) << " && " <<  baseVar.str() << "<" <<  gauss->GetParameter(1) + 4.0*gauss->GetParameter(2) ;
  TCut photopeakCorrSig = cutPeakSig.str().c_str();

  TH1F *clone = new TH1F("clone","clone",1000,0,35000);
  
  clone->SetFillStyle(3001);
  clone->SetFillColor(3);
  var.str("");
  var << baseVar.str()  <<  ">> clone";
  tree->Draw(var.str().c_str(),broadCut+trigger+CutZX->GetName()+CutZY->GetName()+photopeakCorrSig,"same");
//   clone->Draw("same");
  var.str("");
  
  multi->cd(7);
  
  TH1D* spectrum2d_1_copy =  (TH1D*) spectrum2d_1->Clone();
  spectrum2d_1_copy->Draw();
  spectrum2d_1_copy->GetYaxis()->SetRangeUser(10200,11400);
  TF1 *linearCrystal = new TF1("linearCrystal",  "[0]*x + [1]",0.35,0.42);
  spectrum2d_1_copy->Fit(linearCrystal,"RQ");
  
  multi->cd(8);
  baseVar.str("");
  baseVar << "(( ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 ) * ( (" << linearCrystal->GetParameter(0) * medianPoint + linearCrystal->GetParameter(1) << ") / (" << linearCrystal->GetParameter(0) << "* FloodZ + "<< linearCrystal->GetParameter(1) << " ))) >> hCorrLinear";
//   baseVar << "((   ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15   ) - ( (" << w << " - " <<  d - (1/c)*TMath::Log(1 + (b/(a+(b/2) -a))) << " ) * ( " << linearCrystal->GetParameter(0) << ") )) >> hCorrLinear";
  tree->Draw(baseVar.str().c_str(),broadCut+trigger+CutZX->GetName()+CutZY->GetName());
  
//   TF1 *gauss2 = new TF1("gauss2","[0]*exp(-0.5*((x-[1])/[2])**2)",10000,12000);
//   gauss2->SetParameter(0,1600);
//   gauss2->SetParameter(1,10000);
//   gauss2->SetParameter(2,600);
//   hCorrLinear->Fit(gauss2,"RQ");
//   multi->cd(7);
//   tree->Draw("FloodZ >> hFloodZ",crystal+photopeak);
//   TF1 *fit1 = new TF1("fit1",thetaFunction,0,1,3);
//   fit1->SetParameter( 0, 0.36); // on this w histo, the first bin above 20% max
//   fit1->SetParameter( 1, 0.54);  // on this w histo, the last bin above 20% max
//   fit1->SetParameter( 2, 1200);
//   hFloodZ->Fit(fit1,"R");
  multi->cd(6);
  tree->Draw("ch6/(ch1+ch2+ch3+ch5+ch6+ch7+ch9+ch10+ch11) >> hNear",broadCut+ photopeakCorrSig+CutZX->GetName()+CutZY->GetName());
  TF1 *fit2 = new TF1("fit2",thetaFunction,0,1,3);
  fit2->SetParameter( 0, 0.36); // on this w histo, the first bin above 20% max
  fit2->SetParameter( 1, 0.54);  // on this w histo, the last bin above 20% max
  fit2->SetParameter( 2, 1200);
  hNear->Fit(fit2,"RQ");
  multi->cd(9);
  tree->Draw("ch6/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15) >> hAll",broadCut+photopeakCorrSig+CutZX->GetName()+CutZY->GetName());
  TF1 *fit3 = new TF1("fit1",thetaFunction,0,1,3);
  fit3->SetParameter( 0, 0.1); // on this w histo, the first bin above 20% max
  fit3->SetParameter( 1, 0.8);  // on this w histo, the last bin above 20% max
  fit3->SetParameter( 2, 1200);
  hAll->Fit(fit3,"RQ");
  
  std::cout <<  linearCrystal->GetParameter(0) << " " << linearCrystal->GetParameter(1) << std::endl;
  std::cout <<  fit2->GetParameter(0) << " " << fit2->GetParameter(1) << std::endl;
  std::cout <<  fit3->GetParameter(0) << " " << fit3->GetParameter(1) << std::endl;
  std::cout << "Energy Resolution Corrected (Sigmoid) = " << gauss->GetParameter(2)*2.355 /  gauss->GetParameter(1) << std::endl;
  std::cout << "Near Channels - Delta w = " << fit2->GetParameter(1) - fit2->GetParameter(0) << std::endl;
  std::cout << "All  Channels - Delta w = " << fit3->GetParameter(1) - fit3->GetParameter(0) << std::endl;
  
  TFile *fOut= new TFile("afile.root","RECREATE");
  fOut->cd();
  multi->Write();
  fOut->Close();
  
  return 0;
  
}
















