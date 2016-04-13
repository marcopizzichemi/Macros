// g++ -o modulePlot modulePlot.cpp `root-config --cflags --glibs` 
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
#include <string>
#include <iomanip>      // std::setprecision
#include <TCanvas.h>
#include <TH2F.h>
#include <TF1.h>
#include <TFile.h>
#include <TH1F.h>
#include <TGraphErrors.h>
#include <THStack.h>
#include <TLegend.h>


int main(int argc, char * argv[])
{
  TFile *f = new TFile(argv[1]);
  f->cd("Module 0.0");
  
  //-----------------------------//
  //      En Res Corrected
  //-----------------------------//
  c = (TCanvas*) gDirectory->Get("Corrected Energy res FWHM vs. i,j");
  
  spectrum2d = new TH2F();
  spectrum2d = (TH2F*)c->GetPrimitive("Corrected Energy res FWHM vs. i,j");
  all = new TH1F("Corrected Energy res FWHM","",100,0,1);
  all->GetXaxis()->SetTitle("Corrected Energy Resolution FWHM");
  all->SetName("Corrected Energy Resolution FWHM");
//   all->GetYaxis()->SetTitle("N");
  all->SetStats(0);
  central = new TH1F("central","",100,0,1);
//   all->GetYaxis()->SetTitle("N");
  central->SetStats(0);
  central->SetFillStyle(3001);
  central->SetFillColor(kRed);
  double sumEN_corr = 0;
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
  double stdEN_corr = 0;
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
  cenergyResolution2DCorrected->Print("energyResolution2DCorrected.png");
  cc->Print("energyResolutionCorrected.png");
  delete spectrum2d;
  delete all;
  delete central;
  delete hs;
  delete cc;
  delete legend;
  
  
  std::cout << "En. Res Central Channels [FWHM] = " << sumEN/16 << " +/- " << stdEN << std::endl;
  std::cout << "En. Res Corrected Central Channels [FWHM] = " << sumEN_corr/16 << " +/- " << stdEN_corr << std::endl;
  std::cout << "Light Output Central Channels [ADC Ch.] = " << sumLO/16 << " +/- " << stdLO << std::endl;
  
  
  
  return 0;
//   TFile *f = new TFile("farSourceGlassTop.root");
//   f->cd("Module 0.0");
// 
//   //-----------------------------//
//   //         En Res
//   //-----------------------------//
//   //   TCanvas *c = (TCanvas*) gDirectory->Get("Flood Histogram 2D");
//   TCanvas *c = (TCanvas*) gDirectory->Get("Energy res FWHM vs. i,j");
//   
//   TH2F *spectrum2d = new TH2F();
//   spectrum2d = (TH2F*)c->GetPrimitive("Energy res FWHM vs. i,j");
//   TH1F *all = new TH1F("Energy res FWHM","",100,0,1);
//   all->GetXaxis()->SetTitle("Energy Resolution FWHM");
//   all->SetName("Energy Resolution FWHM");
//   //   all->GetYaxis()->SetTitle("N");
//   all->SetStats(0);
//   TH1F *central = new TH1F("central","",100,0,1);
// //   all->GetYaxis()->SetTitle("N");
//   central->SetStats(0);
//   central->SetFillStyle(3001);
//   central->SetFillColor(kRed);
//   double sumEN = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	central->Fill(spectrum2d->GetBinContent(i+1,j+1));
// 	sumEN += spectrum2d->GetBinContent(i+1,j+1);
//       }
//       all->Fill(spectrum2d->GetBinContent(i+1,j+1));
//     } 
//   }
//   THStack *hs = new THStack("hs","");
//   double averageEN = sumEN/16 ;
//   double stdEN = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	stdEN += (1.0/16.0)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageEN,2));
//       }
//     }
//   }
//   TCanvas *cc = new TCanvas("cc","",1200,800);
//   all->SetFillStyle(3001);
//   all->SetFillColor(kBlue);
//   cc->cd();
//   hs->Add(all);
//   hs->Add(central);
// //   all->Draw();
// //   central->Draw("same");
//   TLegend *legend = new TLegend(0.5,0.62,0.893,0.89,"");
//   legend->SetFillStyle(0);
//   legend->AddEntry(all,"All channels","f");
//   legend->AddEntry(central,"Central channels","f");
//   hs->Draw("nostack");
//   hs->GetXaxis()->SetTitle("Energy Resolution FWHM");
//   hs->GetXaxis()->SetTitleOffset(1);
//   hs->GetXaxis()->SetTitleSize(0.045);
//   hs->GetXaxis()->SetLabelSize(0.045);
//   hs->GetYaxis()->SetLabelSize(0.045);
//   legend->Draw();
//   
//   spectrum2d->GetYaxis()->SetTitleOffset(1.8);
//   spectrum2d->GetXaxis()->SetTitleOffset(1.5);
//   spectrum2d->GetZaxis()->SetTitleOffset(1.7);
//   spectrum2d->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d->GetZaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d->GetZaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
//   spectrum2d->GetYaxis()->SetNdivisions(8,2,0, kTRUE);
//   
//   TCanvas *cenergyResolution2D = new TCanvas("cenergyResolution2D","cenergyResolution2D",800,800);
//   cenergyResolution2D->SetLeftMargin(0.16);
//   spectrum2d->Draw("LEGO2");  
//   cenergyResolution2D->Print("energyResolution2D.png");
//   cc->Print("energyResolution.png");
//   delete spectrum2d;
//   delete all;
//   delete central;
//   delete hs;
//   delete cc;
//   delete c;
//   delete legend;  
//   
//   
//   
//   //-----------------------------//
//   //         Photopeaks
//   //-----------------------------//
//   TCanvas *c2 = (TCanvas*) gDirectory->Get("Photopeak positions vs. i,j");
//   
//   spectrum2d = (TH2F*)c2->GetPrimitive("Photopeak positions vs. i,j");
//   all = new TH1F("all","",100,0,12000);
//   all->GetXaxis()->SetTitle("Photopeak Position [ADC Channels]");
//   //all->GetXaxis()->SetLabelSize(0.5);
//   all->SetName("511 KeV Photopeak Position");
// //   all->GetYaxis()->SetTitle("N");
//   all->SetStats(0);
//   central = new TH1F("central","",100,0,12000);
// //   all->GetYaxis()->SetTitle("N");
//   central->SetStats(0);
//   central->SetFillStyle(3001);
//   central->SetFillColor(kRed);
//   double sumLO = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	central->Fill(spectrum2d->GetBinContent(i+1,j+1));
// 	sumLO += spectrum2d->GetBinContent(i+1,j+1);
//       }
//       all->Fill(spectrum2d->GetBinContent(i+1,j+1));
//     } 
//   }
//   hs = new THStack("hs","");
//   double averageLO = sumLO/16 ;
//   double stdLO = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	stdLO += (1.0/16.0)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageLO,2));
//       }
//     }
//   }
//   
//   cc = new TCanvas("cc","",1200,800);
//   all->SetFillStyle(3001);
//   all->SetFillColor(kBlue);
//   cc->cd();
//   hs->Add(all);
//   hs->Add(central);
// //   all->Draw();
// //   central->Draw("same");
//   legend = new TLegend(0.15,0.62,0.54,0.89,"");
//   legend->SetFillStyle(0);
//   legend->AddEntry(all,"All channels","f");
//   legend->AddEntry(central,"Central channels","f");
//   hs->Draw("nostack");
//   hs->GetXaxis()->SetTitle("ADC Channels");
//   hs->GetXaxis()->SetTitleOffset(1);
//   hs->GetXaxis()->SetTitleSize(0.045);
//   hs->GetXaxis()->SetLabelSize(0.045);
//   hs->GetYaxis()->SetLabelSize(0.045);
//   //hs->GetXaxis()->SetTickSize(0.08);
//   legend->Draw();
//   
//   spectrum2d->GetYaxis()->SetTitleOffset(1.8);
//   spectrum2d->GetXaxis()->SetTitleOffset(1.5);
//   spectrum2d->GetZaxis()->SetTitleOffset(2.1);
//   spectrum2d->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d->GetZaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d->GetZaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
//   spectrum2d->GetYaxis()->SetNdivisions(8,2,0, kTRUE);
//   
//   TCanvas *clightOutput2D = new TCanvas("clightOutput2D","clightOutput2D",800,800);
//   clightOutput2D->SetLeftMargin(0.18);
//   spectrum2d->Draw("LEGO2");   
//   clightOutput2D->Print("lightOutput2D.png");
//   cc->Print("lightOutput.png");
//   delete spectrum2d;
//   delete all;
//   delete central;
//   delete hs;
//   delete cc;
//   delete c2;
//   delete legend;
//   
//   
//   //-----------------------------//
//   //      En Res Corrected
//   //-----------------------------//
//   c = (TCanvas*) gDirectory->Get("Corrected Energy res FWHM vs. i,j");
//   
//   spectrum2d = new TH2F();
//   spectrum2d = (TH2F*)c->GetPrimitive("Corrected Energy res FWHM vs. i,j");
//   all = new TH1F("Corrected Energy res FWHM","",100,0,1);
//   all->GetXaxis()->SetTitle("Corrected Energy Resolution FWHM");
//   all->SetName("Corrected Energy Resolution FWHM");
// //   all->GetYaxis()->SetTitle("N");
//   all->SetStats(0);
//   central = new TH1F("central","",100,0,1);
// //   all->GetYaxis()->SetTitle("N");
//   central->SetStats(0);
//   central->SetFillStyle(3001);
//   central->SetFillColor(kRed);
//   double sumEN_corr = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	central->Fill(spectrum2d->GetBinContent(i+1,j+1));
// 	sumEN_corr += spectrum2d->GetBinContent(i+1,j+1);
//       }
//       all->Fill(spectrum2d->GetBinContent(i+1,j+1));
//     } 
//   }
//   hs = new THStack("hs","");
//   double averageEN_corr = sumEN_corr/16 ;
//   double stdEN_corr = 0;
//   for(int i = 0 ; i < 8 ; i++)
//   {
//     for(int j = 0 ; j < 8 ; j++)
//     {
//       if(i > 1 && i < 6 && j > 1 && j < 6)
//       {
// 	stdEN_corr += (1.0/16.0)*sqrt(pow(spectrum2d->GetBinContent(i+1,j+1)-averageEN_corr,2));
//       }
//     }
//   }
//   cc = new TCanvas("cc","",1200,800);
//   all->SetFillStyle(3001);
//   all->SetFillColor(kBlue);
//   cc->cd();
//   hs->Add(all);
//   hs->Add(central);
// //   all->Draw();
// //   central->Draw("same");
//   legend = new TLegend(0.5,0.62,0.893,0.89,"");
//   legend->SetFillStyle(0);
//   legend->AddEntry(all,"All channels","f");
//   legend->AddEntry(central,"Central channels","f");
//   hs->Draw("nostack");
//   hs->GetXaxis()->SetTitle("Energy Resolution FWHM");
//   hs->GetXaxis()->SetTitleOffset(1);
//   hs->GetXaxis()->SetTitleSize(0.045);
//   hs->GetXaxis()->SetLabelSize(0.045);
//   hs->GetYaxis()->SetLabelSize(0.045);
//   legend->Draw();
//   
//   spectrum2d->GetYaxis()->SetTitleOffset(1.8);
//   spectrum2d->GetXaxis()->SetTitleOffset(1.5);
//   spectrum2d->GetZaxis()->SetTitleOffset(1.7);
//   spectrum2d->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d->GetZaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d->GetZaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetNdivisions(8,2,0, kTRUE);
//   spectrum2d->GetYaxis()->SetNdivisions(8,2,0, kTRUE);
//   
//   TCanvas *cenergyResolution2DCorrected = new TCanvas("cenergyResolution2DCorrected","cenergyResolution2DCorrected",800,800);
//   cenergyResolution2DCorrected->SetLeftMargin(0.16);
//   spectrum2d->Draw("LEGO2");
//   cenergyResolution2DCorrected->Print("energyResolution2DCorrected.png");
//   cc->Print("energyResolutionCorrected.png");
//   delete spectrum2d;
//   delete all;
//   delete central;
//   delete hs;
//   delete cc;
//   delete legend;
//   
//   
//   std::cout << "En. Res Central Channels [FWHM] = " << sumEN/16 << " +/- " << stdEN << std::endl;
//   std::cout << "En. Res Corrected Central Channels [FWHM] = " << sumEN_corr/16 << " +/- " << stdEN_corr << std::endl;
//   std::cout << "Light Output Central Channels [ADC Ch.] = " << sumLO/16 << " +/- " << stdLO << std::endl;
//   
//   
//   // crystal 14 - example plots
//   f->cd("Module 0.0");
//   gDirectory->cd("MPPC D1 - 0.0-0.3");
//   gDirectory->cd("Crystal 14");
//   
//   //ADCvsW
//   c = new TCanvas();
//   c = (TCanvas*) gDirectory->Get("Complete ADC channels vs. W - Crystal 14");
//   spectrum2d = new TH2F();
//   spectrum2d = (TH2F*)c->GetPrimitive("Complete ADC channels vs. W - Crystal 14");
//   spectrum2d->SetTitle("");
//   spectrum2d->GetYaxis()->SetTitleOffset(1.5);
//   spectrum2d->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetTitleSize(0.045);
//   TCanvas *cADCvsW = new TCanvas("cADCvsW","cADCvsW",1200,800);
//   cADCvsW->SetRightMargin(0.13);
//   cADCvsW->SetLeftMargin(0.13);
//   spectrum2d->Draw("COLZ");
//   cADCvsW->Print("ADCvsW.png");
//   delete spectrum2d;
//   delete c;
//   
// //   f->cd("Module 0.0");
// //   gDirectory->cd("MPPC D1 - 0.0-0.3");
// //   gDirectory->cd("Crystal 14");
//   
//   //ADCvsW_2D
//   c = (TCanvas*) gDirectory->Get("ADC channels vs. W - Crystal 14");
//   spectrum2d = new TH2F();
//   spectrum2d = (TH2F*)c->GetPrimitive("ADC channels vs. W - Crystal 14");
//   c->SetLeftMargin(0.13);
//   spectrum2d->SetTitle("");
//   spectrum2d->GetYaxis()->SetTitleOffset(2.6);
//   spectrum2d->GetXaxis()->SetTitleOffset(1.9);
//   spectrum2d->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d->GetZaxis()->SetLabelSize(0.045);
//   spectrum2d->GetYaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d->GetXaxis()->SetNdivisions(4,2,0, kTRUE);
//   spectrum2d->GetYaxis()->SetNdivisions(4,2,0, kTRUE);
//   spectrum2d->GetXaxis()->SetRangeUser(0.4,0.56);
//   spectrum2d->GetYaxis()->SetRangeUser(5000,11000);
//   spectrum2d->Draw("LEGO2");
//   c->Print("ADCvsW_2D.png");
//   delete spectrum2d;
//   delete c;
//   
//   
//   //PeaksVsW
//   TCanvas *cP = (TCanvas*) gDirectory->Get("spectrum2d_1");
//   TH1D *spectrum2d_1 = new TH1D();
//   TF1 *fit = new TF1();
//   spectrum2d_1 = (TH1D*)cP->GetPrimitive("spectrum2d_1");
//   fit = (TF1*) cP->GetPrimitive("linearCrystal");
//   spectrum2d_1->SetTitle("");
//   spectrum2d_1->GetYaxis()->SetTitle("ADC channels");
//   spectrum2d_1->GetYaxis()->SetTitleOffset(1.1);
//   spectrum2d_1->GetXaxis()->SetTitleOffset(1);
//   spectrum2d_1->GetXaxis()->SetTitleSize(0.045);
//   spectrum2d_1->GetYaxis()->SetTitleSize(0.045);
//   spectrum2d_1->GetXaxis()->SetLabelSize(0.045);
//   spectrum2d_1->GetYaxis()->SetLabelSize(0.045);
//   spectrum2d_1->GetXaxis()->SetRangeUser(0.42,0.55);
//   spectrum2d_1->GetYaxis()->SetRangeUser(7200,8800);
//   spectrum2d_1->Draw();
//   fit->Draw("same");
//   cP->Print("PeaksVsW.png");
//   
//   //SpectraComparison
//   
//   TCanvas *cSpOriginal = new TCanvas();
//   TCanvas *cSpCorrected = new TCanvas();
//   TH1F *spOriginal = new TH1F();
//   TH1F *spCorrected = new TH1F();
//   TF1 *fitOriginal = new TF1();
//   TF1 *fitCorrected = new TF1();
//   TCanvas *cComparison= new TCanvas();
//   double enresOriginal,enresCorrected;
//   
//   // cry 14
//   cSpOriginal  = (TCanvas*) gDirectory->Get("Charge Spectrum - Crystal 14 - MPPC D1");
//   cSpCorrected = (TCanvas*) gDirectory->Get("Charge Spectrum Corrected - Crystal 14");
//   spOriginal   = (TH1F*) cSpOriginal->GetPrimitive("Charge Spectrum - Crystal 14 - MPPC D1");
//   spCorrected  = (TH1F*) cSpCorrected->GetPrimitive("Charge Spectrum Corrected - Crystal 14");
//   fitOriginal  = (TF1*) cSpOriginal->GetPrimitive("gauss");
//   fitCorrected = (TF1*) cSpCorrected->GetPrimitive("gauss_corr");
//   cComparison = new TCanvas("cComparison","cComparison",1200,800);
//   cComparison->cd();
//   spOriginal->Scale(1/spOriginal->GetBinContent(spOriginal->GetMaximumBin()));
//   spCorrected->Scale(1/spCorrected->GetBinContent(spCorrected->GetMaximumBin()));
//   spOriginal->SetFillStyle(3001);
//   spOriginal->SetFillColor(kBlue);
//   spOriginal->SetTitle("");
//   spOriginal->GetYaxis()->SetTitle("A. U.");
//   spOriginal->GetYaxis()->SetTitleOffset(0.7);
//   spOriginal->GetXaxis()->SetTitleOffset(1);
//   spOriginal->GetXaxis()->SetTitleSize(0.045);
//   spOriginal->GetYaxis()->SetTitleSize(0.045);
//   spOriginal->GetXaxis()->SetLabelSize(0.045);
//   spOriginal->GetYaxis()->SetLabelSize(0.045);
//   spCorrected->SetFillStyle(3001);
//   spCorrected->SetFillColor(kRed);
//   enresOriginal = 100.0 * (2.355 * fitOriginal->GetParameter(2) / fitOriginal->GetParameter(1));
//   enresCorrected = 100.0 * (2.355 * fitCorrected->GetParameter(2) / fitCorrected->GetParameter(1));
//   std::stringstream legendOriginal,legendCorrected;
//   legendOriginal <<  "En. Res. (FWHM) = " << std::setprecision(3) << enresOriginal << "%";
//   legendCorrected << "En. Res. (FWHM) = "<< std::setprecision(3) << enresCorrected << "%";
//   legend = new TLegend(0.45,0.7,0.893,0.89,"");
//   legend->SetFillStyle(0);
//   legend->AddEntry(spOriginal,legendOriginal.str().c_str(),"f");
//   legend->AddEntry(spCorrected,legendCorrected.str().c_str(),"f");
//   spOriginal->Draw();
//   spCorrected->Draw("same");
//   legend->Draw("same");
//   cComparison->Print("SpectraComparison14.png");
//   cSpOriginal->Print("SampleSpectrum14.png");
//   delete legend;
//   delete cSpOriginal,cSpCorrected,spOriginal,spCorrected,fitOriginal,fitCorrected,cComparison;
//   legendOriginal.str("");
//   legendCorrected.str("");
//   
//   //cry 27
//   // crystal 14 - example plots
//   f->cd("Module 0.0");
//   gDirectory->cd("MPPC B2 - 0.0-1.1");
//   gDirectory->cd("Crystal 27");
//   cSpOriginal  = (TCanvas*) gDirectory->Get("Charge Spectrum - Crystal 27 - MPPC B2");
//   cSpCorrected = (TCanvas*) gDirectory->Get("Charge Spectrum Corrected - Crystal 27");
//   spOriginal   = (TH1F*) cSpOriginal->GetPrimitive("Charge Spectrum - Crystal 27 - MPPC B2");
//   spCorrected  = (TH1F*) cSpCorrected->GetPrimitive("Charge Spectrum Corrected - Crystal 27");
//   fitOriginal  = (TF1*) cSpOriginal->GetPrimitive("gauss");
//   fitCorrected = (TF1*) cSpCorrected->GetPrimitive("gauss_corr");
//   cComparison = new TCanvas("cComparison","cComparison",1200,800);
//   cComparison->cd();
//   spOriginal->Scale(1/spOriginal->GetBinContent(spOriginal->GetMaximumBin()));
//   spCorrected->Scale(1/spCorrected->GetBinContent(spCorrected->GetMaximumBin()));
//   spOriginal->SetFillStyle(3001);
//   spOriginal->SetFillColor(kBlue);
//   spOriginal->SetTitle("");
//   spOriginal->GetYaxis()->SetTitle("A. U.");
//   spOriginal->GetYaxis()->SetTitleOffset(0.7);
//   spOriginal->GetXaxis()->SetTitleOffset(1);
//   spOriginal->GetXaxis()->SetTitleSize(0.045);
//   spOriginal->GetYaxis()->SetTitleSize(0.045);
//   spOriginal->GetXaxis()->SetLabelSize(0.045);
//   spOriginal->GetYaxis()->SetLabelSize(0.045);
//   spCorrected->SetFillStyle(3001);
//   spCorrected->SetFillColor(kRed);
//   enresOriginal = 100.0 * (2.355 * fitOriginal->GetParameter(2) / fitOriginal->GetParameter(1));
//   enresCorrected = 100.0 * (2.355 * fitCorrected->GetParameter(2) / fitCorrected->GetParameter(1));
// //   std::stringstream legendOriginal,legendCorrected;
//   legendOriginal <<  "En. Res. (FWHM) = " << std::setprecision(3) << enresOriginal << "%";
//   legendCorrected << "En. Res. (FWHM) = "<< std::setprecision(3) << enresCorrected << "%";
//   legend = new TLegend(0.45,0.7,0.893,0.89,"");
//   legend->SetFillStyle(0);
//   legend->AddEntry(spOriginal,legendOriginal.str().c_str(),"f");
//   legend->AddEntry(spCorrected,legendCorrected.str().c_str(),"f");
//   spOriginal->Draw();
//   spCorrected->Draw("same");
//   legend->Draw("same");
//   cComparison->Print("SpectraComparison27.png");
//   cSpOriginal->Print("SampleSpectrum27.png");
//   delete legend;
//   
//   
//   
//   //floodTop
//   f->cd("Module 0.0");
//   TCanvas *cFlood = (TCanvas*) gDirectory->Get("Flood Histogram 2D");
//   cFlood->SetRightMargin(0.15);
//   cFlood->Print("floodTop.png");
//   
//   return 0;  
}