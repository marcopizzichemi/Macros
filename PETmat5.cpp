//----------------------------------------------------------//
//                                                          //
//  PROGRAM FOR ANALYSING THE "MATRIX 5" OUTPUT (SIM)       //
//                                                          //
//----------------------------------------------------------//

// compile with 
// g++ -o PETmat5 PETmat5.cpp `root-config --cflags --glibs` -lSpectrum -lMLP -lTreePlayer

// The program takes the output of the simulation, creates a TTree (adding some brach to make analysis easier)
// and tries to identify the x y z coordinates of the gamma interaction


//TODO



#include "TROOT.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
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



//----------------------------------------------------------//
//  Main program                                            //
//----------------------------------------------------------//
int main (int argc, char** argv)
{
  
  //----------------------------------------------------------//
  //  Input args                                              //
  //----------------------------------------------------------//
  std::string ListFileName = argv[1];
  int channel = 5;
  //define names for output
  std::stringstream SoutputFile;
  SoutputFile << ListFileName.substr(0,ListFileName.length()-4) << "_ch" << channel << ".root"; 
  std::string fileRoot = ListFileName.substr(0,ListFileName.length()-4) + ".root";
  std::string outputFile = SoutputFile.str();
  
  //----------------------------------------------------------//
  //  Create the TTree                                        //
  //----------------------------------------------------------//
  long long int counter = 0;
  TTree *t1 = new TTree("adc","adc"); //creates the tree
  //leaf names
  int charge[16];
  std::stringstream snames[16];
  std::stringstream stypes[16];
  std::string names[16];
  std::string types[16];
  float x,y,z;
  float energy;
  float doi;
  int TriggerChannel;
  float floodx,floody,rowsum,columnsum,total;
  //mppc positions
  float xmppc[16] = {4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8,4.8,1.6,-1.6,-4.8};
  float ymppc[16] = {4.8,4.8,4.8,4.8,1.6,1.6,1.6,1.6,-1.6,-1.6,-1.6,-1.6,-4.8,-4.8,-4.8,-4.8};
  for (int i = 0 ; i < 16 ; i++)
  {
    charge[i] = 0;
    snames[i] << "ch" << i;
    stypes[i] << "ch" << i << "/I";
    names[i] = snames[i].str();
    types[i] = stypes[i].str();
    t1->Branch(names[i].c_str(),&charge[i],types[i].c_str());
  }
  t1->Branch("x",&x,"x/F");
  t1->Branch("y",&y,"y/F");
  t1->Branch("z",&z,"z/F");
  t1->Branch("energy",&energy,"energy/F");
  t1->Branch("doi",&doi,"doi/F");
  t1->Branch("TriggerChannel",&TriggerChannel,"TriggerChannel/I");
  t1->Branch("floodx",&floodx,"floodx/F");
  t1->Branch("floody",&floody,"floody/F");
  ifstream file(ListFileName.c_str());
  while(!file.eof())
  {
    total = 0;
    rowsum = 0;
    columnsum = 0;
    int maxCharge = 0;
    for (int i = 0 ; i < 16 ; i++)
    {
      file >> charge[i];
      if (charge[i] > maxCharge)
      {
	maxCharge = charge[i];
	TriggerChannel = i;
      }
      total += charge[i];
      rowsum += charge[i]*xmppc[i];
      columnsum += charge[i]*ymppc[i];
      //cout << charge[i] << " ";
    }
    file >> x >> y >> z >> energy >> doi;
    floodx=(rowsum/total);
    floody=(columnsum/total);
    //cout << x << y << z << energy << doi << TriggerChannel << endl;
    if(!file.eof())
    {
      t1->Fill();//fills the tree with the 4 values of charge deposited
      counter++;
      if( (counter % 1000) == 0 )
	std::cout << counter << std::endl;
    }
  }
  
  
  //----------------------------------------------------------//
  //  Analysis pt.1: find the interaction x and y             //
  //----------------------------------------------------------//
  //channel/crystal mapping
  // 3  2  1  0
  // 7  6  5  4
  // 11 10 9  8
  // 15 14 13 12
  //bottom part of the matrix, columns are merged
  //upper part of the matrix, rows are merged
  
  //first, the evergreen anger logic plot
  TH2F *Anger = new TH2F("Anger","Flood histogram of all events",200,-2.8,2.8,200,-2.8,2.8);
  Anger->GetXaxis()->SetTitle("Anger X [mm]");
  Anger->GetYaxis()->SetTitle("Anger Y [mm]");
  t1->Draw("floody:floodx >> Anger");
  //then, see it's x "profile"
  TH1F *Check2 = new TH1F("Check2","Flood X",400,-7,7);
  //TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
  Check2->GetXaxis()->SetTitle("Flood X [mm]");
  Check2->GetYaxis()->SetTitle("Counts");
  t1->Draw("floodx >> Check2");
  
  //now find the 4 peaks in this histo  
  TSpectrum *SpectrumAnger;
  SpectrumAnger = new TSpectrum(4);
  Int_t Peaks = SpectrumAnger->Search(Check2,1,"",0.5); //TODO pass to "goff"
  Float_t *PeaksPos = SpectrumAnger->GetPositionX();
  Float_t *PeaksH = SpectrumAnger->GetPositionY();
  for(int i = 0; i < Peaks ; i ++)
  {
    std::cout << PeaksPos[i] << std::endl;
  }
  //sort the peaks
  float swap,swapY;
  for ( int jSwap = 0 ; jSwap < Peaks - 1 ; jSwap++)
  {
    for(int kSwap = (jSwap+1); kSwap < Peaks; kSwap++)   // rest of the elements
    {
      if (PeaksPos[jSwap] > PeaksPos[kSwap])             // asc order
      {
	swap = PeaksPos[jSwap];
	swapY = PeaksH[jSwap];
	PeaksPos[jSwap] = PeaksPos[kSwap];
	PeaksH[jSwap] = PeaksH[kSwap];
	PeaksPos[kSwap] = swap;
	PeaksH[kSwap] = swapY;
      }
    }
  }
  //check ordering
  for(int i = 0; i < Peaks ; i ++)
  {
    std::cout << PeaksPos[i] << std::endl;
  }
  
  
//   TF1 *FourGauss = new TF1("FourGauss",  "gaus(0)+gaus(3)+gaus(6)+gaus(9)",-2,2);
//   for(int i = 0; i < Peaks ; i ++)
//   {
//     FourGauss->SetParameter(i*3,PeaksH[i]);
//     FourGauss->SetParameter((i+1)*3,PeaksPos[i]);
// //     FourGauss->SetParameter(i*3,PeaksH[i]);
//   }
//   Check2->Fit("FourGauss","","",-2,2);
  
  //fit them and store mean and sigma in a vector
  std::vector<float> mean,sigma;
  float distance;
  for(int i = 0; i < Peaks ; i ++)
  {
    if(i < (Peaks-1))
    {
      distance = (PeaksPos[i] - PeaksPos[i+1])/2.0;
      TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",PeaksPos[i]-0.25*distance,PeaksPos[i]+0.25*distance);
      gauss->SetParameter(0,1500);
      gauss->SetParameter(1,PeaksPos[i]);
      gauss->SetParameter(2,0.2);
      Check2->Fit("gauss","R+");
      mean.push_back(gauss->GetParameter(1));
      sigma.push_back(gauss->GetParameter(2));
    }
    else
    {
      distance = (PeaksPos[i-1] - PeaksPos[i])/2.0;
      TF1 *gauss = new TF1("gauss",  "[0]*exp(-0.5*((x-[1])/[2])**2)",PeaksPos[i]-0.25*distance,PeaksPos[i]+0.25*distance);
      gauss->SetParameter(0,1500);
      gauss->SetParameter(1,PeaksPos[i]);
      gauss->SetParameter(2,0.2);
      Check2->Fit("gauss","R+");
      mean.push_back(gauss->GetParameter(1));
      sigma.push_back(gauss->GetParameter(2));
    }
  }
  
  TCut CutCrystalD[16];
  TCut CutCrystalU[16];
  //now for the 4 columns, we would need to find and separate the 4 crystals
  //for the moment, being this a simulation, we can plot the wieghted sums on y and see if there is a correlation with the real y
  TCanvas *C_ycorr = new TCanvas("C_ycorr","",1200,800);
  C_ycorr->Divide(2,2);
  for(int i = 0; i < Peaks ; i ++)
  {
    //careful: i is the column starting 0 from left,
    // so i is ok for the mean and sigma order, but is reversed for the ch num
    int MPPCcolNum = (3 - i);
    C_ycorr->cd(i+1);
    TCut columnCut;
    //streams and string to create the variables
    std::stringstream columnCutStream,columnTitleStream,columnVarStream;
    std::string columnCutString,columnTitleString,columnVarString;
    //histogram title
    columnTitleStream << "Weighted Sum on Column MPPCs vs. Real Y - Column  " << i;
    columnTitleString = columnTitleStream.str();
    //variable. the mppcs are the ones in the Nth column, starting from the right
    columnVarStream << "(";
    for(int j=0 ; j<4 ; j++)
    {
      columnVarStream << "ch" << (j*4)+MPPCcolNum << "*" << ymppc[(j*4)+MPPCcolNum];
      if (j!=3) columnVarStream << "+";
    }
    columnVarStream << ")/(";
    for(int j=0 ; j<4 ; j++)
    {
      columnVarStream << "ch" << (j*4)+MPPCcolNum ;
      if (j!=3) columnVarStream << "+";
    }
    columnVarStream << ")";
    //we use them as base for the "crystal" identification
    //CutCrystalD[(j*4)+MPPCcolNum] = columnVarStream.str();
    //CutCrystalU[(j*4)+MPPCcolNum] = columnVarStream.str();
    
    columnVarStream << ": y >>" << columnTitleString.c_str();
    columnVarString = columnVarStream.str();
    
    //Cut
    columnCutStream << "floodx > " << mean[i] - 2.0*sigma[i] << " && " << "floodx < " << mean[i] + 2.0*sigma[i];
    columnCutString = columnCutStream.str();
    
    std::cout << columnTitleString << std::endl;
    std::cout << columnVarString << std::endl;
    std::cout << columnCutString << std::endl;
    std::cout << std::endl;
    
    TH2F *WSum1 = new TH2F(columnTitleString.c_str(),columnTitleString.c_str(),200,-7,7,200,-7,7);
    WSum1->GetXaxis()->SetTitle("Real Y [mm]");
    WSum1->GetYaxis()->SetTitle("W Sum [mm]");
    t1->Draw(columnVarString.c_str(),columnCutString.c_str());
    
    TLine *line[5];
    line[0] = new TLine(-6.4,-6.4,-6.4,6.4);
    line[1] = new TLine(-3.2,-6.4,-3.2,6.4);
    line[2] = new TLine(0,-6.4,0,6.4);
    line[3] = new TLine(3.2,-6.4,3.2,6.4);
    line[4] = new TLine(6.4,-6.4,6.4,6.4);
    
    
    for ( int iLine = 0 ; iLine < 5 ; iLine++)
    {
      line[iLine]->SetLineColor(kRed);
      line[iLine]->SetLineWidth(2);
      line[iLine]->Draw();
    }
    
    
  }
  
  //from this we can have the "crystal" separation points
  
//   CutCrystalD[] = ;
//   CutCrystalU[] =;
  
  
  
  
  
  
  new TCanvas;
  
  
  
  
  
  
  
  
  
  
  new TCanvas;
  
  //P.S. anyway, we can separate the columns even just by the distance.
  TH2F *WSum = new TH2F("WSum","Weighted Sum on Column MPPCs vs. Real Y - column3",200,-7,7,200,-7,7);
  //TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
  WSum->GetXaxis()->SetTitle("Real Y [mm]");
  WSum->GetYaxis()->SetTitle("W Sum [mm]");
  t1->Draw("(ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13) : y >> WSum","floodx >0.2 && floodx < 0.8");
  //WSum->Draw(); 
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  
  //Interaction crystal
  TCanvas *Canvas = new TCanvas("Canvas","",1200,800);
  Canvas->Divide(3,2);
  //first, plot the sum of all channels
  Canvas->cd(1);
  TH1F *SumPlot = new TH1F("SumPlot","Energy spectrum - Sum of all channels",200,0,16000);
  SumPlot->GetXaxis()->SetTitle("Ch Sum");
  SumPlot->GetYaxis()->SetTitle("Counts");
  t1->Draw("ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 >> SumPlot");
  SumPlot->Draw();  
  //then, the sum vs doi (and there is no dependance, that's ok
  Canvas->cd(2);
  TH2F *SumVsDoi = new TH2F("SumVsDoi","DOI vs. Sum of all channels",200,0,16000,200,0,16);
  SumVsDoi->GetXaxis()->SetTitle("Ch Sum");
  SumVsDoi->GetYaxis()->SetTitle("DOI [mm]");
  TH2F *SumVsDoiCut = new TH2F("SumVsDoiCut","DOI vs. Sum of all channels (Cut)",200,0,16000,200,0,16);
  t1->Draw("doi : ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 >> SumVsDoi");
  t1->Draw("doi : ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 >> SumVsDoiCut","energy >0.66","COLZ");
  SumVsDoi->Draw();
  SumVsDoiCut->Draw("COLZ same");
  //so the cut can be defined on the total charge of the 16 ch
  TCut energyCutFromTotal = "ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch14+ch15 > 13000";
  //now, the anger logic of all the events
  Canvas->cd(3);
  //TH2F *Anger = new TH2F("Anger","Flood histogram of all events",200,-2.8,2.8,200,-2.8,2.8);
  TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
  //Anger->GetXaxis()->SetTitle("Anger X [mm]");
  //Anger->GetYaxis()->SetTitle("Anger Y [mm]");
  //t1->Draw("floody:floodx >> Anger");
  t1->Draw("floody:floodx >> AngerCut","floodx >0.2 && floodx < 0.8","COLZ same ");
  Anger->Draw();
  AngerCut->Draw("COLZ same");
  Canvas->cd(4);
  TH1F *Check22 = new TH1F("Check22","Flood X",200,-7,7);
  //TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
  Check22->GetXaxis()->SetTitle("Flood X [mm]");
  Check22->GetYaxis()->SetTitle("Counts");
  t1->Draw("floodx >> Check22");
  Check22->Draw();  
  //Weighted sum for the y coordinate
  //now, check that events are in the right place
  Canvas->cd(5);
  TH1F *Check = new TH1F("Check","Flood Y for column 3",200,-7,7);
  //TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
  Check->GetXaxis()->SetTitle("Flood Y [mm]");
  Check->GetYaxis()->SetTitle("Counts");
  t1->Draw("floody >> Check","floodx >0.2 && floodx < 0.8");
  Check->Draw();  
  //correlation on x?
  Canvas->cd(6);
//   TH2F *WSum = new TH2F("WSum","Weighted Sum on Column MPPCs vs. Real Y - column3",200,-7,7,200,-7,7);
//   //TH2F *AngerCut = new TH2F("AngerCut","Flood histogram of all events (cut)",200,-2.8,2.8,200,-2.8,2.8);
//   WSum->GetXaxis()->SetTitle("Real Y [mm]");
//   WSum->GetYaxis()->SetTitle("W Sum [mm]");
//   t1->Draw("(ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13) : y >> WSum","floodx >0.2 && floodx < 0.8");
   WSum->Draw();  
  
  
  
  TCut crystal;
  TCut crystal5 = "(0 < x ) && ( x < 3.2 ) && (0  < y) && ( y < 3.2 )";
  
  TCut PEenergy = "energy > 0.66";
  TCut DOIup =  "doi > 8.25";
  TCut DOIdown = "doi < 8.25";
  
  //define useful combined TCuts
  TCut cry5PE = crystal5 + PEenergy;
  TCut cry5PEup = cry5PE + DOIup;
  TCut cry5PEdown = cry5PE + DOIdown;
   
  
  
  
  
  //look for channel 5
  TCanvas *Csearch = new TCanvas("Csearch","",1200,800);
  Csearch->Divide(2);
  Csearch->cd(1);
  TH2F *WSumCut = new TH2F("WSumCut","Weighted Sum on Column MPPCs vs. Real Y (Cut)",200,-7,7,200,-7,7);
  t1->Draw("(ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13) : y >> WSumCut","floodx >0.2 && floodx < 0.8" + crystal5);
  WSum->Draw();
  WSumCut->Draw("COLZ same");

  //you get from 0.11 to 1.8 from this plot  
  TCut CutCol3 = "floodx >0.2 && floodx < 0.8";
  TCut yCut = "(((ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13)) > - 0.07) && (((ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13)) < 1.67)";
  TCut crystalIdentified = CutCol3 + yCut;
  
  //and check
  Csearch->cd(2);
  TH2F *WSumRealCut = new TH2F("WSumRealCut","Weighted Sum on Column MPPCs vs. Real Y (Real Cut)",200,-7,7,200,-7,7);
  t1->Draw("(ch1*6.4+ch5*3.2-ch9*3.2-ch13*6.4)/(ch1+ch5+ch9+ch13) : y >> WSumRealCut",crystalIdentified);
  WSum->Draw();
  WSumRealCut->Draw("COLZ same");
  
  
  //anger logic doesn't give much. It can easily distinguish the column, 
  //so let's say that we found out third col (crystal 5), now we can analyze only the other 3 cols
  
//   TCut CutCol3 = " (floodx > 0.25 ) && (floodx < 0.9 )";
//   TCanvas *flood = new TCanvas("flood","",1200,800);
//   flood->Divide(2,3);
//   flood->cd(1);
//   TH2F *scatter = new TH2F("scatter","Anger Logic",200,-3,3,200,-3,3);
//   t1->Draw("floody:floodx >> scatter","","COLZ");
//   flood->cd(3);
//   TH1F *row1 = new TH1F("row1","Row 1",200,0,1); 
//   TH1F *row2 = new TH1F("row2","Row 2",200,0,1); 
//   TH1F *row3 = new TH1F("row3","Row 3",200,0,1); 
//   TH1F *row4 = new TH1F("row4","Row 4",200,0,1); 
//   t1->Draw("(ch0+ch2+ch3)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row1",CutCol3);
//   flood->cd(4);
//   t1->Draw("(ch7+ch6+ch4)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row2",CutCol3);
//   flood->cd(5);
//   t1->Draw("(ch11+ch10+ch8)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15 )>> row3",CutCol3);
//   flood->cd(6);
//   t1->Draw("(ch12+ch14+ch15)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row4",CutCol3);
//   flood->cd(2);
  //TCut CutCol3 = " (floodx > 0.25 ) && (floodx < 0.9 )";
  TCanvas *flood = new TCanvas("flood","",1200,800);
  flood->Divide(2,3);
  flood->cd(1);
  TH2F *scatter = new TH2F("scatter","Anger Logic",200,-7,7,200,-7,7);
  t1->Draw("floody:floodx >> scatter", crystalIdentified,"COLZ");
  TLine *line[10];
  line[0] = new TLine(-6.4,-6.4,-6.4,6.4);
  line[1] = new TLine(-3.2,-6.4,-3.2,6.4);
  line[2] = new TLine(0,-6.4,0,6.4);
  line[3] = new TLine(3.2,-6.4,3.2,6.4);
  line[4] = new TLine(6.4,-6.4,6.4,6.4);
  line[5] = new TLine(-6.4,-6.4,6.4,-6.4);
  line[6] = new TLine(-6.4,-3.2,6.4,-3.2);
  line[7] = new TLine(-6.4,0,6.4,0);
  line[8] = new TLine(-6.4,3.2,6.4,3.2);
  line[9] = new TLine(-6.4,6.4,6.4,6.4);
  
  for ( int iLine = 0 ; iLine < 10 ; iLine++)
  {
    line[iLine]->SetLineColor(kRed);
    line[iLine]->SetLineWidth(2);
    line[iLine]->Draw();
  }
  
  
  flood->cd(3);
  //TH2F *h1b = new TH2F("h1b","points",4,-7,7,4,-7,7);
  
  TH1F *row1 = new TH1F("row1","Row 1",200,0,1); 
  TH1F *row2 = new TH1F("row2","Row 2",200,0,1); 
  TH1F *row3 = new TH1F("row3","Row 3",200,0,1); 
  TH1F *row4 = new TH1F("row4","Row 4",200,0,1); 
  t1->Draw("(ch1)/(ch1+ch5+ch9+ch13) >> row1",CutCol3);
  flood->cd(4);
  t1->Draw("(ch5)/(ch1+ch5+ch9+ch13) >> row2",CutCol3);
  flood->cd(5);
  t1->Draw("(ch9)/(ch1+ch5+ch9+ch13)>> row3",CutCol3);
  flood->cd(6);
  t1->Draw("(ch13)/(ch1+ch5+ch9+ch13) >> row4",CutCol3);
  flood->cd(2);
  
  TLegend *legend0 = new TLegend(0.13,0.7,0.39,0.89,"");
  legend0->SetFillStyle(0);
  
  row1->SetLineColor(2);
  row1->SetFillColor(2);
  row1->SetFillStyle(3001);
  legend0->AddEntry(row1,"Row 1","f");
  row1->Draw();
  
  row2->SetLineColor(3);
  row2->SetFillColor(3);
  row2->SetFillStyle(3001);
  legend0->AddEntry(row2,"Row 2","f");
  row2->Draw("same");
  
  row3->SetLineColor(4);
  row3->SetFillColor(4);
  row3->SetFillStyle(3001);
  legend0->AddEntry(row3,"Row 3","f");
  row3->Draw("same");
  
  row4->SetLineColor(5);
  row4->SetFillColor(5);
  row4->SetFillStyle(3001);
  legend0->AddEntry(row4,"Row 4","f");
  row4->Draw("same");
  legend0->Draw();
  
  //same but using only PE events
//   TCut CutCol3Plus = CutCol3 + crystal5 + PEenergy;
//   
//   TCanvas *floodPE = new TCanvas("floodPE","",1200,800);
//   floodPE->Divide(2,3);
//   floodPE->cd(1);
//   TH2F *scatterPE = new TH2F("scatterPE","Anger Logic - Only PE",200,-3,3,200,-3,3);
//   t1->Draw("floody:floodx >> scatterPE",CutCol3Plus,"COLZ");
//   floodPE->cd(3);
//   TH1F *row1PE = new TH1F("row1PE","Row 1  - Only PE",200,0,1); 
//   TH1F *row2PE = new TH1F("row2PE","Row 2  - Only PE",200,0,1); 
//   TH1F *row3PE = new TH1F("row3PE","Row 3  - Only PE",200,0,1); 
//   TH1F *row4PE = new TH1F("row4PE","Row 4  - Only PE",200,0,1); 
//   t1->Draw("(ch0+ch2+ch3)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row1PE",CutCol3Plus);
//   floodPE->cd(4);
//   t1->Draw("(ch7+ch6+ch4)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row2PE",CutCol3Plus);
//   floodPE->cd(5);
//   t1->Draw("(ch11+ch10+ch8)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15 )>> row3PE",CutCol3Plus);
//   floodPE->cd(6);
//   t1->Draw("(ch12+ch14+ch15)/(ch0+ch1+ch2+ch3+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> row4PE",CutCol3Plus);
//   floodPE->cd(2);
//   
//   
//   TLegend *legend0PE = new TLegend(0.13,0.7,0.39,0.89,"");
//   legend0PE->SetFillStyle(0);
//   
//   row1PE->SetLineColor(2);
//   row1PE->SetFillColor(2);
//   row1PE->SetFillStyle(3001);
//   legend0PE->AddEntry(row1PE,"Row 1","f");
//   row1PE->Draw();
//   
//   row2PE->SetLineColor(3);
//   row2PE->SetFillColor(3);
//   row2PE->SetFillStyle(3001);
//   legend0PE->AddEntry(row2PE,"Row 2","f");
//   row2PE->Draw("same");
//   
//   row3PE->SetLineColor(4);
//   row3PE->SetFillColor(4);
//   row3PE->SetFillStyle(3001);
//   legend0PE->AddEntry(row3PE,"Row 3","f");
//   row3PE->Draw("same");
//   
//   row4PE->SetLineColor(5);
//   row4PE->SetFillColor(5);
//   row4PE->SetFillStyle(3001);
//   legend0PE->AddEntry(row4PE,"Row 4","f");
//   row4PE->Draw("same");
//   legend0PE->Draw();
  
  
  
  TCanvas *IntCry = new TCanvas("IntCry","",1200,800);
  IntCry->Divide(2,2);
  TH1F *bareSpectrum = new TH1F("bareSpectrum","Spectrum of Channel 5",200,0,5000); 
  TH1F *RealCr5 = new TH1F("RealCr5","Real Interactions in Channel 5",200,0,5000);
  TH1F *OnlyPE = new TH1F("OnlyPE","Selecting only PE in bare spectrum",200,0,5000); 
  
  //TH1F *colDoiDown = new TH1F("colDoiDown","Tr.Ch/Col",100,0,1);
  
  IntCry->cd(1);
  t1->Draw("ch5 >> bareSpectrum");
  IntCry->cd(3);
  t1->Draw("ch5 >> RealCr5",crystal5);
  IntCry->cd(2);
  TCut photopeak = "ch5 > 1000";
  t1->Draw("ch5 >> OnlyPE",photopeak);
  TCut intersect = photopeak + crystal5;
  IntCry->cd(4);
  t1->Draw("ch5",intersect);
//  
  
  
  //DOI
  //Using the info from the lines
  TH1F *rowDoiUp = new TH1F("rowDoiUp","Tr.Ch/Row",100,0,1); 
  TH1F *rowDoiDown = new TH1F("rowDoiDown","Tr.Ch/Row",100,0,1);
  TH1F *colDoiUp = new TH1F("colDoiUp","Tr.Ch/Col",100,0,1); 
  TH1F *colDoiDown = new TH1F("colDoiDown","Tr.Ch/Col",100,0,1);
  //Using the info from lines / all
  TH1F *rowDivAllDoiUp = new TH1F("rowDivAllDoiUp","Row/All",100,0,1); 
  TH1F *rowDivAllDoiDown = new TH1F("rowDivAllDoiDown","Row/All",100,0,1);
  TH1F *colDivAllDoiUp = new TH1F("colDivAllDoiUp","Column/All",100,0,1); 
  TH1F *colDivAllDoiDown = new TH1F("colDivAllDoiDown","Column/All",100,0,1);
  
  
  
  TCanvas *LinesCanvas = new TCanvas("Lines","",1200,800);
  LinesCanvas->Divide(2,2);
  
  LinesCanvas->cd(1);
  TLegend *legend1 = new TLegend(0.13,0.7,0.39,0.89,"");
  legend1->SetFillStyle(0);
  rowDoiDown->SetLineColor(2);
  rowDoiDown->SetFillColor(2);
  rowDoiDown->SetFillStyle(3001);
  legend1->AddEntry(rowDoiDown,"Bottom DOI","f");
  rowDoiUp->SetLineColor(3);
  rowDoiUp->SetFillColor(3);
  rowDoiUp->SetFillStyle(3001);
  legend1->AddEntry(rowDoiUp,"Top DOI","f");
  t1->Draw("ch5/(ch4+ch5+ch6+ch7) >> rowDoiDown",cry5PEup);
  t1->Draw("ch5/(ch4+ch5+ch6+ch7) >> rowDoiUp",cry5PEdown,"same");
  legend1->Draw();
  
  LinesCanvas->cd(3);
  TLegend *legend2 = new TLegend(0.13,0.7,0.39,0.89,"");
  legend2->SetFillStyle(0);
  colDoiDown->SetLineColor(2);
  colDoiDown->SetFillColor(2);
  colDoiDown->SetFillStyle(3001);
  legend2->AddEntry(colDoiDown,"Bottom DOI","f");
  colDoiUp->SetLineColor(3);
  colDoiUp->SetFillColor(3);
  colDoiUp->SetFillStyle(3001);
  legend2->AddEntry(colDoiUp,"Top DOI","f");
  t1->Draw("ch5/(ch1+ch5+ch9+ch13) >> colDoiUp",cry5PEup);
  t1->Draw("ch5/(ch1+ch5+ch9+ch13) >> colDoiDown",cry5PEdown,"same");
  legend2->Draw();

//   TCanvas *Canvas2 = new TCanvas("Two","",1200,800);
//   Canvas2->Divide(2,2);
  
  LinesCanvas->cd(4);
  TLegend *legend3 = new TLegend(0.13,0.7,0.39,0.89,"");
  legend3->SetFillStyle(0);
  rowDivAllDoiDown->SetLineColor(2);
  rowDivAllDoiDown->SetFillColor(2);
  rowDivAllDoiDown->SetFillStyle(3001);
  legend3->AddEntry(rowDivAllDoiDown,"Bottom DOI","f");
  rowDivAllDoiUp->SetLineColor(3);
  rowDivAllDoiUp->SetFillColor(3);
  rowDivAllDoiUp->SetFillStyle(3001);
  legend3->AddEntry(rowDivAllDoiUp,"Top DOI","f");
  t1->Draw("(ch4+ch5+ch6+ch7)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> rowDivAllDoiUp",cry5PEup);
  t1->Draw("(ch4+ch5+ch6+ch7)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> rowDivAllDoiDown",cry5PEdown,"same");
  legend3->Draw();
  
  
  LinesCanvas->cd(2);
  TLegend *legend4 = new TLegend(0.13,0.7,0.39,0.89,"");
  legend4->SetFillStyle(0);
  colDivAllDoiDown->SetLineColor(2);
  colDivAllDoiDown->SetFillColor(2);
  colDivAllDoiDown->SetFillStyle(3001);
  legend4->AddEntry(colDivAllDoiDown,"Bottom DOI","f");
  colDivAllDoiUp->SetLineColor(3);
  colDivAllDoiUp->SetFillColor(3);
  colDivAllDoiUp->SetFillStyle(3001);
  legend4->AddEntry(colDivAllDoiUp,"Top DOI","f");
  t1->Draw("(ch1+ch5+ch9+ch13)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> colDivAllDoiUp",cry5PEup);
  t1->Draw("(ch1+ch5+ch9+ch13)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) >> colDivAllDoiDown",cry5PEdown,"same");
  legend4->Draw();
  
  
  TCanvas *due = new TCanvas("due","",1200,800);
  due->Divide(2,2);
  
  TH2F *rowDivAll = new TH2F("rowDivAll","Row/All vs. DOI",100,0,16,100,0.1,0.6);
  TH2F *colDivAll = new TH2F("colDivAll","Column/All vs. DOI",100,0,16,100,0.1,0.6);
  TH2F *rowTimesColDivAll = new TH2F("rowTimesColDivAll","(Row*Column)/(All^2) vs. DOI",100,0,16,100,0,0.4);
  //TH2F *rowDivAll new TH2F("rowDivAll","Row/All vs. DOI",100,0,1,100,0,1);
  
  due->cd(1);
  t1->Draw("(ch4+ch5+ch6+ch7)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) : doi >> rowDivAll",crystalIdentified + energyCutFromTotal);
  due->cd(3);
  t1->Draw("(ch1+ch5+ch9+ch13)/(ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15) : doi >> colDivAll",crystalIdentified + energyCutFromTotal);
  due->cd(2);
  t1->Draw("(ch1+ch5+ch9+ch13)*(ch4+ch5+ch6+ch7)/((ch0+ch1+ch2+ch3+ch4+ch5+ch6+ch7+ch8+ch9+ch10+ch11+ch12+ch13+ch15)^2) : doi >> rowTimesColDivAll",crystalIdentified + energyCutFromTotal);
  
  TFile* fOut = new TFile(outputFile.c_str(),"recreate");
  fOut->cd();
  Anger->Write();
  C_ycorr->Write();
  Check2->Write();
  Canvas->Write();
  Csearch->Write();
  due->Write();
  flood->Write();
  IntCry->Write();
  LinesCanvas->Write();
  fOut->Close();
  
  TFile* fTree = new TFile(fileRoot.c_str(),"recreate");
  fTree->cd();
  t1->Write();
  fTree->Close();
  
}