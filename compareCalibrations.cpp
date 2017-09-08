// g++ -o ../build/compareCalibrations compareCalibrations.cpp `root-config --cflags --glibs`
// -Wl,--no-as-needed -lHist -lCore -lMathCore

//FIXME offset could be inserted directly from command line?

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

//class of input points from doi tag bench
class inputDoi_t
{
public:
  int i;
  int j;
  double m;
  double q;
  double doires;
  double avgs;
  std::vector<double> w;
  std::vector<double> sw;
  std::vector<double> sqrt_nentries;
  std::vector<double> z;
  std::vector<double> sz;

  int pointsFromDoi;
  inputDoi_t(int a){ pointsFromDoi = a;};
  void clear()
  {
    w.clear();
    sw.clear();
    sqrt_nentries.clear();
    z.clear();
    sz.clear();
  };

  friend std::istream& operator>>(std::istream& input, inputDoi_t& s)
  {
    input >> s.i;
    input >> s.j;
    input >> s.m;
    input >> s.q;
    input >> s.doires;
    input >> s.avgs;
    for(int p = 0; p < s.pointsFromDoi; p++)
    {
      double wValue,swValue,sqrtValue;
      input >> wValue >> swValue >> sqrtValue; /*>> zValue >> szValue*/
      s.w.push_back(wValue);
      s.sw.push_back(swValue);
      s.sqrt_nentries.push_back(sqrtValue);
      //       s.z.push_back(zValue);
      //       s.sz.push_back(szValue);
    }
    return input;
  }
};

// class of real z positions of tag points (after check of alignment)
class inputZ_t
{
public:
  int pointsFromDoi;
  int i;
  int j;
  std::vector<double> z;
  inputZ_t(int a){ pointsFromDoi = a;};
  void clear(){z.clear();};

  friend std::istream& operator>>(std::istream& input, inputZ_t& s)
  {
    input >> s.i; //read i
    input >> s.j;           //read
    for(int p = 0; p < s.pointsFromDoi; p++)
    {
      double zValue;
      input >> zValue;
      s.z.push_back(zValue);
    }
    return input;
  }

};

// double z[pointsFromDoi] = {10.8,8,5.2,2.4};
// double sz[pointsFromDoi] = {0.5,0.5,0.5,0.5};


struct data_t
{
  std::string file;
  std::string title;
  int color;
  data_t (std::string a, std::string b, int c):file(a),title(b),color(c){}
};

struct cry_res_t
{


};


int main(int argc, char **argv)
{
  if(argc < 6)
  {
    std::cout << "USAGE:\t\t compareCalibrations calibration_params.txt zPositions.txt pointsFromDoi moduleCalibrationFile1.root title1 [moduleCalibrationFile2.root title2] ... [moduleCalibrationFileN.root titleN]" << std::endl;
    std::cout << std::endl;
    return 1;
  }

  //doi tag file
  std::vector<inputDoi_t> inputDoi;
  std::ifstream fDoiTag;
  fDoiTag.open(argv[1],std::ios::in);

  std::ifstream fZPos;
  fZPos.open(argv[2],std::ios::in);
  std::vector<inputZ_t> inputZ;

  int pointsFromDoi = atoi(argv[3]);

  std::vector<data_t> dataset;
  int color = 2;
  for (int i = 4 ; i < argc ; i = i+2)
  {
    std::string file,title;
    file = argv[i];
    title = argv[i+1];

    data_t a(file,title,color);
    dataset.push_back(a);
    color++;
  }
  //   std::cout << "aaa" << std::endl;
  //   data_t a("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/background_good/Run_2016_13_04_11_37_39/nearChannels.root","background",2);
  //   data_t b("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/lateral/Run_2016_05_04_08_48_47/nearChannels.root","lateral",3);
  //   data_t c("/home/marco/Universita/NewClearPEM/Experiments/HardDepo2/calibration_paper/top/Run_2016_05_04_11_59_20/nearChannels.root","top",4);
  //
  //   dataset.push_back(a);
  //   dataset.push_back(b);
  //   dataset.push_back(c);

  std::vector<TFile*> file;
  for(int i = 0 ; i < dataset.size(); i++)
  {
    TFile *temp = new TFile(dataset[i].file.c_str());
    file.push_back(temp);
  }


  //   while(!fZPos.eof()) //file z_positions_etc
  //   {
  //     int a,b;
  //     inputZ_t tempInput(pointsFromDoi);
  //
  //     fZPos >> tempInput.i >> tempInput.j;
  //     for(int k = 0 ; k < pointsFromDoi ; k++)
  //     {
  //       double value;
  //       fZPos >> value;
  //       tempInput.z.push_back(value);
  //     }
  //
  // //     if(!fZPos.eof())
  // //     {
  //       inputZ.push_back(tempInput);
  // //     }
  //   }
  inputZ_t tempInput(pointsFromDoi);
  while(fZPos >> tempInput)
  {
    inputZ.push_back(tempInput);
    tempInput.clear();
  }

  //   double TagOffset = 1.4;
  //   while(!fDoiTag.eof()) //file calibration_params_etc.. One line per crystal
  //   {
  //     int a,b;
  //     inputDoi_t tempInput;
  //
  //
  //     fDoiTag >> tempInput.i >> tempInput.j >>tempInput.m >> tempInput.q >> tempInput.doires >> tempInput.avgs;
  //
  //     //     tempInput.q = tempInput.q /*+ TagOffset*/; //until we properly characterize the alignment of the tagging setup and modify the z positions at source, we use this temporary evaluation of the offset between the aligment of the x-y-z stage scale that we thought was correct and the one that we derive from the fine analisys of the tagging setup
  //
  //     for(int k = 0 ; k < pointsFromDoi ; k++)
  //     {
  //       double w,sw,sqrt_nentries;
  //       fDoiTag >> w >> sw >> sqrt_nentries;
  //       tempInput.w.push_back(w);
  //       tempInput.sw.push_back(sw);
  //       tempInput.sqrt_nentries.push_back(sqrt_nentries);
  //       //tempInput.z[k] = z[k] + TagOffset; //until we properly characterize the alignment of the tagging setup and modify the z positions at source, we use this temporary evaluation of the offset between the aligment of the x-y-z stage scale that we thought was correct and the one that we derive from the fine analisys of the tagging setup
  //       //tempInput.sz[k] = sz[k];
  //     }
  //
  //
  //     for(int k = 0; k < inputZ.size(); k++)
  //     {
  //       if(tempInput.i == inputZ[k].i && tempInput.j == inputZ[k].j)
  //       {
  // 	for(int h = 0 ; h < pointsFromDoi ; h++)
  // 	{
  // 	  tempInput.z.push_back(inputZ[k].z[h]);
  // 	  tempInput.sz.push_back(0.2); //FIXME what is the correct value?
  // 	}
  //       }
  //     }
  //
  //
  //     if(!fDoiTag.eof())
  //     {
  //       inputDoi.push_back(tempInput);
  //     }
  //   }

  //

  inputDoi_t tempInputDoi(pointsFromDoi);
  while(fDoiTag >> tempInputDoi)
  {
    for(int k = 0; k < inputZ.size(); k++)
    {
      if(tempInputDoi.i == inputZ[k].i && tempInputDoi.j == inputZ[k].j)
      {
        for(int h = 0 ; h < pointsFromDoi ; h++)
        {
          tempInputDoi.z.push_back(inputZ[k].z[h]);
          tempInputDoi.sz.push_back(0.2); //FIXME what is the correct value?
        }
      }
    }
    inputDoi.push_back(tempInputDoi);
    tempInputDoi.clear();
  }

  fDoiTag.close();




  //DEBUG
  std::cout << "InputZ -------------" << std::endl;
  for(int i = 0; i < inputZ.size(); i++)
  {
    std::cout << inputZ[i].i << " " << inputZ[i].j << " ";
    for(int j = 0 ; j < inputZ[i].z.size(); j++)
    std::cout << inputZ[i].z[j] << " ";
    std::cout << std::endl;
  }
  std::cout << "-------------" << std::endl;
  std::cout << std::endl;
  std::cout << "InputDoi ------------" << std::endl;
  for(int i = 0; i < inputDoi.size(); i++)
  {
    std::cout << inputDoi[i].i << " " << inputDoi[i].j << " " << inputDoi[i].m << " " << inputDoi[i].q << " " << inputDoi[i].doires << " " << inputDoi[i].avgs << " ";
    for(int j = 0 ; j < inputDoi[i].z.size(); j++)
    std::cout << inputDoi[i].w[j] << " "<< inputDoi[i].sw[j] << " " << inputDoi[i].sqrt_nentries[j] << " "<< inputDoi[i].z[j] << " " << inputDoi[i].sz[j] << " ";
    std::cout << std::endl;
  }
  std::cout << "-------------" << std::endl;
  //----DEBUG

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
  //   int nCurves = nmodulex*nmoduley*nmppcx*nmppcy*ncrystalsx*ncrystalsy;
  int nCurves = inputDoi.size();
  //   std::vector<TCanvas *> C_curves;
  //   std::vector<TMultiGraph *> curves;
  //   for(int i = 0 ; i < nCurves ; i++)
  //   {
  //
  //   }

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



  //   int curveNumber = 0;
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


              TGraphErrors *points = NULL;
              TMultiGraph *mg = NULL;
              TCanvas *temp = NULL;
              //draw the points
              for(int k = 0; k < inputDoi.size(); k++)
              {
                int i = inputDoi[k].i;
                int j = inputDoi[k].j;
                if(iCrystal == i && jCrystal == j)
                {
                  std::cout << i << " " << j << " " << crystalNumber<< std::endl;
                  //prepare the temp canvas and multigraph
                  std::stringstream title;
                  title << "Curves - Crystal " << crystalNumber;
                  temp = new TCanvas(title.str().c_str(),title.str().c_str(),1200,800);
                  // 		  C_curves.push_back(temp);
                  mg = new TMultiGraph();
                  // 		  curves.push_back(mg);
                  points = new TGraphErrors(pointsFromDoi,&inputDoi[k].w[0],&inputDoi[k].z[0],&inputDoi[k].sw[0],&inputDoi[k].sz[0]);
                  points->SetLineColor(kBlack);
                  points->SetLineWidth(2);
                  points->SetMarkerStyle(21);
                  points->SetMarkerSize(1);
                  points->SetMarkerColor(kBlack);
                  mg->Add(points,"P");
                  legend->AddEntry(points,"tagging bench","lep");
                }
              }
              // 	      std::cout << "aaaaaaa" << std::endl;
              //if no points where found, stop
              if(points != NULL)
              {
                //draw the curves
                for(int i = 0 ; i < dataset.size(); i++)
                {

                  file[i]->cd();
                  //go to the crystal
                  bool dir_exists =file[i]->cd(directory.str().c_str());
                  TCanvas* C_graph;
                  TGraph *calibGraph;

                  if(dir_exists)
                  {
                    //take the graph
                    std::stringstream stream;
                    stream << "Calibration Plot - Crystal " << crystalNumber;
                    C_graph = (TCanvas*) gDirectory->Get(stream.str().c_str());
                    if(C_graph)
                    calibGraph = (TGraph*) C_graph->GetPrimitive(stream.str().c_str());

                    for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
                    {

                      double xPoint;
                      double yPoint;
                      points->GetPoint(pointNum,xPoint,yPoint);
                      // 		  std::cout << "aaa" << std::endl;
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
                    mg->Add(calibGraph,"L");
                    legend->AddEntry(calibGraph,dataset[i].title.c_str(),"l");

                    //numb of events
                    // 		    stream.str("");
                    // 		    stream << "Charge Spectrum Corrected - Crystal " << crystalNumber;
                    // 		    TCanvas *C_canvas;
                    // 		    TH1F* histo;
                    // 		    C_canvas = (TCanvas*) gDirectory->Get(stream.str().c_str());
                    // 		    if(C_graph)
                    // 		      histo = (TH1F*) C_graph->GetPrimitive(stream.str().c_str());
                    // 		    if(histo)
                    // 		    {
                    //
                    // 		    }
                    //numb of events in p.e. peak

                    //energy res

                    //peak position



                  }

                }
                fOut->cd();
                temp->cd();
                // 		std::cout << directory.str() << std::endl;
                // 		gPad->Update();


                mg->Draw("a");
                std::stringstream stream;
                stream << "Calibrations comparison - Crystal " << crystalNumber;


                mg->SetTitle(stream.str().c_str());
                mg->GetXaxis()->SetTitle("W");
                mg->GetYaxis()->SetTitle("Z [mm]");
                legend->Draw();


                temp->Write();


                //comparisons part
                for(int i = 0 ; i < dataset.size() - 1; i++)
                {
                  for(int j = i+1 ; j < dataset.size(); j++)
                  {
                    std::cout << dataset[i].title << " - " << dataset[j].title << std::endl;
                  }
                }

                // 		curveNumber++;
              }
            }
          }
        }
      }
    }
  }

  //   std::cout << "aaaa" << std::endl;

  //   for(int i = 0 ; i < nCurves ; i++)
  //   {
  //     curves[i]->Write();
  //   }
  //


  TGraph *mean = new TGraph();

  fOut->cd();
  for(int i = 0 ; i < dataset.size() ; i++)
  {

    histDoiPrecision[i]->Write();
    histDoiPrecision_central[i]->Write();
    mean->SetPoint(i,i,histDoiPrecision_central[i]->GetMean());
  }
  mean->Write();

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
