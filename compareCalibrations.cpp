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

struct data_t
{
  std::string file;
  std::string title;
  int color;
  data_t (std::string a, std::string b, int c):file(a),title(b),color(c){}
};

void usage()
{
  std::cout << "\t\t" << "[ -k <calibration_params file> ] " << std::endl
            << "\t\t" << "[ -z <z_positions file> ] " << std::endl
            << "\t\t" << "[ --fileTop <ModuleCalibraton output file for TOP irradiation> ] " << std::endl
            << "\t\t" << "[ --fileLat <ModuleCalibraton output file for LATERAL irradiation> ] " << std::endl
            << "\t\t" << "[ --fileBg  <ModuleCalibraton output file for BACKGROUND irradiation> ] " << std::endl
            << "\t\t" << "[ --points <points from doi scan> ] " << std::endl
            << "\t\t" << "[ --length <length of crystals in mm> ] " << std::endl
            << "\t\t" << "[ --nmodulex <number of modules in x> ] " << std::endl
            << "\t\t" << "[ --nmoduley <number of modules in y> ] " << std::endl
            << "\t\t" << "[ --nmppcx <number of mppc in x PER MODULE> ] " << std::endl
            << "\t\t" << "[ --nmppcy <number of mppc in y PER MODULE> ] " << std::endl
            << "\t\t" << "[ --ncrystalsx <number of crystals in x PER MPPC> ] " << std::endl
            << "\t\t" << "[ --ncrystalsy <number of crystals in y PER MPPC> ] " << std::endl
            << "\t\t" << "[ --doiPrecisionMin  <doi precision histograms: min> ] " << std::endl
            << "\t\t" << "[ --doiPrecisionMax  <doi precision histograms: max> ] " << std::endl
            << "\t\t" << "[ --doiPrecisionBins <doi precision histograms: bins> ] " << std::endl
            << "\t\t" << "[ --zPrecision <precision in z bench positioning, in mm> ] " << std::endl
            << "\t\t" << std::endl;
}


int main(int argc, char **argv)
{
  if(argc < 2)
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
    usage();
    return 1;
  }

  std::ifstream fDoiTag;  // doi tag calibration_params file
  std::ifstream fZPos;    // z_positions file

  std::string doiTagFileName;
  std::string zPosFileName;
  std::string topFileName;
  std::string latFileName;
  std::string bgFileName;


  //default values
  float doiPrecisionMin = -3.0;
  float doiPrecisionMax = 3.0;
  int doiPrecisionBins = 30;
  float zPrecision = 0.2;
  int inputFiles = 0;
  int pointsFromDoi = 1;
  double crystaLenght = 15.0;
  int nmodulex = 1;
  int nmoduley = 1;
  int nmppcx = 4;
  int nmppcy = 4;
  int ncrystalsx = 1;
  int ncrystalsy = 1;
  std::string letter[14] = {"A","B","C","D","E","F","G","H","I","J","K","L","M","N"};  //standard ordering
  std::string number[14] = {"1","2","3","4","5","6","7","8","9","10","11","12","13","14"};  //standard ordering

  //input module calibration files
  std::vector<data_t> dataset;
  // int color = 2;
  // for (int i = 5 ; i < argc ; i = i+2)
  // {
  //   std::string file,title;
  //   file = argv[i];
  //   title = argv[i+1];
  //
  //   data_t a(file,title,color);
  //   dataset.push_back(a);
  //   color++;
  // }
  int color[] = {2,4,8,28,30,46};
  int colorCounter = 0;
  static struct option longOptions[] =
  {
      { "fileTop", required_argument, 0, 0 },
      { "fileLat", required_argument, 0, 0 },
      { "fileBg", required_argument, 0, 0 },
      { "points", required_argument, 0, 0 },
			{ "length", required_argument, 0, 0 },
      { "nmodulex", required_argument, 0, 0 },
      { "nmoduley", required_argument, 0, 0 },
      { "nmppcx", required_argument, 0, 0 },
      { "nmppcy", required_argument, 0, 0 },
      { "ncrystalsx", required_argument, 0, 0 },
      { "ncrystalsy", required_argument, 0, 0 },
      { "doiPrecisionMin", required_argument, 0, 0 },
      { "doiPrecisionMax", required_argument, 0, 0 },
      { "doiPrecisionBins", required_argument, 0, 0 },
      { "zPrecision", required_argument, 0, 0},
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "k:z:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'k'){
			doiTagFileName = (char *)optarg;
    }
		else if (c == 'z'){
      zPosFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 0){
      topFileName = (char *)optarg;
      data_t a(topFileName,"Top irradiation",color[colorCounter]);
      dataset.push_back(a);
      colorCounter++;
      // inputFiles++;
    }
    else if (c == 0 && optionIndex == 1){
      latFileName = (char *)optarg;
      data_t a(latFileName,"Lateral irradiation",color[colorCounter]);
      dataset.push_back(a);
      colorCounter++;
      // inputFiles++;
    }
    else if (c == 0 && optionIndex == 2){
      bgFileName = (char *)optarg;
      data_t a(bgFileName,"Lu background",color[colorCounter]);
      dataset.push_back(a);
      colorCounter++;
      // inputFiles++;
    }
		else if (c == 0 && optionIndex == 3){
      pointsFromDoi = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 4){
      crystaLenght = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 5){
      nmodulex = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 6){
      nmoduley = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      nmppcx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      nmppcy = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      ncrystalsx = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      ncrystalsy = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      doiPrecisionMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      doiPrecisionMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      doiPrecisionBins = atoi((char *)optarg);
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}

  if(dataset.size() == 0)
  {
    std::cout	<< "Usage: " << argv[0] << std::endl;
    usage();
    return 1;
  }

  fDoiTag.open(doiTagFileName.c_str(),std::ios::in);
  fZPos.open(zPosFileName.c_str(),std::ios::in);


  std::vector<TFile*> file;
  for(int i = 0 ; i < dataset.size(); i++)
  {
    TFile *temp = new TFile(dataset[i].file.c_str());
    file.push_back(temp);
  }

  //input z positions file
  std::vector<inputZ_t> inputZ;
  inputZ_t tempInput(pointsFromDoi);
  while(fZPos >> tempInput)
  {
    inputZ.push_back(tempInput);
    tempInput.clear();
  }


  //doi tag file
  std::vector<inputDoi_t> inputDoi;
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
          tempInputDoi.sz.push_back(zPrecision); //FIXME what is the correct value?
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


  // const int nmodulex = 1;
  // const int nmoduley = 1;
  // const int nmppcx = 4;
  // const int nmppcy = 4;
  // const int ncrystalsx = 1;
  // const int ncrystalsy = 1;
  // std::string letter[4] = {"A","B","C","D"};
  // std::string number[4] = {"1","2","3","4"};


  int nCurves = inputDoi.size();


  //doi precision
  std::vector<TH1F*> histDoiPrecision;
  std::vector<TH1F*> histDoiPrecision_central;
  TH1F* histLinePrecision = new TH1F("Line precision","Line precision",doiPrecisionBins,doiPrecisionMin,doiPrecisionMax);
  histLinePrecision->SetLineColor(kBlack);
  histLinePrecision->SetFillColor(kBlack);
  histLinePrecision->SetFillStyle(3001);
  histLinePrecision->SetLineWidth(0);
  TH1F* histLinePrecision_central= new TH1F("Line precision - Central Crystals","Line precision - Central Crystals",doiPrecisionBins,doiPrecisionMin,doiPrecisionMax);
  histLinePrecision_central->SetLineWidth(0);
  histLinePrecision_central->SetLineColor(kBlack);
  histLinePrecision_central->SetFillColor(kBlack);
  histLinePrecision_central->SetFillStyle(3001);
  for(int i = 0 ; i < dataset.size(); i++)
  {
    std::stringstream title;
    title << "Calibration Precision - " << dataset[i].title;
    TH1F* temp = new TH1F(title.str().c_str(),title.str().c_str(),doiPrecisionBins,doiPrecisionMin,doiPrecisionMax);
    temp->SetLineWidth(3);
    temp->SetLineColor(dataset[i].color);
    histDoiPrecision.push_back(temp);


    title.str("");
    title << "Calibration Precision, Central Crystals - " << dataset[i].title;
    temp = new TH1F(title.str().c_str(),title.str().c_str(),doiPrecisionBins,doiPrecisionMin,doiPrecisionMax);
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
              int crystalNumber = (cryI*ncrystalsy*nmppcy + cryJ);
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
              legend = new TLegend(0.5,0.55,0.89,0.89,"");
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
                  points->SetMarkerSize(1.2);
                  points->SetMarkerColor(kBlack);
                  TF1 *lineFit = new TF1("lineFit","[0]*x + [1]",0,1);
                  points->Fit(lineFit,"RNQ");
                  for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
                  {
                    double xPoint;
                    double yPoint;
                    points->GetPoint(pointNum,xPoint,yPoint);
                    histLinePrecision->Fill(lineFit->Eval(xPoint)- yPoint);
                    if(iCrystal > (ncrystalsx -1) && iCrystal < (ncrystalsx*nmppcx - ncrystalsx) && jCrystal  > (ncrystalsy-1) && jCrystal < (ncrystalsy*nmppcy - ncrystalsy))
                    {
                      histLinePrecision_central->Fill(lineFit->Eval(xPoint)-yPoint);
                    }

                  }
                  // points->GetPoint(pointNum,xPoint,yPoint);
                  // histLinePrecision->Fill();
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
                      histDoiPrecision[i]->Fill(calibGraph->Eval(xPoint)-yPoint);
                      if(iCrystal > (ncrystalsx -1) && iCrystal < (ncrystalsx*nmppcx - ncrystalsx) && jCrystal  > (ncrystalsy-1) && jCrystal < (ncrystalsy*nmppcy - ncrystalsy))
                      {
                        histDoiPrecision_central[i]->Fill(calibGraph->Eval(xPoint)-yPoint);
                      }

                    }


                    //inclinations and offsets in case of "line"
                    // double crystaLenght = 15.0;
                    double minBound = -1;
                    double maxBound = -1;
                    for(int step = 0 ; step < 10000 ; step ++)
                    {
                      if(minBound == -1)
                      {
                        if(calibGraph->Eval((crystaLenght/10000)*step) < (crystaLenght - 1))
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
                    // std::cout << minBound << " " << maxBound << std::endl;
                    TF1 *lineFitSigma = new TF1("lineFitSigma","[0]*x + [1]",minBound,maxBound);
                    lineFitSigma->SetLineColor(1);
                    // 		  calibGraph->Draw("AL");
                    calibGraph->Fit(lineFitSigma,"RNQ");


                    // for(int pointNum = 0; pointNum < pointsFromDoi ; pointNum++)
                    // {
                    //
                    //   double xPoint;
                    //   double yPoint;
                    //   points->GetPoint(pointNum,xPoint,yPoint);
                    //   histLinePrecision[i]->Fill(lineFitSigma->Eval(xPoint)-yPoint);
                    //   if(iCrystal > 1 && iCrystal < 6 && jCrystal  > 1 && jCrystal < 6)
                    //   {
                    //     histLinePrecision_central[i]->Fill(lineFitSigma->Eval(xPoint)-yPoint);
                    //   }
                    //
                    // }


                    //add the graph to the multigraph
                    calibGraph->SetLineColor(dataset[i].color);
                    calibGraph->SetLineWidth(2);
                    mg->Add(calibGraph,"L");
                    legend->AddEntry(calibGraph,dataset[i].title.c_str(),"l");
                  }
                }
                fOut->cd();
                temp->cd();
                // 		std::cout << directory.str() << std::endl;
                // 		gPad->Update();


                mg->Draw("a");
                std::stringstream stream;
                stream << "Calibrations comparison - Crystal " << crystalNumber;


                // mg->SetTitle(stream.str().c_str());
                mg->GetXaxis()->SetTitle("W");
                mg->GetYaxis()->SetTitle("Z [mm]");
                mg->GetXaxis()->SetLabelSize(0.05);
                mg->GetYaxis()->SetLabelSize(0.05);
                mg->GetXaxis()->SetTitleSize(0.06);
                mg->GetYaxis()->SetTitleSize(0.06);
                mg->GetXaxis()->SetTitleOffset(0.8);
                mg->GetYaxis()->SetTitleOffset(0.7);

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
              }
            }
          }
        }
      }
    }
  }

  TGraph *mean = new TGraph();

  fOut->cd();
  for(int i = 0 ; i < dataset.size() ; i++)
  {

    histDoiPrecision[i]->Write();
    histDoiPrecision_central[i]->Write();
    mean->SetPoint(i,i,histDoiPrecision_central[i]->GetMean());
  }

  histLinePrecision->Write();
  histLinePrecision_central->Write();
  mean->Write();

  TCanvas *C_allHisto = new TCanvas("Precision Comparison","Precision Comparison",1200,800);
  TLegend* legend;
  legend = new TLegend(0.65,0.7,0.89,0.89,"");

  legend->SetFillStyle(0);

  THStack *hs = new THStack("hs","");
  hs->Add(histLinePrecision);
  legend->AddEntry(histLinePrecision,"line fit","l");
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


  TCanvas *C_allHisto_central = new TCanvas("Precision Comparison_central","",1200,800);
  TLegend* legend_central;
  legend_central = new TLegend(0.65,0.7,0.89,0.89,"");
  legend_central->SetFillStyle(0);
  THStack *hs_central = new THStack("hs_central","");
  hs_central->Add(histLinePrecision_central);
  legend_central->AddEntry(histLinePrecision_central,"line fit","f");
  for(int i = 0 ; i < dataset.size() ; i++)
  {
    hs_central->Add(histDoiPrecision_central[i]);
    legend_central->AddEntry(histDoiPrecision_central[i],dataset[i].title.c_str(),"l");
  }

  // hs_central->SetTitle("Precision Comparison_central");
  hs_central->Draw("nostack");
  hs_central->GetXaxis()->SetTitle("delta [mm]");
  hs_central->GetXaxis()->SetLabelSize(0.05);
  hs_central->GetYaxis()->SetLabelSize(0.05);
  hs_central->GetXaxis()->SetTitleSize(0.06);
  hs_central->GetYaxis()->SetTitleSize(0.06);
  hs_central->GetXaxis()->SetTitleOffset(0.8);
  hs_central->GetYaxis()->SetTitleOffset(0.7);
  legend_central->Draw();
  fOut->cd();
  C_allHisto_central->Write();



  fOut->Close();
  return 0;
}
