// compile with
// g++ -o ../build/timeAnalysis timeAnalysis.cpp `root-config --cflags --glibs` -Wl,--no-as-needed -lHist -lCore -lMathCore -lTree -lTreePlayer -lgsl -lgslcblas

// small program to extract timing calibration and data

#include "TROOT.h"
#include "TFile.h"
#include "TStyle.h"
#include "TSystem.h"
#include "TLegend.h"
#include "TCanvas.h"
#include "TH1F.h"
#include "TH2F.h"
#include "TH3I.h"
#include "TString.h"
#include "TApplication.h"
#include "TLegend.h"
#include "TTree.h"
#include "TF2.h"
#include "TGraph2D.h"
#include "TGraph.h"
#include "TSpectrum.h"
#include "TSpectrum2.h"
#include "TTreeFormula.h"
#include "TMath.h"
#include "TChain.h"
#include "TCut.h"
#include "TLine.h"
#include "TError.h"
#include "TEllipse.h"
#include "TFormula.h"
#include "TGraphErrors.h"
#include "TGraph2DErrors.h"
#include "TMultiGraph.h"
#include "TCutG.h"
#include "TGaxis.h"
#include "TPaveStats.h"
#include "TProfile.h"
#include "TH1D.h"
#include "TPaveText.h"
#include "TGraphDelaunay.h"
#include "TVector.h"
#include "TNamed.h"
#include "TPaveLabel.h"
#include "THStack.h"
#include "TFitResult.h"
#include "TMatrixD.h"


#include <iostream>
#include <fstream>
#include <string>
#include <sstream>
#include <stdio.h>
#include <getopt.h>
#include <algorithm>    // std::sort
#include <gsl/gsl_matrix_double.h>
#include <gsl/gsl_linalg.h>

#include <sys/types.h>
#include <dirent.h>

// typedef std::vector<std::string> stringvec;
// list files in directory
// taken from
// http://www.martinbroadhurst.com/list-the-files-in-a-directory-in-c.html
void read_directory(const std::string& name, std::vector<std::string> &v)
{
    DIR* dirp = opendir(name.c_str());
    struct dirent * dp;
    while ((dp = readdir(dirp)) != NULL) {
        v.push_back(dp->d_name);
    }
    closedir(dirp);
}


gsl_matrix *
invert_a_matrix(gsl_matrix *matrix,size_t size)
{
    gsl_permutation *p = gsl_permutation_alloc(size);
    int s;

    // Compute the LU decomposition of this matrix
    gsl_linalg_LU_decomp(matrix, p, &s);

    // Compute the  inverse of the LU decomposition
    gsl_matrix *inv = gsl_matrix_alloc(size, size);
    gsl_linalg_LU_invert(matrix, p, inv);

    gsl_permutation_free(p);

    return inv;
}

struct slice_t
{
  Float_t wmin;
  Float_t wmax;
  Float_t wmean;
  Float_t werr;
  long long int entries;
  std::vector<int> tChannel;
  std::vector<Float_t> averageDelay;
  std::vector<Float_t> averageDeltaT;
  std::vector<Float_t> varianceDeltaT;
  std::vector<long long int> nVarianceDeltaT;
  std::vector<long long int> nDeltaT;
  Float_t **covariance;
  Float_t **inverse_covariance;
  long long int **entries_covariance;
  Float_t **normalized_covariance;

  // TGraphErrors* delay;
  // TF1* delay_line;
  // TH2F* deltaTscatter;
  // TGraphErrors* deltaTgraph;
  // TF1* deltaTline;
  // std::vector<TH1F*> deltaTslice;
  // std::vector<double> wmean;
  // std::vector<double> werr;
  //
  //
  // std::vector<double> dx;
  // std::vector<double> dy;
  // std::vector<double> dex;
  // std::vector<double> dey;
};


struct Crystal_t
{
  int number;
  int detectorChannel;
  int timingChannel;
  std::vector<int> relevantForW;
  std::vector<int> delayTimingChannels;
  TCut *CrystalCut;
  TCut *CrystalCutWithoutCutG;
  TCut *PhotopeakEnergyCut;
  std::vector<TCutG*> cutg;
  TGraph *calibrationGraph;
  TGraph *wz;
  TH1F *simpleCTR;
  TH1F *centralCTR;
  TH1F *allCTR;
  TH1F *poliCorrCTR;
  TH1F *likeCTR;
  TH1F *hybridCTR;
  std::vector<double> vSimple;
  std::vector<double> vCentral;
  std::vector<double> vAll;
  std::vector<double> vPoli;
  std::vector<double> vLike;
  std::vector<double> vhybrid;
  TH1F *simpleCTR_norm;
  TH1F *centralCTR_norm;
  TH1F *allCTR_norm;
  TH1F *poliCorrCTR_norm;
  TH1F *hybridCTR_norm;
  TTreeFormula *Formula;
  TH1F *likeCTR_norm;

  std::vector<double> z;
  TGraphErrors* tw_correction;
  TGraphErrors* rms_tw_correction;
  std::vector<TGraphErrors*> delay;
  std::vector<TGraphErrors*> rms_delay;
  TF1 *tw_correction_line;
  TF1 *rms_tw_correction_line;
  std::vector<TF1*> delay_line;
  std::vector<TF1*> rms_delay_line;
  const char* path;
  bool accepted;
  bool polishedCorrection;
  std::vector<int>    tChannelsForPolishedCorrection;
  std::vector<double> meanForPolishedCorrection;
  std::vector<double> fwhmForPolishedCorrection;

  std::vector<slice_t> slice;

  TGraph *** inverse_covariance_element; // matrix of TGraph, one for each element of the inverse covariance element s^{-1}_i,j(w) that is a function of doi...
  TF1 *** inverse_covariance_element_line;


  // Float_t *variance;
  // long long int *entries_variance;
  // std::vector<Float_t**> covariance;
  // std::vector<long long int**> entries_covariance;
  // std::vector<Float_t**> inverse_covariance;
  //
  // TH2F** deltaTscatter;

  TH1F* lightCentralHisto;
  TH1F* lightAllHisto;
  TH1F* basicCTRhisto;

  // TCanvas
};

struct detector_t
{
  int digitizerChannel;
  float saturation;
  float pedestal;
};




void extractFromHisto(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double* res)
{
  // TF1 *gexp;
  // TF1 *cb;

  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // resctrict the fitting range of gauss function

  gaussDummy->SetLineColor(kRed);
  double fitGaussMin = histo->GetMean()-2.0*histo->GetRMS();
  double fitGaussMax = histo->GetMean()+2.0*histo->GetRMS();
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  if(fitGaussMin < f1min)
  {
    fitGaussMin = f1min;
  }
  if(fitGaussMax > f1max)
  {
    fitGaussMax = f1max;
  }
  TFitResultPtr rGauss = histo->Fit(gaussDummy,"QNS","",fitGaussMin,fitGaussMax);
  Int_t fitStatusGauss= rGauss;

  //NB fit results converted to int gives the fit status:
  // fitStatusGauss == 0 -> fit OK
  // fitStatusGauss != 0 -> fit FAILED

  double fitMin;
  double fitMax;
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    fitMin = fitGaussMin;
    fitMax = fitGaussMax;
  }
  else
  {
    // use fit values
    fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
    fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  }

  // chech that they are not outside the limits defined by user
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  //fit with crystalball
  TF1 *cb  = new TF1("cb","crystalball",f1min,f1max);
  cb->SetLineColor(kBlue);
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    cb->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS(),1,3);
  }
  else
  {
    // use fit values
    cb->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  }
  TFitResultPtr rCb = histo->Fit(cb,"QNS","",fitMin,fitMax);

  //fit with gauss + exp
  TF1* gexp  = new TF1("gexp","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"Mean");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    gexp->SetParameter(0,histo->GetEntries());
    gexp->SetParameter(1,histo->GetMean());
    gexp->SetParameter(2,histo->GetRMS());
    gexp->SetParameter(3,histo->GetRMS()); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  else
  {
    // use fit values
    gexp->SetParameter(0,gaussDummy->GetParameter(0));
    gexp->SetParameter(1,gaussDummy->GetParameter(1));
    gexp->SetParameter(2,gaussDummy->GetParameter(2));
    gexp->SetParameter(3,gaussDummy->GetParameter(2)); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  TFitResultPtr rGexp = histo->Fit(gexp,"QNS","",fitMin,fitMax);

  Int_t fitStatusCb = rCb;
  Int_t fitStatusGexp = rGexp;

  double chi2gexp;
  double chi2cb;

  if(fitStatusGexp == 0) // if Gexp worked
  {
    chi2gexp = rGexp->Chi2();
  }
  if(fitStatusCb == 0)// if cb worked
  {
    chi2cb   = rCb->Chi2();
  }


  //set function to measure ctr etc...
  TF1 *f1;
  // std::cout << histo->GetName() << std::endl;
  // std::cout << fitStatusGexp << " "
  //           << fitStatusCb << " "
  //           << fitStatusGauss << " "
  //           << std::endl;
  if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss != 0)) // all fit didn't work, just set everything to 0
  {
    // std::cout << "None" << std::endl;
    res[0] = 0;
    res[1] = 0;
    res[2] = 0;
    res[3] = 0;
  }
  else
  {
    if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss == 0)) // only gauss worked
    {
      // std::cout << "Gauss" << std::endl;
      f1 = gaussDummy;
      histo->Fit(f1,"Q","",fitGaussMin,fitGaussMax);
      res[0] = f1->GetParameter(1);
      res[1] = f1->GetParameter(2);
      res[2] = f1->GetParError(1);
      res[3] = f1->GetParError(2);
      delete gexp;
      delete cb;
    }
    else // one between gexp and cb worked
    {
      if((fitStatusGexp  != 0) && (fitStatusCb == 0)) // only cb worked
      {
        // std::cout << "cb" << std::endl;
        f1 = cb;
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete gexp;

      }
      else if((fitStatusGexp  == 0) && (fitStatusCb != 0)) // only gexp worked
      {
        // std::cout << "gexp" << std::endl;
        f1 = gexp;
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete cb;
      }
      else // both worked
      {
        if(chi2gexp > chi2cb)
        {
          // std::cout << "cb better than gexp" << std::endl;
          f1 = cb;
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete gexp;
        }
        else
        {
          // std::cout << "gexp better than cb" << std::endl;
          f1 = gexp;
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete cb;
        }
      }
      // f1->SetLineColor(kRed);
      res[0] = f1->GetParameter(1);
      res[1] = f1->GetParameter(2);
      res[2] = f1->GetParError(1);
      res[3] = f1->GetParError(2);
    }
  }
  delete cTemp;
}

/*** find width at half max ***/
// float ComputeFWHM(TH1F* histo)
// {
//    float max = histo->GetMaximum();
//    float halfMax = max / 2.0;
//    int binMin = histo->FindFirstBinAbove(halfMax);
//    int binMax = histo->FindLastBinAbove(halfMax);
//    float down = histo->GetBinCenter(binMin);
//    float up   = histo->GetBinCenter(binMax);
//    float ret = up -down;
//    return ret;
// }

void extractWithGaussAndExp(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double tagFwhm, double* res)
{
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  histo->Fit(gaussDummy,"QN");

  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();


  TF1* f1  = new TF1("f1","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))");
  f1->SetLineColor(kBlack);
  f1->SetParName(0,"N");
  f1->SetParName(1,"Mean");
  f1->SetParName(2,"Sigma");
  f1->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  f1->SetParameter(0,gaussDummy->GetParameter(0));
  f1->SetParameter(1,gaussDummy->GetParameter(1));
  f1->SetParameter(2,gaussDummy->GetParameter(2));
  f1->SetParameter(3,5e-9);
  double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  histo->Fit(f1,"","",fitMin,fitMax);
  std::cout << "-----------------------------------------------------------------------" << std::endl;
  std::cout << histo->GetName() << std::endl;
  std::cout << f1->GetMaximum(fitMin,fitMax)    << "\t"
            // << fitMin              << "\t"
            // << fitMax              << "\t"
            << f1->GetParameter(0) << "\t"
            << f1->GetParameter(1) << "\t"
            << f1->GetParameter(2) << "\t"
            << f1->GetParameter(3) << "\t"
            << std::endl;
  double min,max,min10,max10;
  // int divs = 3000;
  double step = (fitMin-fitMax)/divs;
  double funcMax = f1->GetMaximum(fitMin,fitMax);
  for(int i = 0 ; i < divs ; i++)
  {
    if( (f1->Eval(f1min + i*step) < funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/2.0) )
    {
      min = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/2.0) )
    {
      max = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) < funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/10.0) )
    {
      min10 = f1min + (i+0.5)*step;
    }
    if( (f1->Eval(f1min + i*step) > funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/10.0) )
    {
      max10 = f1min + (i+0.5)*step;
    }
  }
  // res[0] = f1->GetParameter(1);  // res[0] is mean
  // res[1] = max-min;              // res[1] is FWHM
  res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
  res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));
  std::cout <<"----> " << res[0] << " " << res[1]<< std::endl;
  delete cTemp;

}

void extractCTR(TH1F* histo,double fitPercMin,double fitPercMax, int divs, double tagFwhm, double* res, double* fitRes)
{
  //first, dummy gaussian fit
  // TCanvas *cTemp  = new TCanvas("temp","temp");
  // TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // histo->Fit(gaussDummy,"QN");
  //
  // double f1min = histo->GetXaxis()->GetXmin();
  // double f1max = histo->GetXaxis()->GetXmax();
  // // std::cout << f1min << " " << f1max << std::endl;
  // TF1* f1  = new TF1("f1","crystalball");
  // f1->SetLineColor(kBlack);
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  // double fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
  // double fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  // if(fitMin < f1min)
  // {
  //   fitMin = f1min;
  // }
  // if(fitMax > f1max)
  // {
  //   fitMax = f1max;
  // }
  // histo->Fit(f1,"Q","",fitMin,fitMax);


  // preliminary gauss fit
  TCanvas *cTemp  = new TCanvas("temp","temp");
  TF1 *gaussDummy = new TF1("gaussDummy","gaus");
  // resctrict the fitting range of gauss function

  gaussDummy->SetLineColor(kRed);
  double fitGaussMin = histo->GetMean()-2.0*histo->GetRMS();
  double fitGaussMax = histo->GetMean()+2.0*histo->GetRMS();
  double f1min = histo->GetXaxis()->GetXmin();
  double f1max = histo->GetXaxis()->GetXmax();
  if(fitGaussMin < f1min)
  {
    fitGaussMin = f1min;
  }
  if(fitGaussMax > f1max)
  {
    fitGaussMax = f1max;
  }
  TFitResultPtr rGauss = histo->Fit(gaussDummy,"QNS","",fitGaussMin,fitGaussMax);
  Int_t fitStatusGauss= rGauss;

  //NB fit results converted to int gives the fit status:
  // fitStatusGauss == 0 -> fit OK
  // fitStatusGauss != 0 -> fit FAILED

  double fitMin;
  double fitMax;
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    fitMin = fitGaussMin;
    fitMax = fitGaussMax;
  }
  else
  {
    // use fit values
    fitMin = gaussDummy->GetParameter(1) - fitPercMin*(gaussDummy->GetParameter(2));
    fitMax = gaussDummy->GetParameter(1) + fitPercMax*(gaussDummy->GetParameter(2));
  }

  // chech that they are not outside the limits defined by user
  if(fitMin < f1min)
  {
    fitMin = f1min;
  }
  if(fitMax > f1max)
  {
    fitMax = f1max;
  }

  //fit with crystalball
  TF1 *cb  = new TF1("cb","crystalball",f1min,f1max);
  cb->SetLineColor(kBlue);
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    cb->SetParameters(histo->GetEntries(),histo->GetMean(),histo->GetRMS(),1,3);
  }
  else
  {
    // use fit values
    cb->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  }
  TFitResultPtr rCb = histo->Fit(cb,"QNS","",fitMin,fitMax);

  //fit with gauss + exp
  TF1* gexp  = new TF1("gexp","[0]/sqrt(2)*exp([2]^2/2/[3]^2-(x-[1])/[3])*(1-TMath::Erf(([1]-x+[2]^2/[3])/(sqrt(2*[2]^2))))",f1min,f1max);
  gexp->SetLineColor(kGreen);
  gexp->SetParName(0,"N");
  gexp->SetParName(1,"Mean");
  gexp->SetParName(2,"Sigma");
  gexp->SetParName(3,"tau");
  // f1->SetParameters(gaussDummy->GetParameter(0),gaussDummy->GetParameter(1),gaussDummy->GetParameter(2),1,3);
  if(fitStatusGauss != 0) // gauss fit didn't work
  {
    // use the histogram values
    gexp->SetParameter(0,histo->GetEntries());
    gexp->SetParameter(1,histo->GetMean());
    gexp->SetParameter(2,histo->GetRMS());
    gexp->SetParameter(3,histo->GetRMS()); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  else
  {
    // use fit values
    gexp->SetParameter(0,gaussDummy->GetParameter(0));
    gexp->SetParameter(1,gaussDummy->GetParameter(1));
    gexp->SetParameter(2,gaussDummy->GetParameter(2));
    gexp->SetParameter(3,gaussDummy->GetParameter(2)); // ROOT really needs all parameters initialized, and a "good" guess for tau is the sigma of the previous fit...
  }
  TFitResultPtr rGexp = histo->Fit(gexp,"QNS","",fitMin,fitMax);

  Int_t fitStatusCb = rCb;
  Int_t fitStatusGexp = rGexp;

  double chi2gexp;
  double chi2cb;

  if(fitStatusGexp == 0) // if Gexp worked
  {
    chi2gexp = rGexp->Chi2();
  }
  if(fitStatusCb == 0)// if cb worked
  {
    chi2cb   = rCb->Chi2();
  }


  //set function to measure ctr etc...
  TF1 *f1;
  if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss != 0)) // all fit didn't work, just set everything to 0
  {
    res[0] = 0;
    res[1] = 0;

    fitRes[0] = 0;
    fitRes[1] = 0;
    fitRes[2] = 0;
    // res[2] = 0;
    // res[3] = 0;
  }
  else
  {
    if((fitStatusGexp  != 0) && (fitStatusCb != 0) && (fitStatusGauss == 0)) // only gauss worked
    {
      f1 = gaussDummy;
      f1->SetLineColor(kRed);
      histo->Fit(f1,"Q","",fitGaussMin,fitGaussMax);
      res[0] = sqrt(2)*sqrt(pow((2.355*f1->GetParameter(2)),2)-pow(tagFwhm,2));
      res[1] = sqrt(2)*sqrt(pow((4.29*f1->GetParameter(2)),2)-pow((tagFwhm/2.355)*4.29,2));

      fitRes[0] = f1->GetChisquare();
      fitRes[1] = f1->GetNDF();
      fitRes[2] = f1->GetProb();

      delete gexp;
      delete cb;
    }
    else
    {
      if((fitStatusGexp  != 0) && (fitStatusCb == 0)) // only cb worked
      {
        f1 = cb;
        f1->SetLineColor(kRed);
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete gexp;

      }
      else if((fitStatusGexp  == 0) && (fitStatusCb != 0)) // only gexp worked
      {
        f1 = gexp;
        f1->SetLineColor(kRed);
        histo->Fit(f1,"Q","",fitMin,fitMax);
        delete cb;
      }
      else // both worked
      {
        if(chi2gexp > chi2cb)
        {
          f1 = cb;
          f1->SetLineColor(kRed);
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete gexp;
        }
        else
        {
          f1 = gexp;
          f1->SetLineColor(kRed);
          histo->Fit(f1,"Q","",fitMin,fitMax);
          delete cb;
        }
      }

      double min,max,min10,max10;
      // int divs = 3000;
      double step = (f1max-f1min)/divs;
      double funcMax = f1->GetMaximum(fitMin,fitMax);
      for(int i = 0 ; i < divs ; i++)
      {
        if( (f1->Eval(f1min + i*step) < funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/2.0) )
        {
          min = f1min + (i+0.5)*step;
        }
        if( (f1->Eval(f1min + i*step) > funcMax/2.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/2.0) )
        {
          max = f1min + (i+0.5)*step;
        }
        if( (f1->Eval(f1min + i*step) < funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) > funcMax/10.0) )
        {
          min10 = f1min + (i+0.5)*step;
        }
        if( (f1->Eval(f1min + i*step) > funcMax/10.0) && (f1->Eval(f1min + (i+1)*step) < funcMax/10.0) )
        {
          max10 = f1min + (i+0.5)*step;
        }
      }
      res[0] = sqrt(2)*sqrt(pow((max-min),2)-pow(tagFwhm,2));
      res[1] = sqrt(2)*sqrt(pow((max10-min10),2)-pow((tagFwhm/2.355)*4.29,2));

      fitRes[0] = f1->GetChisquare();
      fitRes[1] = f1->GetNDF();
      fitRes[2] = f1->GetProb();
      // std::cout << f1->GetChisquare()/f1->GetNDF() << std::endl;
      delete cTemp;
    }
  }
}


//**** per std::vector -- non binnata
double FindSmallestInterval(double& mean,
                            double& meanErr,
                            double& min,
                            double& max,
                            std::vector<double>& vals,
                            const double& fraction,
                            const bool& verbosity)
{
   if( verbosity )
     std::cout << ">>>>>> FindSmallestInterval" << std::endl;


   std::sort(vals.begin(),vals.end());

   unsigned int nPoints = vals.size();
   unsigned int maxPoints = (unsigned int)(fraction * nPoints);

   unsigned int minPoint = 0;
   unsigned int maxPoint = 0;
   double delta = 999999.;
   for(unsigned int point = 0; point < nPoints-maxPoints; ++point)
   {
     double tmpMin = vals.at(point);
     double tmpMax = vals.at(point+maxPoints-1);
     if( tmpMax-tmpMin < delta )
     {
       delta = tmpMax - tmpMin;
       min = tmpMin;
       max = tmpMax;
       minPoint = point;
       maxPoint = point + maxPoints - 1;
     }
   }
   return delta;
}


/*** find effective sigma ***/
void FindSmallestInterval(double* retValues, TH1F* histo, const float&
fraction, const bool& verbosity, double tagFwhm)
{
  float ret[4];
   float integralMax = fraction * histo->Integral();

   int N = histo -> GetNbinsX();
   std::vector<float> binCenters(N);
   std::vector<float> binContents(N);
   std::vector<float> binIntegrals(N);
   for(int bin1 = 0; bin1 < N; ++bin1)
   {
     binCenters[bin1] = histo->GetBinCenter(bin1+1);
     binContents[bin1] = histo->GetBinContent(bin1+1);

     for(int bin2 = 0; bin2 <= bin1; ++bin2)
       binIntegrals[bin1] += binContents[bin2];
   }

   float min = 0.;
   float max = 0.;
   float delta = 999999.;
   for(int bin1 = 0; bin1 < N; ++bin1)
   {
     for(int bin2 = bin1+1; bin2 < N; ++bin2)
     {
       if( (binIntegrals[bin2]-binIntegrals[bin1]) < integralMax ) continue;

       float tmpMin = histo -> GetBinCenter(bin1);
       float tmpMax = histo -> GetBinCenter(bin2);

       if( (tmpMax-tmpMin) < delta )
       {
         delta = (tmpMax - tmpMin);
         min = tmpMin;
         max = tmpMax;
       }

       break;
     }
   }

   TH1F* smallHisto = (TH1F*)( histo->Clone("smallHisto") );
   for(int bin = 1; bin <= smallHisto->GetNbinsX(); ++bin)
   {
     if( smallHisto->GetBinCenter(bin) < min )
       smallHisto -> SetBinContent(bin,0);

     if( smallHisto->GetBinCenter(bin) > max )
       smallHisto -> SetBinContent(bin,0);
   }
   smallHisto -> SetFillColor(kYellow);

   float mean = smallHisto -> GetMean();
   float meanErr = smallHisto -> GetMeanError();

   ret[0] = mean;
   ret[1] = meanErr;
   ret[2] = min;
   ret[3] = max;

   //mean is the smallest interval containing the 68% (fraction) of data. this would be from -1 sigma to +1 sigma, so 2 sigmas. therefore we get teh "fwhm" of this distro by
   double fwhm = 2.355* ((max-min) / 2.0);

   retValues[0] = sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
   retValues[1] = 0;
}


bool compareByNumber(const Crystal_t &a,const Crystal_t  &b)
{
  return a.number < b.number;
}

void usage()
{
  std::cout << "\t\t" << "[-i|--input] <file_prefix>  [-o|--output] <output.root> [-c|--calibration] calibration.root [--coincidence] coincidence.root [OPTIONS]" << std::endl
            << "\t\t" << "<file_prefix>                                      - prefix of TTree files to analyze"   << std::endl
            << "\t\t" << "<output.root>                                      - output file name"   << std::endl
            << "\t\t" << "<calibration.root>                                 - calibration file " << std::endl
            // << "\t\t" << "<coincidence.root>                                 - time calibration file " << std::endl
            << "\t\t" << "--simulation                                       - the datast is from a simulation (therefore the tagging photopeak is ignored)" << std::endl
            << "\t\t" << "--length <value>                                   - crystal length in mm, default = 15.0"  << std::endl
            << "\t\t" << "--doiFraction <value>                              - fraction of DOI length towards which the time stamps are corrected (from 0 to 1)"  << std::endl
            << "\t\t" << "                                                   - 0 = front of the crystal (DOI close to detector) "  << std::endl
            << "\t\t" << "                                                   - 1 = back of the crystal (DOI far from detector) "  << std::endl
            << "\t\t" << "--tagFwhm <value>                                  - FWHM timing resolution of reference board, in sec - default = 70e-12"  << std::endl
            << "\t\t" << "--rmsLow <value>                                   - lower bound of CTR fit -> mean - rmsLow*mean - default = 1.75"  << std::endl
            << "\t\t" << "--rmsHigh <value>                                  - upper bound of CTR fit -> mean + rmsHigh*mean - default = 1.75"  << std::endl
            << "\t\t" << "--histoMin <value>                                 - lower limit of CTR spectra, in sec - default = -15e-9"  << std::endl
            << "\t\t" << "--histoMax <value>                                 - upper limit of CTR spectra, in sec - default = 15e-9"  << std::endl
            << "\t\t" << "--histoBins <value>                                - n of bins for CTR spectra - default = 500"  << std::endl
            << "\t\t" << "--smooth <value>                                   - n of iteration in CTR histograms smoothing - default = 0 (no smoothing)"  << std::endl
            << "\t\t" << "--fitPercMin <value>                               - time fit min is set to ((gauss fit mean) - fitPercMin*(gauss fit sigma))  - default = 5"  << std::endl
            << "\t\t" << "--fitPercMax <value>                               - time fit max is set to ((gauus fit mean) - fitPercMax*(gauss fit sigma))  - default = 6" << std::endl
            << "\t\t" << "--divs <value>                                     - n of divisions when looking for FWHM - default = 10000"  << std::endl
            << "\t\t" << "--bins <value>                                     - n of bins in summary CTR histograms - deafult 40"  << std::endl
            << "\t\t" << "--func <value>                                     - function for fitting (default = 0)"  << std::endl
            << "\t\t" << "                                                   - 0 = crystalball "  << std::endl
            << "\t\t" << "                                                   - 1 = gauss+exp "  << std::endl
            << "\t\t" << "--unbinned                                         - use also the unbinned method to calculate CTR - default = 0 (false)"  << std::endl
            << "\t\t" << "--fitCorrection                                    - use line fit to perform correction   - default = not given (false)"  << std::endl
            << "\t\t" << "--exclude-channels                                 - channels to exclude from time correction, comma separated - default = "" "  << std::endl
            << "\t\t" << "--start-time                                       - acq time from which events are accepted [h]  - default = 0"  << std::endl
            << "\t\t" << "--sliced                                           - if given, it's a slice acq                   - default = not given"  << std::endl
            << "\t\t" << "--likelihood                                       - if given, perform likelihood correction                   - default = not given"  << std::endl
            << "\t\t" << "--likeMin <value>                                  - lower limit of likelihood spectra, in sec - default = -5e-9"  << std::endl
            << "\t\t" << "--likeMax <value>                                  - upper limit of likelihood spectra, in sec - default = 5e-9"  << std::endl
            << "\t\t" << "--likeBins <value>                                 - n of bins for likelihood spectra - default = 500"  << std::endl
            << "\t\t" << "--wBins <value>                                    - number of bins in w slicing for likelihood - default = 10"  << std::endl
            << "\t\t" << "--margin <value>                                   - margin in w slicing [mm] - default = 0.1"  << std::endl
            << "\t\t" << "--basicLikelihood                                  - likelihood without line fits   - default = not given (false)"  << std::endl
            << "\t\t" << "--likelihoodLine                                   - using line fit of inverse covariance matrix, instead of tgraph eval. valid only if basicLikelihood is false - default = not given (false)"  << std::endl
            << "\t\t" << "--hybridCorrection                                 - performing hybrid correction - default = not given (false)"  << std::endl
            << "\t\t" << std::endl;
}

int main (int argc, char** argv)
{
  if(argc < 2)
  {
    std::cout << argv[0] << std::endl;
    usage();
    return 1;
  }

  std::stringstream streamCommand;
  for(int i=0 ; i < argc; i++)
  {
    streamCommand << argv[i] << " ";
  }

  gStyle->SetOptStat(1111);
  gStyle->SetOptFit(1111);

  std::string inputFileName = "";
  std::string outputFileName = "";
  std::string calibrationFileName = "";
  // std::string coincidenceCalibrationFileName = "";

  std::string exclude_channels = "";
  bool exclude = false;
  bool simulation = false;
  Float_t length = 15.0; //mm
  Float_t doiFraction = 0.5;
  Float_t tagFwhm = 70.0e-12; //s
  Float_t rmsLow = 1.75;
  Float_t rmsHigh = 1.75;
  Float_t histoMin = -15e-9;//s
  Float_t histoMax = 15e-9;//s
  Float_t fitPercMin = 5;
  Float_t fitPercMax = 6;
  Float_t likeMin = -15e-9;//s
  Float_t likeMax = 15e-9;//s
  int likeBins = 500;
  int divs       = 10000;
  int histoBins = 500;
  int smooth = 0; //
  int bins = 40;
  double minCTR = 100;
  double maxCTR = 500;
  int func = 0;
  bool unbinned = false;
  bool fitCorrection = false;
  bool basicLikelihood = false;
  bool hybridCorrection = false;
  double start_time = 0;
  bool sliced = false;
  bool likelihood = false;
  bool likelihoodLine = false;
  int WrangeBinsForTiming = 10;
  float marginWZgraph = 0.1;

  // parse arguments
  static struct option longOptions[] =
  {
			{ "input", required_argument, 0, 0 },
      { "output", required_argument, 0, 0 },
      { "calibration", required_argument, 0, 0 },
      { "simulation", no_argument, 0, 0 },
      { "length", required_argument, 0, 0 },
      { "doiFraction", required_argument, 0, 0 },
      // { "coincidence", required_argument, 0, 0 },
      { "tagFwhm", required_argument, 0, 0 },
      { "rmsLow", required_argument, 0, 0 },
      { "rmsHigh", required_argument, 0, 0 },
      { "histoMin", required_argument, 0, 0 },
      { "histoMax", required_argument, 0, 0 },
      { "histoBins", required_argument, 0, 0 },
      { "smooth", required_argument, 0, 0 },
      { "fitPercMin", required_argument, 0, 0 },
      { "fitPercMax", required_argument, 0, 0 },
      { "divs", required_argument, 0, 0 },
      { "bins", required_argument, 0, 0 },
      { "func", required_argument, 0, 0 },
      { "unbinned", no_argument, 0, 0 },
      { "fitCorrection", no_argument, 0, 0 },
      { "exclude-channels", required_argument, 0, 0 },
      { "start-time", required_argument, 0, 0 },
      { "sliced", no_argument, 0, 0 },
      { "likelihood", no_argument, 0, 0 },
      { "likeMin", required_argument, 0, 0 },
      { "likeMax", required_argument, 0, 0 },
      { "likeBins", required_argument, 0, 0 },
      { "wBins", required_argument, 0, 0 },
      { "margin", required_argument, 0, 0 },
      { "basicLikelihood", no_argument, 0, 0 },
      { "likelihoodLine", no_argument, 0, 0 },
      { "hybridCorrection", no_argument, 0, 0 },
			{ NULL, 0, 0, 0 }
	};

  while(1) {
		int optionIndex = 0;
		int c = getopt_long(argc, argv, "i:o:c:", longOptions, &optionIndex);
		if (c == -1) {
			break;
		}
		if (c == 'i'){
			inputFileName = (char *)optarg;
    }
		else if (c == 'o'){
      outputFileName = (char *)optarg;
    }
    else if (c == 'c'){
      calibrationFileName = (char *)optarg;
    }
		else if (c == 0 && optionIndex == 0){
      inputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 1){
      outputFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 2){
      calibrationFileName = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 3){
      std::cout << "Dataset from simulation " << std::endl;
      simulation = true;
    }
    else if (c == 0 && optionIndex == 4){
      length = atof((char *)optarg);;
    }
    else if (c == 0 && optionIndex == 5){
      doiFraction = atof((char *)optarg);;
    }
    // else if (c == 0 && optionIndex == 6){
    //   coincidenceCalibrationFileName = (char *)optarg;
    // }
    else if (c == 0 && optionIndex == 6){
      tagFwhm = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 7){
      rmsLow = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 8){
      rmsHigh = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 9){
      histoMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 10){
      histoMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 11){
      histoBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 12){
      smooth = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 13){
      fitPercMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 14){
      fitPercMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 15){
      divs = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 16){
      bins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 17){
      func = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 18){
      unbinned = true;
    }
    else if (c == 0 && optionIndex == 19){
      fitCorrection = true;
    }
    else if (c == 0 && optionIndex == 20){
      exclude = true;
      exclude_channels = (char *)optarg;
    }
    else if (c == 0 && optionIndex == 21){
      start_time = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 22){
      sliced = true;
    }
    else if (c == 0 && optionIndex == 23){
      likelihood = true;
    }
    else if (c == 0 && optionIndex == 24){
      likeMin = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 25){
      likeMax = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 26){
      likeBins = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 27){
      WrangeBinsForTiming = atoi((char *)optarg);
    }
    else if (c == 0 && optionIndex == 28){
      marginWZgraph = atof((char *)optarg);
    }
    else if (c == 0 && optionIndex == 29){
      basicLikelihood = true;
    }
    else if (c == 0 && optionIndex == 30){
      likelihoodLine = true;
    }
    else if (c == 0 && optionIndex == 31){
      hybridCorrection = true;
    }
		else {
      std::cout	<< "Usage: " << argv[0] << std::endl;
			usage();
			return 1;
		}
	}


  std::vector<int> forbidden_channels;
  if(exclude)
  {
    // std::vector<int> vect;
    std::stringstream ss(exclude_channels);
    int i;
    while (ss >> i)
    {
      forbidden_channels.push_back(i);
      if (ss.peek() == ',')
        ss.ignore();
    }
    std::cout << "Channels excluded from time correction (for depolished): " << std::endl;
    for (i=0; i< forbidden_channels.size(); i++)
        std::cout << forbidden_channels.at(i)<<std::endl;
  }

  std::cout << "Chosen (length * doiFraction) = " << length * doiFraction << std::endl;
  if(fitCorrection)
  {
    std::cout << "Using linear fits to perform time correction" << std::endl;
  }
  if(likelihood)
  {
    std::cout << "Performing likelihood correction " << std::endl;
  }
  if(basicLikelihood)
  {
    std::cout << "Likelihood correction uses sliced arrays" << std::endl;
  }
  else
  {
    if(likelihoodLine)
    {
      std::cout << "Likelihood correction interpolates arrays with lines" << std::endl;
    }
    else
    {
      std::cout << "Likelihood correction interpolates arrays with TGraphs" << std::endl;
    }
  }
  if(hybridCorrection)
  {
    std::cout << "Performing hybrid correction" << std::endl;
  }



  // read file in dir
  std::vector<std::string> v;
  read_directory(".", v);
  // std::copy(v.begin(), v.end(),std::ostream_iterator<std::string>(std::cout, "\n"));
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;

  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(v[i]);
    }
  }


  //prepare output text file
  std::string textFileName = outputFileName.substr(0,outputFileName.size()-5);
  textFileName += ".txt";
  // std::cout << textFileName << std::endl;

  std::ofstream textfile;
  textfile.open (textFileName.c_str(),std::ofstream::out);




  // TFile *treeFile = new TFile(inputFileName.c_str());
  // TTree* tree = (TTree*) treeFile->Get("adc");

  //----------------------------------------------------------//
  //  Get TChain of input files                               //
  //----------------------------------------------------------//
  TChain* tree = new TChain("adc");  // create the input tchain and the analysis ttree
  for(unsigned int i = 0 ; i < listInputFiles.size(); i++)
  {
    std::cout << "Adding file " << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }

  // Add several TTreeFormula to the list;
  TList formulas;

  //variables for the input TChain
  // ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  // ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous event
  // Int_t        *ChainAdcChannel;
  // Short_t      *ChainDesktopAdcChannel;                              // input TChain data for desktop digitizers - data is int_16
  // UShort_t     *ChainVMEadcChannel;                                  // input TChain data for VME digitizers - data is uint_16
  // Float_t      *ChainTimeStamp;
  // Float_t      *TDCBinning;
  // // Short_t      *ChainPetirocChannel;                                 //FIXME temporary data type of petiroc charge input - ask
  // Float_t       RealX;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       RealY;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       RealZ;                                               // "real" gamma interaction positions (from simulation data)
  // Float_t       simTaggingCharge;
  // Float_t       simTaggingTime;
  // Short_t       CrystalsHit;                                         // "real" number of crystals hit in the event (from simulation data)
  // Short_t       NumbOfInteractions;                                  // "real" number of interaction (energy depositions) in the event (from simulation data)
  //
  // //branches for the input TChain
  // TBranch      *bChainExtendedTimeTag;                               // branches for above data
  // TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  // TBranch     **bChainAdcChannel;                                    // branches for above data
  // TBranch     **bChainTimeStamp;
  // TBranch      *bRealX;                                              // branches for above data
  // TBranch      *bRealY;                                              // branches for above data
  // TBranch      *bRealZ;                                              // branches for above data
  // TBranch      *bsimTaggingCharge;                                              // branches for above data
  // TBranch      *bsimTaggingTime;                                              // branches for above data
  // TBranch      *bCrystalsHit;                                        // branches for above data
  // TBranch      *bNumbOfInteractions;                                 // branches for above data
  // TBranch      *bTotalCryEnergy;                                     //
  //
  // if(std::string(argv[1]) == std::string("-c")) // first argument is -c, then the config file name is passed by command line
  // {
  //   for (int i = 3; i < argc ; i++) // run on the remaining arguments to add all the input files
  //   {
  //     std::cout << "Adding file " << argv[i] << std::endl;
  //     tree->Add(argv[i]);
  //   }
  // }
  // else // the config file was indeed the default one
  // {
  //   for (int i = 1; i < argc ; i++) // run on the remaining arguments to add all the input files
  //   {
  //     std::cout << "Adding file " << argv[i] << std::endl;
  //     tree->Add(argv[i]);
  //   }
  // }

  std::vector<int> detector_channels;

  TObjArray *leavescopy = tree->GetListOfLeaves();
  int nLeaves = leavescopy->GetEntries();
  std::vector<std::string> leavesName;
  // fill a vector with the leaves names
  for(int i = 0 ; i < nLeaves ; i++)
  {
    leavesName.push_back(leavescopy->At(i)->GetName());
  }
  // count the entries that start with "ch"
  int numOfCh = 0;
  // int numOfCry = 0;
  std::string ch_prefix("ch");
  std::string t_prefix("t");

  // std::string cry_prefix("cry");
  for(int i = 0 ; i < nLeaves ; i++)
  {
    //     leavesName.push_back(leavescopy->At(i)->GetName());
    if (!leavesName[i].compare(0, ch_prefix.size(), ch_prefix))
    {
      numOfCh++;
      detector_channels.push_back(atoi( (leavesName[i].substr(ch_prefix.size(),leavesName[i].size()-ch_prefix.size())).c_str() )) ;
    }
    // if (!leavesName[i].compare(0, cry_prefix.size(), cry_prefix))
    // numOfCry++;
  }
  //the string "cry" appears 4 times per crystal..
  // numOfCry = numOfCry / 4;
  std::cout << "Detector Channels \t= " << numOfCh << std::endl;
  // std::cout << "Number of Crystals \t= "<< numOfCry << std::endl;


  // first, create the adc channels variables and branches
  // ChainAdcChannel        = new Int_t [numOfCh];
  // ChainDesktopAdcChannel = new Short_t [numOfCh]; // input from ADC desktop
  // ChainVMEadcChannel     = new UShort_t [numOfCh]; // input from VME
  // ChainTimeStamp         = new Float_t[numOfCh];
  // // TDCBinning             = new Float_t[numOfCh];
  // // DigitizerChannelOn     = new bool[adcChannels];
  // bChainAdcChannel       = new TBranch* [numOfCh];
  // bChainTimeStamp        = new TBranch* [numOfCh];
  ULong64_t     ChainExtendedTimeTag;                                // extended time tag
  ULong64_t     ChainDeltaTimeTag;                                   // delta tag from previous
  UShort_t      *charge;
  Float_t      *timeStamp;
  TBranch      *bChainExtendedTimeTag;                               // branches for above data
  TBranch      *bChainDeltaTimeTag;                                  // branches for above data
  TBranch      **bCharge;
  TBranch      **btimeStamp;

  charge = new UShort_t[numOfCh];
  timeStamp = new Float_t[numOfCh];
  bCharge = new TBranch*[numOfCh];
  btimeStamp = new TBranch*[numOfCh];

  // set branches for reading the input files
  tree->SetBranchAddress("ExtendedTimeTag", &ChainExtendedTimeTag, &bChainExtendedTimeTag);
  tree->SetBranchAddress("DeltaTimeTag", &ChainDeltaTimeTag, &bChainDeltaTimeTag);
  // if(usingRealSimData)
  // {
  //   tree->SetBranchAddress("RealX", &RealX, &bRealX);
  //   tree->SetBranchAddress("RealY", &RealY, &bRealY);
  //   tree->SetBranchAddress("RealZ", &RealZ, &bRealZ);
  //   // fchain->SetBranchAddress("Tagging", &simTaggingCharge, &bsimTaggingCharge);
  //   // fchain->SetBranchAddress("TaggingTimeStamp", &simTaggingTime, &bsimTaggingTime);
  //   tree->SetBranchAddress("CrystalsHit",&CrystalsHit, &bCrystalsHit);
  //   tree->SetBranchAddress("NumbOfInteractions",&NumbOfInteractions, &bNumbOfInteractions);
  //   // fchain->SetBranchAddress("TotalCryEnergy",&TotalCryEnergy, &bTotalCryEnergy);
  // }
  for (int i = 0 ; i < detector_channels.size() ; i++)
  {
    //empty the stringstreams
    std::stringstream sname;
    sname << "ch" << detector_channels[i];
    // stype << "ch" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&charge[detector_channels[i]],&bCharge[detector_channels[i]]);
    sname.str("");
    // stype.str("");

    sname << "t" << detector_channels[i];
    // stype << "t" << i << "/F";
    tree->SetBranchAddress(sname.str().c_str(),&timeStamp[detector_channels[i]],&btimeStamp[detector_channels[i]]);
    sname.str("");
    // stype.str("");
  }



  // find the crystals with complete calibration data
  // this means looking into TWO files:
  // 1. the standard calibration file, but produced with active timing part
  //    this is produced NOT in coincidence and will give all but two plot
  //    for the crystals calibration. the only two missing are the only two that HAVE
  //    to be measured in coincidence, i.e. the tw_correction plots (and it's RMS) and the tag photopeak cut
  // 2. a configuration run performed in coincidence, from which the tw_correction are taken
  //
  // afterwards the crystals will be accepted if they are present in both calibration files

  std::vector<Crystal_t> crystal;
  std::vector<detector_t> detectorSaturation;

  // STANDARD CALIBRATION FILE
  TFile *calibrationFile = new TFile(calibrationFileName.c_str());
  calibrationFile->cd("Module 0.0");
  TList *listModule = gDirectory->GetListOfKeys();
  int nKeysMod = listModule->GetEntries();
  std::vector<std::string> keysModName;
  // fill a vector with the leaves names
  std::string mppc_prefix("MPPC");
  for(int i = 0 ; i < nKeysMod ; i++){
    keysModName.push_back(listModule->At(i)->GetName());
  }

  TCut *taggingPhotopeakCut;
  int taggingCrystalTimingChannel;
  std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  std::string taggingCrystalTimingChannel_prefix("taggingCrystalTimingChannel");
  // std::string det_prefix("channels");
  // std::string saturation_prefix("saturation");

  std::vector<int> *pChannels;
  std::vector<float> *pSaturation;
  std::vector<float> *pPedestal;
  gDirectory->GetObject("channels",pChannels);
  gDirectory->GetObject("saturation",pSaturation);
  gDirectory->GetObject("pedestal",pPedestal);

  std::vector<int> DetChannels = pChannels[0];
  std::vector<float> saturation = pSaturation[0];
  std::vector<float> pedestal = pPedestal[0];

  for(unsigned int iSat = 0; iSat < DetChannels.size(); iSat++)
  {
    detector_t tempDetector;
    tempDetector.digitizerChannel = DetChannels[iSat];
    tempDetector.saturation = saturation[iSat];
    tempDetector.pedestal = pedestal[iSat];
    detectorSaturation.push_back(tempDetector);
  }

  std::vector<std::string> MPPCfolders;
  for(unsigned int i = 0 ; i < keysModName.size() ; i++)
  {
    if (!keysModName[i].compare(0, mppc_prefix.size(), mppc_prefix))
    {
      MPPCfolders.push_back(keysModName[i]);
    }
    if(!keysModName[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
    {
//       std::cout << keysCryName[i] << std::endl;
      taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
//       if(cut)
//         temp_crystal.CrystalCut = cut;
    }
    if(!keysModName[i].compare(0,taggingCrystalTimingChannel_prefix.size(),taggingCrystalTimingChannel_prefix)) // find tcut
    {
      std::stringstream snameCh;
      snameCh << ((TNamed*) gDirectory->Get(keysModName[i].c_str()))->GetTitle();
      taggingCrystalTimingChannel = atoi(snameCh.str().c_str());
    }

    // if(!keysModName[i].compare(0,det_prefix.size(),det_prefix)) // find tcut
    // {
    //   std::stringstream snameCh;
    //   std::vector<int> *v;
    //   gDirectory->GetObject("channelsNumRelevantForW",v);
    //   temp_crystal.relevantForW = v[0];
    // }
    //
    // if(!keysModName[i].compare(0,saturation_prefix.size(),saturation_prefix)) // find tcut
    // {
    //   std::stringstream snameCh;
    //   std::vector<int> *v;
    //   gDirectory->GetObject("channelsNumRelevantForW",v);
    //   temp_crystal.relevantForW = v[0];
    // }


  }

  std::stringstream sformulaname;
  sformulaname << "FormulaTag";
  TCut taggingPhotopeakCutName;
  taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // std::cout << "FormulaTag ------------- \n" << taggingPhotopeakCutName << std::endl;
  TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);
  formulas.Add(FormulaTag);



  for(unsigned int iMppc = 0 ; iMppc < MPPCfolders.size() ; iMppc++)
  {
    // std::cout << MPPCfolders[iMppc] << std::endl;
    gDirectory->cd(MPPCfolders[iMppc].c_str());
    TList *listMppc = gDirectory->GetListOfKeys();
    int nKeysMppc = listMppc->GetEntries();
    std::vector<std::string> keysMppcName;
    // fill a vector with the leaves names
    std::string crystal_prefix("Crystal");
    // std::string det_prefix("digitizerChannel");
    // std::string saturation_prefix("saturation");
    for(int i = 0 ; i < nKeysMppc ; i++){
      keysMppcName.push_back(listMppc->At(i)->GetName());
    }

    std::vector<std::string> CrystalFolders;

    for(unsigned int i = 0 ; i < keysMppcName.size() ; i++)
    {
      if (!keysMppcName[i].compare(0, crystal_prefix.size(), crystal_prefix))
      {
        CrystalFolders.push_back(keysMppcName[i]);
      }
      // if (!keysMppcName[i].compare(0, det_prefix.size(), det_prefix))
      // {
      //   std::stringstream snameCh;
      //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
      //   tempDetector.digitizerChannel = atoi(snameCh.str().c_str());
      // }
      // if (!keysMppcName[i].compare(0, saturation_prefix.size(), saturation_prefix))
      // {
      //   std::stringstream snameCh;
      //   snameCh << ((TNamed*) gDirectory->Get(keysMppcName[i].c_str()))->GetTitle();
      //   tempDetector.saturation = atof(snameCh.str().c_str());
      // }
    }
    // detectorSaturation.push_back(tempDetector);

    for(unsigned int iCry = 0 ; iCry < CrystalFolders.size() ; iCry++)
    {
      //  std::cout << CrystalFolders[iCry] << std::endl;
       gDirectory->cd(CrystalFolders[iCry].c_str());

       Crystal_t temp_crystal;
       temp_crystal.number = -1;
       temp_crystal.CrystalCut = NULL;
       temp_crystal.CrystalCutWithoutCutG = NULL;
       temp_crystal.PhotopeakEnergyCut = NULL;
       temp_crystal.calibrationGraph = NULL;
       temp_crystal.simpleCTR = NULL;
       temp_crystal.centralCTR = NULL;
       temp_crystal.allCTR = NULL;
       temp_crystal.poliCorrCTR = NULL;
       temp_crystal.simpleCTR_norm = NULL;
       temp_crystal.centralCTR_norm = NULL;
       temp_crystal.allCTR_norm = NULL;
       temp_crystal.poliCorrCTR_norm = NULL;
       temp_crystal.likeCTR = NULL;
       temp_crystal.likeCTR_norm = NULL;
       temp_crystal.hybridCTR = NULL;
       temp_crystal.hybridCTR_norm = NULL;
       temp_crystal.wz = NULL;
       temp_crystal.accepted = true;
       temp_crystal.tw_correction = NULL;
       temp_crystal.polishedCorrection = false;


       //get crystal number
       temp_crystal.number = atoi((CrystalFolders[iCry].substr(crystal_prefix.size()+1,CrystalFolders[iCry].size()-crystal_prefix.size()-1)).c_str());

       TList *listCry = gDirectory->GetListOfKeys();
       int nKeysCry = listCry->GetEntries();
       std::vector<std::string> keysCryName;
       if(nKeysCry) //if directory not empty
       {

         for(int i = 0 ; i < nKeysCry ; i++){
           keysCryName.push_back(listCry->At(i)->GetName());
         }
         std::string CalibName;
         std::string CutName;
         std::vector<std::string> cutgNames;
         std::string cutG_prefix("cutg");
         std::string calibration_prefix("Calibration");
         std::string crystalCut_prefix("CrystalCut");
         std::string crystalCutWithoutCutG_prefix("CrystalCutWithoutCutG");
         std::string photopeakEnergyCut_prefix("PhotopeakEnergyCut");
         std::string channel_prefix("digitizerChannel");
         std::string w_channels_prefix("channelsNumRelevantForW");
         std::string timing_channel_prefix("timingChannel");
         std::string wz_prefix("w(z)");
         std::string t_channels_for_poli_prefix("tChannelsForPolishedCorrection");
         std::string mean_for_poli_prefix("meanForPolishedCorrection");
         std::string fwhm_for_poli_prefix("fwhmForPolishedCorrection");

         std::string lightCentral_prefix("Light collected in trigger crystal");
         std::string lightAll_prefix("Sum spectrum highlighted");
         std::string basicCTR_prefix("Basic CTR histogram");
        //  std::string correction_prefix("Correction Graph");
        //  std::string correction_rms_prefix("RMS Correction Graphs");
         for(unsigned int i = 0 ; i < keysCryName.size() ; i++)
         {
           if(!keysCryName[i].compare(0,calibration_prefix.size(),calibration_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.calibrationGraph = calibGraph;
           }

           if(!keysCryName[i].compare(0,wz_prefix.size(),wz_prefix)) //find calibration graph
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCanvas* C_graph = NULL;
             TGraph *calibGraph = NULL;
             C_graph = (TCanvas*) gDirectory->Get(keysCryName[i].c_str());
             if(C_graph)
               calibGraph = (TGraph*) C_graph->GetPrimitive(keysCryName[i].c_str());
             if(calibGraph)
               temp_crystal.wz = calibGraph;
           }

           if(!keysCryName[i].compare(0,crystalCut_prefix.size(),crystalCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCut = cut;
           }
           if(!keysCryName[i].compare(0,crystalCutWithoutCutG_prefix.size(),crystalCutWithoutCutG_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.CrystalCutWithoutCutG = cut;
           }
           if(!keysCryName[i].compare(0,photopeakEnergyCut_prefix.size(),photopeakEnergyCut_prefix)) // find tcut
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
             if(cut)
               temp_crystal.PhotopeakEnergyCut = cut;
           }

           if(!keysCryName[i].compare(0,cutG_prefix.size(),cutG_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TCutG* cutg = (TCutG*) gDirectory->Get(keysCryName[i].c_str());

             temp_crystal.cutg.push_back(cutg);
           }

           if(!keysCryName[i].compare(0,channel_prefix.size(),channel_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.detectorChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,lightCentral_prefix.size(),lightCentral_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightCentralHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,lightAll_prefix.size(),lightAll_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.lightAllHisto = aHisto;
           }

           if(!keysCryName[i].compare(0,basicCTR_prefix.size(),basicCTR_prefix)) // find tcutgs
           {
            //  std::cout << keysCryName[i] << std::endl;
             TH1F* aHisto = (TH1F*) gDirectory->Get(keysCryName[i].c_str());
             temp_crystal.basicCTRhisto = aHisto;
           }







           if(!keysCryName[i].compare(0,timing_channel_prefix.size(),timing_channel_prefix)) // find timing channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.timingChannel = atoi(snameCh.str().c_str());
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,w_channels_prefix.size(),w_channels_prefix)) // find detector channel
           {
             //  std::cout << keysCryName[i] << std::endl;
             std::stringstream snameCh;
             std::vector<int> *v;
             gDirectory->GetObject("channelsNumRelevantForW",v);
             // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
             //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
            //  istringstream()
             temp_crystal.relevantForW = v[0];
             //  std::cout <<temp_crystal.detectorChannel << std::endl;
            //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                      //  << temp_crystal.detectorChannel << std::endl;
           }

           if(!keysCryName[i].compare(0,t_channels_for_poli_prefix.size(),t_channels_for_poli_prefix))
           {
             std::vector<int> *v;
             gDirectory->GetObject("tChannelsForPolishedCorrection",v);
             temp_crystal.polishedCorrection = true;
             temp_crystal.tChannelsForPolishedCorrection = v[0];
           }

           if(!keysCryName[i].compare(0,mean_for_poli_prefix.size(),mean_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("meanForPolishedCorrection",v);
             temp_crystal.meanForPolishedCorrection = v[0];
           }

           if(!keysCryName[i].compare(0,fwhm_for_poli_prefix.size(),fwhm_for_poli_prefix))
           {
             std::vector<double> *v;
             gDirectory->GetObject("fwhmForPolishedCorrection",v);
             temp_crystal.fwhmForPolishedCorrection = v[0];
           }



         }

         bool dirExists;
         std::stringstream sname;
         dirExists = gDirectory->cd("TimeCorrection");
         if(dirExists)
         {
           sname.str("");
           sname << "Central correction - Crystal " << temp_crystal.number;
           temp_crystal.centralCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           sname.str("");
           sname << "Full correction - Crystal " << temp_crystal.number;
           temp_crystal.allCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
           sname.str("");

           if(likelihood)
           {
             sname.str("");
             sname << "Likelihood correction - Crystal " << temp_crystal.number;
             temp_crystal.likeCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
             sname.str("");
           }
           if(hybridCorrection)
           {
             sname.str("");
             sname << "Hybrid correction - Crystal " << temp_crystal.number;
             temp_crystal.hybridCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
             sname.str("");
           }

           temp_crystal.path = gDirectory->GetPath();
           TList *listTcorr = gDirectory->GetListOfKeys();
           int nKeysTcorr = listTcorr->GetEntries();
           std::vector<std::string> keysTcorrName;
           if(nKeysTcorr) //if directory not empty
           {
             for(int i = 0 ; i < nKeysTcorr ; i++){
               keysTcorrName.push_back(listTcorr->At(i)->GetName());
             }

             std::string deltaWGraph_prefix = "DeltaW Graph";
             std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
             std::string graph_delay_prefix = "Graph Delay ch_";
             std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
             std::string delay_timing_ch_prefix = "delayTimingChannels";
             for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
             {
               if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   temp_crystal.tw_correction = calibGraph;
                   //fit with straight line
                   TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                   calibGraph->Fit(line,"Q");
                   temp_crystal.tw_correction_line = line;

                 }

               }

               if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
               {
                TGraphErrors *calibGraph = NULL;
                calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                if(calibGraph)
                {
                  temp_crystal.rms_tw_correction = calibGraph;
                  //fit with straight line
                  TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                  calibGraph->Fit(line,"Q");
                  temp_crystal.rms_tw_correction_line = line;

                }

               }

               if(!keysTcorrName[i].compare(0,graph_delay_prefix.size(),graph_delay_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());
                 if(calibGraph)
                 {
                   // -- check if the channel is not excluded by the user
                   // extract next 2 characters
                   //
                   std::string str2 = keysTcorrName[i].substr (graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                   std::size_t found = str2.find_first_of("_");                                  // find next "_"
                   std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                   // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                   int current_ch = atoi(str3.c_str());                                          // transform in int
                   // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output

                   bool acceptCh = true;
                   for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   {
                     if(current_ch == forbidden_channels[iForb])
                     {
                       acceptCh = false;
                     }
                   }

                   if(acceptCh)               // add graph if the ch is accepted
                   {
                     temp_crystal.delay.push_back(calibGraph);
                     //fit with straight line
                     TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                     calibGraph->Fit(line,"Q");
                     temp_crystal.delay_line.push_back(line);
                   }



                 }

               }

               if(!keysTcorrName[i].compare(0,rms_graph_delay_prefix.size(),rms_graph_delay_prefix))
               {
                 TGraphErrors *calibGraph = NULL;
                 calibGraph = (TGraphErrors*) gDirectory->Get(keysTcorrName[i].c_str());


                 if(calibGraph)
                 {
                   // -- check if the channel is not excluded by the user
                   // extract next 2 characters
                   //
                   std::string str2 = keysTcorrName[i].substr (rms_graph_delay_prefix.size(),6);     // take a string with next 6 characters after the prefix
                   std::size_t found = str2.find_first_of("_");                                  // find next "_"
                   std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                   int current_ch = atoi(str3.c_str());                                          // transform in int
                   // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output


                   bool acceptCh = true;
                   for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                   {
                     if(current_ch == forbidden_channels[iForb])
                     {
                       acceptCh = false;
                     }
                   }

                   if(acceptCh)               // add graph if the ch is accepted
                   {
                     temp_crystal.rms_delay.push_back(calibGraph);
                     //fit with straight line
                     TF1 *line = new TF1("line",  "[0]*x + [1]",0,1);
                     calibGraph->Fit(line,"Q");
                     temp_crystal.rms_delay_line.push_back(line);
                   }
                 }

               }

               if(!keysTcorrName[i].compare(0,delay_timing_ch_prefix.size(),delay_timing_ch_prefix))
               {

                 // // -- check if the channel is not excluded by the user
                 // // extract next 2 characters
                 // //
                 // std::string str2 = keysTcorrName[i].substr (delay_timing_ch_prefix.size(),6);     // take a string with next 6 characters after the prefix
                 // std::size_t found = str2.find_first_of("_");                                  // find next "_"
                 // std::string str3 = str2.substr (0,found);                                     // extract everything before "_"
                 // // std::cout << keysTcorrName[i] << " " << str2 << " " << str3 << std::endl;     // output
                 // int current_ch = atoi(str3.c_str());                                          // transform in int
                 // std::cout << keysTcorrName[i] << "\t" << str2 << "\t" << str3 << "\t" << current_ch << std::endl;     // output
                 //
                 // bool acceptCh = true;
                 // for(int iForb = 0; iForb < forbidden_channels.size(); iForb++)                // check if this ch is in the forbidden_channels list
                 // {
                 //   if(current_ch == forbidden_channels[iForb])
                 //   {
                 //     acceptCh = false;
                 //   }
                 // }
                 //
                 // if(acceptCh)               // add graph if the ch is accepted
                 // {
                   //  std::cout << keysCryName[i] << std::endl;
                   std::stringstream snameCh;
                   std::vector<int> *v;
                   gDirectory->GetObject("delayTimingChannels",v);
                   // snameCh << ((TNamed*) gDirectory->Get(keysCryName[i].c_str()))->GetTitle();
                   //  TCut* cut = (TCut*) gDirectory->Get( keysCryName[i].c_str());
                  //  istringstream()
                   temp_crystal.delayTimingChannels = v[0];
                   //  std::cout <<temp_crystal.detectorChannel << std::endl;
                  //  std::cout << gDirectory->Get(keysCryName[i].c_str())->GetTitle() << "\t"
                            //  << temp_crystal.detectorChannel << std::endl;
                 // }

               }
             }
           }
           gDirectory->cd("..");
         }

         sname << "No correction - Crystal " << temp_crystal.number;
         temp_crystal.simpleCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
         sname.str("");
         sname << "Polished correction - Crystal " << temp_crystal.number;
         temp_crystal.poliCorrCTR = new TH1F(sname.str().c_str(),sname.str().c_str(),histoBins,histoMin,histoMax);
         sname.str("");
         TCut globalCut; // the cut for the formula
         globalCut += temp_crystal.CrystalCutWithoutCutG->GetTitle();     // this is BasicCut (XYZ and taggingPhotopeak) + CutTrigger (TriggerChannel and broadcut)
         globalCut += temp_crystal.PhotopeakEnergyCut->GetTitle();        // this is the cut on photopeak energy of the corrected spectrum for this crystal
         for(unsigned int iCutg = 0; iCutg < temp_crystal.cutg.size(); iCutg++)
         {
           globalCut += temp_crystal.cutg[iCutg]->GetName();              // these are the two cutg for this crystal
         }
         sname.str("");

         sname << "Formula" << temp_crystal.number;
         TTreeFormula* Formula = new TTreeFormula(sname.str().c_str(),globalCut,tree);
         formulas.Add(Formula);
         temp_crystal.Formula = Formula;
         sname.str("");

         // if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut && (temp_crystal.cutg.size() == 2))
         if(temp_crystal.calibrationGraph && temp_crystal.CrystalCutWithoutCutG && temp_crystal.PhotopeakEnergyCut)
         {
           crystal.push_back(temp_crystal);
         }

       }
       gDirectory->cd("..");
    }
    calibrationFile->cd("Module 0.0");
  }

  // //DEBUG
  // for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  // {
  //   std::cout << detectorSaturation[iSat].digitizerChannel << " " << detectorSaturation[iSat].saturation << " " << detectorSaturation[iSat].pedestal << std::endl;
  // }

  // COINCIDENCE CALIBRATION FILE
  //
  // TFile *coincidenceCalibrationFile = new TFile(coincidenceCalibrationFileName.c_str());
  // coincidenceCalibrationFile->cd("Module 0.0");
  // TList *listModuleCoinc = gDirectory->GetListOfKeys();
  // int nKeysModCoinc = listModuleCoinc->GetEntries();
  // std::vector<std::string> keysModNameCoinc;
  // // fill a vector with the leaves names
  // // std::string mppc_prefixCoinc("MPPC");
  // for(int i = 0 ; i < nKeysModCoinc ; i++){
  //   keysModNameCoinc.push_back(listModuleCoinc->At(i)->GetName());
  // }

//   TCut *taggingPhotopeakCut;
//   std::string taggingPhotopeakCut_prefix("taggingPhotopeakCut");
  // std::vector<std::string> MPPCfoldersCoinc;
//   for(unsigned int i = 0 ; i < keysModNameCoinc.size() ; i++)
//   {
    // if (!keysModNameCoinc[i].compare(0, mppc_prefixCoinc.size(), mppc_prefixCoinc))
    // {
    //   MPPCfoldersCoinc.push_back(keysModNameCoinc[i]);
    // }
//     if(!keysModNameCoinc[i].compare(0,taggingPhotopeakCut_prefix.size(),taggingPhotopeakCut_prefix)) // find tcut
//     {
     //  std::cout << keysCryName[i] << std::endl;
//       taggingPhotopeakCut = (TCut*) gDirectory->Get( keysModName[i].c_str());
      // if(cut)
      //   temp_crystal.CrystalCut = cut;
//     }
//   }

  // std::stringstream sformulaname;
  // sformulaname << "FormulaTag";
  // TCut taggingPhotopeakCutName;
  // taggingPhotopeakCutName = taggingPhotopeakCut->GetTitle();
  // TTreeFormula* FormulaTag = new TTreeFormula("FormulaTag",taggingPhotopeakCutName,tree);

  //look only into crystals found in the other calibration file
//   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
//   {
//     bool dirFound = gDirectory->cd(crystal[iCry].path);
//     if(dirFound)
//     {
//       crystal[iCry].accepted = true;
//       TList *listTcorr = gDirectory->GetListOfKeys();
//       int nKeysTcorr = listTcorr->GetEntries();
//       std::vector<std::string> keysTcorrName;
//       if(nKeysTcorr) //if directory not empty
//       {
//         for(int i = 0 ; i < nKeysTcorr ; i++){
//           keysTcorrName.push_back(listTcorr->At(i)->GetName());
//         }
//
//         std::string deltaWGraph_prefix = "DeltaW Graph";
//         std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";
//         // std::string graph_delay_prefix = "Graph Delay ch_";
//         // std::string rms_graph_delay_prefix = "RMS Graph Delay ch_";
//
//         for(unsigned int i = 0 ; i < keysTcorrName.size() ; i++)
//         {
//           if(!keysTcorrName[i].compare(0,deltaWGraph_prefix.size(),deltaWGraph_prefix))
//           {
//             TGraph *calibGraph = NULL;
//             calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//             if(calibGraph)
//               crystal[iCry].tw_correction = calibGraph;
//           }
//
//           if(!keysTcorrName[i].compare(0,rms_deltaWGraph_prefix.size(),rms_deltaWGraph_prefix))
//           {
//            TGraph *calibGraph = NULL;
//            calibGraph = (TGraph*) gDirectory->Get(keysTcorrName[i].c_str());
//            if(calibGraph)
//               crystal[iCry].rms_tw_correction = calibGraph;
//           }
//         }
//       }
//     }
//   }



  // // //BEGIN of DEBUG
  // for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  // {
  //   std::cout << crystal[i].number << "\t"
  //             // << crystal[i].cut->GetTitle() << "\t"
  //             << crystal[i].calibrationGraph->GetName() << "\t";
  //   for(unsigned int j = 0 ; j < crystal[i].cutg.size(); j++)
  //   {
  //     std::cout << crystal[i].cutg[j]->GetName() << "\t";
  //   }
  //   std::cout << std::endl;
  // }
  // // //END of DEBUG


  std::cout << "Calibration data found for crystals: " << std::endl;
  for(unsigned int i = 0 ;  i < crystal.size() ; i++)
  {
    if(crystal[i].accepted)
    {
      std::cout << crystal[i].number << std::endl;
    }
  }

  //notify TTreeFormula(s) to TChain
  tree->SetNotify(&formulas);

  //MAIN LOOP
  long long int nevent = tree->GetEntries();
  // ULong64_t tStart = tree->GetMinimum("ExtendedTimeTag");

  std::cout << "Total number of events = " << nevent << std::endl;
  long int goodEvents = 0;
  long int correlationMatrixEvents = 0;
  long int counter = 0;

  double tStart  = (tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t start in h
  double tEnd    = (tree->GetMaximum("ExtendedTimeTag") - tree->GetMinimum("ExtendedTimeTag"))/(1e9*3600); // t length in h

  double tStart2 = (tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);
  double tEnd2   = (tree->GetMaximum("DeltaTimeTag") - tree->GetMinimum("DeltaTimeTag"))/(1e9*3600);

  if(likelihood)
  {
    std::cout << "LIKELIHOOD CORRECTION" <<  std::endl;
    //create the correlation matrix and its inverse in each accepted crystal
    std::cout << "Preparing arrays... " <<  std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "Crystal " << crystal[iCry].number << std::endl;

        //count the delay plots, add 1, to get the likelihood channels
        int nLike = crystal[iCry].delay.size() + 1;
        float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
        float endW = crystal[iCry].wz->Eval(marginWZgraph);

        std::cout << "Preparing slices..." << std::endl;

        for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
        {
          slice_t temp_slice;

          Float_t wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
          Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
          Float_t wmean = (wmax + wmin) / 2.0;
          Float_t werr = (wmax-wmin)/TMath::Sqrt(12.0);

          temp_slice.wmean = wmean;
          temp_slice.werr = werr;
          temp_slice.wmin = wmin;
          temp_slice.wmax = wmax;
          temp_slice.entries = 0;
          temp_slice.covariance = new Float_t*[nLike];
          temp_slice.inverse_covariance = new Float_t*[nLike];
          temp_slice.entries_covariance = new long long int*[nLike];
          temp_slice.normalized_covariance = new Float_t*[nLike];

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            temp_slice.covariance[iCorr] = new Float_t[nLike];
            temp_slice.inverse_covariance[iCorr] = new Float_t[nLike];
            temp_slice.entries_covariance[iCorr] = new long long int[nLike];
            temp_slice.normalized_covariance[iCorr] = new Float_t[nLike];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              temp_slice.covariance[iCorr][jCorr] = 0;
              temp_slice.inverse_covariance[iCorr][jCorr] = 0;
              temp_slice.entries_covariance[iCorr][jCorr] = 0;
              temp_slice.normalized_covariance[iCorr][jCorr] = 0;
            }
            if(iCorr == 0) // crystal detector
            {
              temp_slice.tChannel.push_back(crystal[iCry].timingChannel);
              temp_slice.averageDelay.push_back(0); // no delay from cry to cry, by definition

            }
            else
            {
              std::string graphName = crystal[iCry].delay[iCorr-1]->GetName();
              std::size_t foundGraph = graphName.find_last_of("_");
              std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
              int graphCh = atoi(tChannelStringFromGraph.c_str() );
              temp_slice.tChannel.push_back(graphCh);
              temp_slice.averageDelay.push_back(crystal[iCry].delay_line[iCorr-1]->Eval(wmean)); // no delay from cry to cry, by definition

            }
            temp_slice.averageDeltaT.push_back(0);
            temp_slice.varianceDeltaT.push_back(0); // inizialize delta T to 0
            temp_slice.nVarianceDeltaT.push_back(0); // inizialize delta T to 0n
            temp_slice.nDeltaT.push_back(0); //inizialize deltaT to 0
          }

          crystal[iCry].slice.push_back(temp_slice);

          // prepare the covariance matrix, its inverse and the entries matrix
          // Float_t **temp_covariance;


        }
      }
    }

    std::cout << "Calculating average delta t for each slice..." << std::endl;
    for (long long int i=0;i<nevent;i++)
    {
      tree->GetEvent(i);              //read complete accepted event in memory
      //skip data if user say so
      bool keepEvent = true;
      if(sliced)
      {
        if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
        {
          keepEvent = false;
        }
      }
      else
      {
        if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
        {
          keepEvent = false;
        }
      }


      if(keepEvent)
      {
        for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
        {
          if(crystal[iCry].accepted)
          {
            if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
            {
              if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
              {
                //find w of this event
                //calculate FloodZ aka w
                Float_t FloodZ;
                float centralChargeOriginal;
                float centralSaturation;
                float centralPedestal;
                Float_t division = 0.0;

                centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                  {
                    centralSaturation = detectorSaturation[iSat].saturation;
                    centralPedestal = detectorSaturation[iSat].pedestal;
                  }
                }
                float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                {
                  // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
                  float originalCh = charge[crystal[iCry].relevantForW[iW]];

                  float saturationCh;
                  float pedestalCorr;
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                    {
                      saturationCh = detectorSaturation[iSat].saturation;
                      pedestalCorr = detectorSaturation[iSat].pedestal;
                    }
                  }
                  // std::cout << originalCh << " "
                  //           << saturationCh << " "
                  //           << pedestalCorr << " "
                  //           << std::endl;
                  division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                }

                FloodZ = centralChargeCorr / division;

                // std::cout << FloodZ << std::endl;



                for(int iSlice = 0 ; iSlice <  crystal[iCry].slice.size() ; iSlice++)
                {
                  if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) )
                  {
                    // this is the slice
                    // run on all detectors
                    int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                    for(int iCorr = 0; iCorr < nLike ; iCorr++)
                    {
                      if((timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] != 0) &&
                         (timeStamp[taggingCrystalTimingChannel] != 0))
                      {
                        Float_t deltaT = (timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[iCorr];
                        crystal[iCry].slice[iSlice].averageDeltaT[iCorr] += deltaT;
                        crystal[iCry].slice[iSlice].nDeltaT[iCorr] += 1;
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }

    std::cout << "Averaging the delta Ts" << std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            crystal[iCry].slice[iSlice].averageDeltaT[iCorr] = crystal[iCry].slice[iSlice].averageDeltaT[iCorr] / crystal[iCry].slice[iSlice].nDeltaT[iCorr];
          }
        }
      }
    }

    // //DEBUG
    std::cout << std::endl;
    std::cout << std::endl;
    std::cout << "|----------------------------------------|" << std::endl;
    std::cout << "|        After averaging deltaTs         |" << std::endl;
    std::cout << "|----------------------------------------|" << std::endl;
    std::cout << std::endl;
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
                      << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
                      << std::endl;
          }
        }
      }
    }
    //END OF DEBUG

// for(unsigned int jCorr = 0; jCorr < nLike ; jCorr++)
          // {
          //   temp_covariance[jCorr] = new Float_t[k];
          // }
          // for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
          // {
          //   for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
          //   {
          //     temp_covariance[iCorr][jCorr] = 0;
          //   }
          // }
    std::cout << "Calculating covariance matrix..." << std::endl;

    for (long long int i=0;i<nevent;i++)
    {
      tree->GetEvent(i);              //read complete accepted event in memory
      //skip data if user say so
      bool keepEvent = true;
      if(sliced)
      {
        if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
        {
          keepEvent = false;
        }
      }
      else
      {
        if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
        {
          keepEvent = false;
        }
      }


      if(keepEvent)
      {
        for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
        {
          if(crystal[iCry].accepted)
          {
            if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
            {
              if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
              {
                //find w of this event
                //calculate FloodZ aka w
                Float_t FloodZ;
                float centralChargeOriginal;
                float centralSaturation;
                float centralPedestal;
                Float_t division = 0.0;

                centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                  {
                    centralSaturation = detectorSaturation[iSat].saturation;
                    centralPedestal = detectorSaturation[iSat].pedestal;
                  }
                }
                float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                {
                  // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
                  float originalCh = charge[crystal[iCry].relevantForW[iW]];

                  float saturationCh;
                  float pedestalCorr;
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                    {
                      saturationCh = detectorSaturation[iSat].saturation;
                      pedestalCorr = detectorSaturation[iSat].pedestal;
                    }
                  }
                  // std::cout << originalCh << " "
                  //           << saturationCh << " "
                  //           << pedestalCorr << " "
                  //           << std::endl;
                  division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                }

                FloodZ = centralChargeCorr / division;

                // std::cout << FloodZ << std::endl;



                for(int iSlice = 0 ; iSlice <  crystal[iCry].slice.size() ; iSlice++)
                {
                  if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) )
                  {
                    // this is the slice
                    // run on all detectors
                    int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                    for(int iCorr = 0; iCorr < nLike ; iCorr++)
                    {
                      if((timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] != 0) &&
                         (timeStamp[taggingCrystalTimingChannel] != 0))
                      {

                        Float_t deltaT_i = (timeStamp[crystal[iCry].slice[iSlice].tChannel[iCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[iCorr];
                        Float_t element_i = deltaT_i - crystal[iCry].slice[iSlice].averageDeltaT[iCorr];
                        //update variance
                        crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] += pow(element_i,2);
                        crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr]++;

                        for(int jCorr = 0; jCorr < nLike ; jCorr++)
                        {
                          if((timeStamp[crystal[iCry].slice[iSlice].tChannel[jCorr]] != 0) &&
                             (timeStamp[taggingCrystalTimingChannel] != 0))
                          {

                            Float_t deltaT_j = (timeStamp[crystal[iCry].slice[iSlice].tChannel[jCorr]] - timeStamp[taggingCrystalTimingChannel]) - crystal[iCry].slice[iSlice].averageDelay[jCorr];
                            Float_t element_j = deltaT_j - crystal[iCry].slice[iSlice].averageDeltaT[jCorr];
                            crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] += element_i * element_j;
                            crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr]++;
                          }
                        }
                      }
                    }
                  }
                }
              }
            }
          }
        }
      }
    }


    std::cout << "Finalizing covariace matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // std::cout << "covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr];
            }
          }
        }
      }
    }


    // //DEBUG
    // for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    // {
    //   if(crystal[iCry].accepted)
    //   {
    //     std::cout << "crystal " << crystal[iCry].number << std::endl;
    //     for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
    //     {
    //       std::cout << "slice = " << iSlice << std::endl;
    //       std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
    //       std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
    //       std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
    //       std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
    //       std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
    //       int nLike = crystal[iCry].slice[iSlice].tChannel.size();
    //
    //       std::cout << "covariance" << std::endl;
    //
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "inverse_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "entries_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //
    //       std::cout << "normalized_covariance" << std::endl;
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         for(int jCorr = 0; jCorr < nLike ; jCorr++)
    //         {
    //           std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
    //         }
    //         std::cout << std::endl;
    //       }
    //       for(int iCorr = 0; iCorr < nLike ; iCorr++)
    //       {
    //         std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
    //                   << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
    //                   << std::endl;
    //       }
    //     }
    //   }
    // }
    // //END OF DEBUG

    std::cout << "Calculating normalized covariance matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // std::cout << "covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }
        }
      }
    }


    //DEBUG
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          // for(int iCorr = 0; iCorr < nLike ; iCorr++)
          // {
          //   std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
          //             << std::endl;
          // }
        }
      }
    }
    //END OF DEBUG


    std::cout << "Inverting covariance matrix..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        // std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          // std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();
          // use TMatrixF/D?
          // TMatrixD matrix(nLike,nLike); // never again, TMatrix simply doesn't work
          // use gls
          gsl_matrix *matrix = gsl_matrix_alloc(nLike, nLike);



          //copy covariance matrix to tmatrix
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              gsl_matrix_set(matrix,iCorr,jCorr,crystal[iCry].slice[iSlice].covariance[iCorr][jCorr]);
              // matrix[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr];
              // crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }

          size_t size = (size_t) nLike;

          gsl_matrix *inverse_matrix = invert_a_matrix(matrix,size);

          //invert tmatrix
          // double det;
          // matrix.Invert(&det);

          //copy matrix into inverted covariance
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            // crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] = crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] / crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr];
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
               crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] = gsl_matrix_get(inverse_matrix,iCorr,jCorr);
              // crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] = crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] / (TMath::Sqrt( crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] * crystal[iCry].slice[iSlice].varianceDeltaT[jCorr]));
            }
          }
        }
      }
    }


    //DEBUG
    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        std::cout << "crystal " << crystal[iCry].number << std::endl;
        for(int iSlice = 0 ; iSlice < WrangeBinsForTiming ; iSlice++)
        {
          std::cout << "slice = " << iSlice << std::endl;
          // std::cout << "wmin    = " << crystal[iCry].slice[iSlice].wmin << std::endl;
          // std::cout << "wmax    = " << crystal[iCry].slice[iSlice].wmax << std::endl;
          // std::cout << "wmean   = " << crystal[iCry].slice[iSlice].wmean << std::endl;
          // std::cout << "werr    = " << crystal[iCry].slice[iSlice].werr << std::endl;
          // std::cout << "entries = " << crystal[iCry].slice[iSlice].entries << std::endl;
          int nLike = crystal[iCry].slice[iSlice].tChannel.size();

          std::cout << "covariance" << std::endl;

          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "inverse_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "entries_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].entries_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          std::cout << "normalized_covariance" << std::endl;
          for(int iCorr = 0; iCorr < nLike ; iCorr++)
          {
            for(int jCorr = 0; jCorr < nLike ; jCorr++)
            {
              std::cout << crystal[iCry].slice[iSlice].normalized_covariance[iCorr][jCorr] << "\t";
            }
            std::cout << std::endl;
          }

          // for(int iCorr = 0; iCorr < nLike ; iCorr++)
          // {
          //   std::cout << crystal[iCry].slice[iSlice].tChannel[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDelay[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].averageDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].varianceDeltaT[iCorr] << "\t"
          //             << crystal[iCry].slice[iSlice].nVarianceDeltaT[iCorr] << "\t"
          //             << std::endl;
          // }
        }
      }
    }
    //END OF DEBUG


    std::cout << "Producing inverse covariance TGraphs..." << std::endl;

    for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
    {
      if(crystal[iCry].accepted)
      {
        //create the 2d matrix of tgraphs
        int nLike = crystal[iCry].slice[0].tChannel.size();
        int nSlice = crystal[iCry].slice.size();
        crystal[iCry].inverse_covariance_element = new TGraph**[nLike];
        crystal[iCry].inverse_covariance_element_line = new TF1**[nLike];
        for(int iCorr = 0; iCorr < nLike ; iCorr++)
        {
          crystal[iCry].inverse_covariance_element[iCorr] = new TGraph*[nLike];
          crystal[iCry].inverse_covariance_element_line[iCorr] = new TF1*[nLike];
          for(int jCorr = 0; jCorr < nLike ; jCorr++)
          {
            std::vector<Float_t> w;
            std::vector<Float_t> inv_s_ij;
            for(int iSlice = 0 ; iSlice < nSlice ; iSlice++)
            {
              w.push_back(crystal[iCry].slice[iSlice].wmean);
              inv_s_ij.push_back(crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr]);
            }
            crystal[iCry].inverse_covariance_element[iCorr][jCorr] = new TGraph(w.size(),&w[0],&inv_s_ij[0]);
            std::stringstream sname;
            sname << "cry_" << crystal[iCry].number <<"_inv_s_" << iCorr << "_" << jCorr;

            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->SetName(sname.str().c_str());
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->SetTitle(sname.str().c_str());
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->GetXaxis()->SetTitle("w");
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->GetYaxis()->SetTitle("inv_cov element");

            crystal[iCry].inverse_covariance_element_line[iCorr][jCorr] = new TF1("line","[0]*x + [1]",0,1);
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Fit(crystal[iCry].inverse_covariance_element_line[iCorr][jCorr],"Q");

          }
        }

      }
    }


  }// end likelihood








  //   std::cout << "Calculating delay(doi) plots... ";
  //   //first, likelihood part if given
  //   for (long long int i=0;i<nevent;i++)
  //   {
  //     // std::cout << "Event " << i << std::endl;
  //     tree->GetEvent(i);              //read complete accepted event in memory
  //     //skip data if user say so
  //     bool keepEvent = true;
  //     if(sliced)
  //     {
  //       if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //     else
  //     {
  //       if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //
  //
  //     if(keepEvent)
  //     {
  //       for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  //       {
  //         if(crystal[iCry].accepted)
  //         {
  //           if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
  //           {
  //             if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
  //             {
  //
  //               correlationMatrixEvents++;
  //
  //               if(crystal[iCry].tw_correction)
  //               {
  //
  //                 //calculate FloodZ aka w
  //                 Float_t FloodZ;
  //                 float centralChargeOriginal;
  //                 float centralSaturation;
  //                 float centralPedestal;
  //                 Float_t division = 0.0;
  //
  //                 centralChargeOriginal = charge[crystal[iCry].detectorChannel];
  //                 for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                 {
  //                   if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
  //                   {
  //                     centralSaturation = detectorSaturation[iSat].saturation;
  //                     centralPedestal = detectorSaturation[iSat].pedestal;
  //                   }
  //                 }
  //                 float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
  //
  //                 for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
  //                 {
  //                   // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
  //                   float originalCh = charge[crystal[iCry].relevantForW[iW]];
  //
  //                   float saturationCh;
  //                   float pedestalCorr;
  //                   for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                   {
  //                     if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
  //                     {
  //                       saturationCh = detectorSaturation[iSat].saturation;
  //                       pedestalCorr = detectorSaturation[iSat].pedestal;
  //                     }
  //                   }
  //                   // std::cout << originalCh << " "
  //                   //           << saturationCh << " "
  //                   //           << pedestalCorr << " "
  //                   //           << std::endl;
  //                   division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  //                 }
  //
  //                 FloodZ = centralChargeCorr / division;
  //
  //                 // std::cout << FloodZ << std::endl;
  //
  //                 unsigned int k = crystal[iCry].LikelihoodChannels.size();
  //                 for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //                 {
  //
  //                   if( (timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] != 0) &&
  //                       (timeStamp[taggingCrystalTimingChannel] != 0))
  //                   {
  //
  //                     Float_t iTimeStamp = timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] - timeStamp[taggingCrystalTimingChannel];
  //                     // std::cout << iTimeStamp << std::endl;
  //                     Float_t iDelay = 0;
  //                     if(iCorr == 0)
  //                     {
  //                       iDelay = 0;
  //                     }
  //                     else
  //                     {
  //                       iDelay = crystal[iCry].LikelihoodChannels[iCorr].delay_line->Eval(FloodZ);
  //                     }
  //                     crystal[iCry].LikelihoodChannels[iCorr].deltaTscatter->Fill(FloodZ,iTimeStamp-iDelay);
  //                   }
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //     counter++;
  //
  //
  //   }
  //   // std::cout << "\r";
  //   std::cout << "done!" << std::endl;
  //   // std::cout << "Events for covariance matrix = " << correlationMatrixEvents << std::endl;
  //
  //
  //
  //   //slice the th2f plots
  //   std::cout << "Creating slice histrograms..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //
  //       float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
  //       float endW = crystal[iCry].wz->Eval(marginWZgraph);
  //       for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //       {
  //         Float_t wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
  //         Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
  //         Float_t wmean = (wmax + wmin) / 2.0;
  //         std::stringstream sname;
  //         sname << "T" << crystal[iCry].LikelihoodChannels[iCorr].tChannel << "_like_slide_cry" <<  crystal[iCry].number << "_w_" << wmin << "_" << wmax;
  //         TH1F* temp_histo = new TH1F(sname.str().c_str(),sname.str().c_str(),likeBins,likeMin,likeMax);
  //         crystal[iCry].LikelihoodChannels[iCorr].deltaTslice.push_back(temp_histo);
  //         crystal[iCry].LikelihoodChannels[iCorr].wmean.push_back(wmean);
  //         crystal[iCry].LikelihoodChannels[iCorr].werr.push_back((wmax-wmin)/TMath::Sqrt(12.0));
  //         // crystal[iCry].LikelihoodChannels[iCorr].dx.push_back(wmean);
  //         // crystal[iCry].LikelihoodChannels[iCorr].dex.push_back((wmax-wmin)/TMath::Sqrt(12.0));
  //       }
  //     }
  //
  //
  //
  //   }
  //
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //     {
  //       unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //       Float_t **temp_covariance;
  //       temp_covariance = new Float_t*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         temp_covariance[jCorr] = new Float_t[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           temp_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].covariance.push_back(temp_covariance);
  //
  //       long long int **temp_entries_covariance;
  //       temp_entries_covariance = new long long int*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         temp_entries_covariance[jCorr] = new long long int[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           temp_entries_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].entries_covariance.push_back(temp_entries_covariance);
  //
  //       Float_t **inverse_covariance;
  //       inverse_covariance = new Float_t*[k];
  //       for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //       {
  //         inverse_covariance[jCorr] = new Float_t[k];
  //       }
  //       for(unsigned int iCorr = 0; iCorr < k ; iCorr++)
  //       {
  //         for(unsigned int jCorr = 0; jCorr < k ; jCorr++)
  //         {
  //           inverse_covariance[iCorr][jCorr] = 0;
  //         }
  //       }
  //       crystal[iCry].inverse_covariance.push_back(inverse_covariance);
  //     }
  //   }
  //
  //
  //
  //   correlationMatrixEvents = 0;
  //   counter = 0;
  //   std::cout << "Filling slice histrograms..." << std::endl;
  //   for (long long int i=0;i<nevent;i++)
  //   {
  //     // std::cout << "Event " << i << std::endl;
  //     tree->GetEvent(i);              //read complete accepted event in memory
  //     //skip data if user say so
  //     bool keepEvent = true;
  //     if(sliced)
  //     {
  //       if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //     else
  //     {
  //       if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
  //       {
  //         keepEvent = false;
  //       }
  //     }
  //
  //
  //     if(keepEvent)
  //     {
  //       for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  //       {
  //         if(crystal[iCry].accepted)
  //         {
  //           if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
  //           {
  //             if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
  //             {
  //
  //               correlationMatrixEvents++;
  //
  //               if(crystal[iCry].tw_correction)
  //               {
  //
  //                 //calculate FloodZ aka w
  //                 Float_t FloodZ;
  //                 float centralChargeOriginal;
  //                 float centralSaturation;
  //                 float centralPedestal;
  //                 Float_t division = 0.0;
  //
  //                 centralChargeOriginal = charge[crystal[iCry].detectorChannel];
  //                 for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                 {
  //                   if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
  //                   {
  //                     centralSaturation = detectorSaturation[iSat].saturation;
  //                     centralPedestal = detectorSaturation[iSat].pedestal;
  //                   }
  //                 }
  //                 float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );
  //
  //                 for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
  //                 {
  //                   // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
  //                   float originalCh = charge[crystal[iCry].relevantForW[iW]];
  //
  //                   float saturationCh;
  //                   float pedestalCorr;
  //                   for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
  //                   {
  //                     if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
  //                     {
  //                       saturationCh = detectorSaturation[iSat].saturation;
  //                       pedestalCorr = detectorSaturation[iSat].pedestal;
  //                     }
  //                   }
  //                   // std::cout << originalCh << " "
  //                   //           << saturationCh << " "
  //                   //           << pedestalCorr << " "
  //                   //           << std::endl;
  //                   division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
  //                 }
  //
  //                 FloodZ = centralChargeCorr / division;
  //
  //
  //                 float beginW = crystal[iCry].wz->Eval(length - marginWZgraph);
  //                 float endW = crystal[iCry].wz->Eval(marginWZgraph);
  //                 for(int iBin = 0; iBin < WrangeBinsForTiming; iBin++) //
  //                 {
  //                   Float_t wmin = beginW + ((iBin*(endW - beginW))/WrangeBinsForTiming);
  //                   Float_t wmax = beginW + (((iBin+1)*(endW - beginW))/WrangeBinsForTiming);
  //
  //                   if((FloodZ >=wmin) && (FloodZ < wmax) )
  //                   {
  //                     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //                     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //                     {
  //                       if( (timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] != 0) &&
  //                       (timeStamp[taggingCrystalTimingChannel] != 0))
  //                       {
  //                         Float_t iTimeStamp = timeStamp[crystal[iCry].LikelihoodChannels[iCorr].tChannel] - timeStamp[taggingCrystalTimingChannel];
  //
  //                         Float_t iDelay = 0;
  //                         if(iCorr == 0)
  //                         {
  //                           iDelay = 0;
  //                         }
  //                         else
  //                         {
  //                           iDelay = crystal[iCry].LikelihoodChannels[iCorr].delay_line->Eval(FloodZ);
  //                         }
  //                         crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iBin]->Fill(iTimeStamp-iDelay);
  //                       }
  //                     }
  //                   }
  //                 }
  //               }
  //             }
  //           }
  //         }
  //       }
  //     }
  //     counter++;
  //
  //     // int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
  //     // if( (perc % 10) == 0 )
  //     // {
  //     //   std::cout << "\r";
  //     //   std::cout << "Calculating covariance matrix... " << perc << "%";
  //     //   //std::cout << counter << std::endl;
  //     // }
  //   }
  //   // std::cout << "\r";
  //   std::cout << "done!" << std::endl;
  //   // std::cout << "Events for covariance matrix = " << correlationMatrixEvents << std::endl;
  //
  //
  //   std::cout << "fitting slice histograms..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //       for(unsigned int iSlice = 0 ; iSlice < crystal[iCry].LikelihoodChannels[iCorr].deltaTslice.size(); iSlice++ )
  //       {
  //         // crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iSlice];
  //         double res[4];
  //         int divisions = 100; //useless
  //         extractFromHisto(crystal[iCry].LikelihoodChannels[iCorr].deltaTslice[iSlice],fitPercMin,fitPercMax,divisions,res);
  //
  //         if(res[0] == 0 && res[1] == 0 && res[2] == 0 && res[3] == 0) //ignore point if fit didn't work
  //         {
  //                                 // skip point
  //         }
  //         else
  //         {
  //
  //           crystal[iCry].LikelihoodChannels[iCorr].dx.push_back(crystal[iCry].LikelihoodChannels[iCorr].wmean[iSlice]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dex.push_back(crystal[iCry].LikelihoodChannels[iCorr].werr[iSlice]);
  //           // crystal[iCry].LikelihoodChannels[iCorr].dex.push_back(res[2]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dy.push_back(res[0]);
  //           crystal[iCry].LikelihoodChannels[iCorr].dey.push_back(res[2]);
  //         }
  //       }
  //     }
  //   }
  //
  //
  //
  //   std::cout << "generating and fitting tgraphs..." << std::endl;
  //   for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++) // for each crystal
  //   {
  //     unsigned int k = crystal[iCry].LikelihoodChannels.size(); // for each 2d plot
  //     for(unsigned int iCorr = 0 ; iCorr < k; iCorr++)
  //     {
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph = new TGraphErrors(
  //         crystal[iCry].LikelihoodChannels[iCorr].dx.size(),
  //         &crystal[iCry].LikelihoodChannels[iCorr].dx[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dy[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dex[0],
  //         &crystal[iCry].LikelihoodChannels[iCorr].dey[0]
  //       );
  //       std::stringstream sname;
  //       sname << "deltaT graph cry " <<  crystal[iCry].number << "_T" <<  crystal[iCry].LikelihoodChannels[iCorr].tChannel;
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->SetName(sname.str().c_str());
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->SetTitle(sname.str().c_str());
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->GetXaxis()->SetTitle("W");
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->GetYaxis()->SetTitle("DeltaT");
  //
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline = new TF1("deltaTline",  "[0]*x + [1]",0,1);
  //       float m = (crystal[iCry].LikelihoodChannels[iCorr].dy[0]-
  //                 crystal[iCry].LikelihoodChannels[iCorr].dy[crystal[iCry].LikelihoodChannels[iCorr].dx.size()-1])/
  //                 (crystal[iCry].LikelihoodChannels[iCorr].dx[0]-
  //                 crystal[iCry].LikelihoodChannels[iCorr].dx[crystal[iCry].LikelihoodChannels[iCorr].dx.size()-1]);
  //       float q = crystal[iCry].LikelihoodChannels[iCorr].dy[0] -
  //                 m*crystal[iCry].LikelihoodChannels[iCorr].dx[0];
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline->SetParameter(0,m);
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTline->SetParameter(1,q);
  //
  //       crystal[iCry].LikelihoodChannels[iCorr].deltaTgraph->Fit(crystal[iCry].LikelihoodChannels[iCorr].deltaTline,"");
  //     }
  //   }
  //
  //   std::cout << "generating covariance scatter plots..." << std::endl;
  //
  //
  //
  //
  //

  counter = 0;











  // for (long long int i=0;i<1000000;i++)
  for (long long int i=0;i<nevent;i++)
  {
    // std::cout << "Event " << i << std::endl;
    tree->GetEvent(i);              //read complete accepted event in memory
    //skip data if user say so
    bool keepEvent = true;
    if(sliced)
    {
      if( ((ChainExtendedTimeTag / (1e9*3600) ) - tStart) < start_time)
      {
        keepEvent = false;
      }
    }
    else
    {
      if( ((ChainDeltaTimeTag    / (1e9*3600) ) - tStart2) < start_time)
      {
        keepEvent = false;
      }
    }


    if(keepEvent)
    {
      for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
      {
        if(crystal[iCry].accepted)
        {
          if(FormulaTag->EvalInstance() || simulation) // if in photopeak of tagging crystal - or if in simulation
          {
            if(crystal[iCry].Formula->EvalInstance())  //if in global cut of crystal
            {

              goodEvents++;

              //temp commented
              Float_t centralcorrection = 0.0;
              Float_t zeroCorrection    = 0.0;
              //no corr
              double simpleCTR = timeStamp[crystal[iCry].timingChannel] - timeStamp[taggingCrystalTimingChannel];
              if((timeStamp[crystal[iCry].timingChannel] != 0) && (timeStamp[taggingCrystalTimingChannel] != 0)) // no zeroes
              {
                crystal[iCry].simpleCTR->Fill(simpleCTR);
                crystal[iCry].vSimple.push_back(simpleCTR);
              }

              if(crystal[iCry].tw_correction)
              {

                //calculate FloodZ...
                Float_t FloodZ;
                float centralChargeOriginal;
                float centralSaturation;
                float centralPedestal;
                Float_t division = 0.0;

                centralChargeOriginal = charge[crystal[iCry].detectorChannel];
                for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                {
                  if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].detectorChannel)
                  {
                    centralSaturation = detectorSaturation[iSat].saturation;
                    centralPedestal = detectorSaturation[iSat].pedestal;
                  }
                }
                float centralChargeCorr = ( -centralSaturation * TMath::Log(1.0 - ( ( (centralChargeOriginal-centralPedestal))/(centralSaturation)) ) );

                for (unsigned int iW = 0; iW < crystal[iCry].relevantForW.size(); iW++)
                {
                  // std::cout << crystal[iCry].relevantForW[iW] << std::endl;
                  float originalCh = charge[crystal[iCry].relevantForW[iW]];

                  float saturationCh;
                  float pedestalCorr;
                  for(unsigned int iSat = 0; iSat < detectorSaturation.size(); iSat++)
                  {
                    if( detectorSaturation[iSat].digitizerChannel  == crystal[iCry].relevantForW[iW])
                    {
                      saturationCh = detectorSaturation[iSat].saturation;
                      pedestalCorr = detectorSaturation[iSat].pedestal;
                    }
                  }
                  // std::cout << originalCh << " "
                  //           << saturationCh << " "
                  //           << pedestalCorr << " "
                  //           << std::endl;
                  division += ( -saturationCh * TMath::Log(1.0 - ( ( (originalCh-pedestalCorr))/(saturationCh)) ) );
                }

                FloodZ = centralChargeCorr / division;

                // std::cout << centralChargeOriginal << " "
                //           << centralSaturation << " "
                //           << centralChargeCorr << " "
                //           << division << " "
                //           << FloodZ << " "
                //           << std::endl;

                // central corr
                // std::string deltaWGraph_prefix = "DeltaW Graph";
                // std::string rms_deltaWGraph_prefix = "RMS DeltaW Graph";


                Float_t averageTimeStamp = 0.0;
                Float_t totalWeight = 0.0;
                // averageTimeStamp += timeStamp[crystal[iCry].detectorChannel];

                double centralCTR;

                if(fitCorrection)
                {
                  centralcorrection = crystal[iCry].tw_correction_line->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction_line->Eval(FloodZ);
                }
                else
                {
                  centralcorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
                }

                centralCTR = (timeStamp[crystal[iCry].timingChannel] + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel];
                crystal[iCry].centralCTR->Fill(centralCTR);
                crystal[iCry].vCentral.push_back(centralCTR);


                if(crystal[iCry].delay.size())
                {
                  // full corr v0
                  // first evalutate all ts as if they where read by the central crystal, i.e. correct all the neighboring channels
                  // using the delta
                  //
                  Float_t weight = 0.0;

                  if(fitCorrection)
                  {
                    weight = pow(crystal[iCry].rms_tw_correction_line->Eval(FloodZ),-2);
                  }
                  else
                  {
                    weight = pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),-2);
                  }


                  averageTimeStamp += weight * timeStamp[crystal[iCry].timingChannel];
                  totalWeight += weight;
                  //
                  //full corr
                  // zeroCorrection = crystal[iCry].tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)) - crystal[iCry].tw_correction->Eval(FloodZ);
                  //
                  // Float_t weight = 0.0;
                  // weight = pow(sqrt(pow(crystal[iCry].rms_tw_correction->Eval(FloodZ),2)+pow(crystal[iCry].rms_tw_correction->Eval(crystal[iCry].wz->Eval(length*doiFraction)),2)),-2);
                  //
                  // averageTimeStamp += weight*(timeStamp[crystal[iCry].timingChannel]- zeroCorrection);
                  // totalWeight += weight;
                  // // std::cout << i << " " << iCry << " " <<  crystal[iCry].number << "\n";
                  for(unsigned int iGraph = 0; iGraph < crystal[iCry].delay.size();iGraph++)
                  {

                    std::string graphName = crystal[iCry].delay[iGraph]->GetName();
                    std::string rmsName = crystal[iCry].rms_delay[iGraph]->GetName();

                    std::size_t foundGraph = graphName.find_last_of("_");
                    std::size_t foundRms   = rmsName.find_last_of("_");

                    std::string tChannelStringFromGraph = graphName.substr(foundGraph+1);
                    std::string tChannelStringFromRms   = rmsName.substr(foundRms  +1);

                    // std::stringstream sname;
                    // sname << "Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
                    // std::string graph_delay_prefix = sname.str();
                    // std::cout << tChannelStringFromGraph << std::endl;
                    // std::cout << tChannelStringFromRms << std::endl;

                    // sname.str("");
                    // sname << "RMS Graph Delay ch_" << crystal[iCry].detectorChannel << "_t_" ;
                    // std::string rms_graph_delay_prefix = sname.str();
                    // std::cout << tChannelStringFromRms << std::endl;
                    // sname.str("");






                    int graphCh = atoi(tChannelStringFromGraph.c_str() );
                    // std::cout << graphCh  << "\n";


                    int rmsCh   = atoi( tChannelStringFromRms.c_str() );
                    if(graphCh != rmsCh)
                    {
                      std::cout << "ERROR! TGraphs of delay and rms are from different timing channels!!!!" << std::endl;
                      break;
                    }
                    // std::cout << rmsCh  << "\n";

                    if(fitCorrection)
                    {
                      weight = pow(crystal[iCry].rms_delay_line[iGraph]->Eval(FloodZ),-2);
                    }
                    else
                    {
                      weight = pow(crystal[iCry].rms_delay[iGraph]->Eval(FloodZ),-2);
                    }

                    totalWeight += weight;
                    if(fitCorrection)
                    {
                      averageTimeStamp += weight*(timeStamp[graphCh] - crystal[iCry].delay_line[iGraph]->Eval(FloodZ));
                    }
                    else
                    {
                      averageTimeStamp += weight*(timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ));
                    }

                    // std::cout << timeStamp[graphCh] - crystal[iCry].delay[iGraph]->Eval(FloodZ) << "\t";
                  }
                  averageTimeStamp = averageTimeStamp/totalWeight;
                  // then correct the average of ts for central correction
                  // double allCTR = (averageTimeStamp + (centralcorrection)) - timeStamp[taggingCrystalTimingChannel];
                  double allCTR = (averageTimeStamp + (0.0)) - timeStamp[taggingCrystalTimingChannel];
                  crystal[iCry].allCTR->Fill(allCTR);
                  crystal[iCry].vAll.push_back(allCTR);


                  if(likelihood)
                  {
                    //combine the detector results using inverse covariance matrix elements
                    // calc the sum of all elements

                    //first check the no timestamp is 0
                    bool noZeroes = true;
                    if(timeStamp[taggingCrystalTimingChannel] == 0)
                    {
                      noZeroes = false;   // and you can already stop
                    }
                    else // look for all the other channels involved
                    {
                      for(unsigned int iCorr = 0; iCorr < crystal[iCry].slice[0].tChannel.size() ; iCorr++)
                      {
                        if(timeStamp[crystal[iCry].slice[0].tChannel[iCorr]] == 0)
                        {
                          noZeroes = false;
                        }
                      }
                    }

                    if(noZeroes)
                    {
                      if(basicLikelihood)
                      {
                        // zero approach: just find the w slice and use that inv cov matrix
                        int nSlice = crystal[iCry].slice.size();
                        for(int iSlice = 0 ; iSlice < nSlice ; iSlice++)
                        {
                          if((FloodZ >= crystal[iCry].slice[iSlice].wmin )&&( FloodZ < crystal[iCry].slice[iSlice].wmax ) ) // this is the w slice
                          {

                              Float_t sum_inv_cov = 0;
                              Float_t Dt_best = 0;
                              int nLike = crystal[iCry].slice[iSlice].tChannel.size();
                              for(int jCorr = 0; jCorr < nLike ; jCorr++)
                              {
                                int tChannel = crystal[iCry].slice[iSlice].tChannel[jCorr];
                                // find delay as function of doi
                                Float_t delay;
                                if(jCorr == 0)
                                {
                                  delay = 0; // no delay for direct detector
                                }
                                else
                                {
                                  delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                                }
                                // find Dt of this detector (t - tR) - delay;
                                Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                                // now calc the weight_i
                                // keeping j fixed, run on i and sum the inv cov elements
                                Float_t weight_i = 0;
                                for(int iCorr = 0; iCorr < nLike ; iCorr++)
                                {
                                  weight_i += crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr];
                                  // weight_i +=  crystal[iCry].inverse_covariance_element[iCorr][jCorr]
                                  // at the same time also sum all the elements of inv cov matrix to get the final division
                                  // sum_inv_cov += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                                  sum_inv_cov += crystal[iCry].slice[iSlice].inverse_covariance[iCorr][jCorr];
                                }
                                Dt_best += weight_i * Dt_j;
                              }

                              Dt_best = Dt_best / sum_inv_cov;
                              double likeCTR = Dt_best + centralcorrection ;
                              crystal[iCry].likeCTR->Fill(likeCTR);
                              crystal[iCry].vLike.push_back(likeCTR);
                          }
                        }
                      }
                      else
                      {
                        if(likelihoodLine)
                        {
                          Float_t sum_inv_cov = 0;
                          Float_t Dt_best = 0;
                          int nLike = crystal[iCry].slice[0].tChannel.size();
                          for(int jCorr = 0; jCorr < nLike ; jCorr++)
                          {
                            int tChannel = crystal[iCry].slice[0].tChannel[jCorr];
                            // find delay as function of doi
                            Float_t delay;
                            if(jCorr == 0)
                            {
                              delay = 0; // no delay for direct detector
                            }
                            else
                            {
                              delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                            }
                            // find Dt of this detector (t - tR) - delay;
                            Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                            // now calc the weight_i
                            // keeping j fixed, run on i and sum the inv cov elements
                            Float_t weight_i = 0;
                            for(int iCorr = 0; iCorr < nLike ; iCorr++)
                            {
                              weight_i += crystal[iCry].inverse_covariance_element_line[iCorr][jCorr]->Eval(FloodZ);
                              sum_inv_cov += crystal[iCry].inverse_covariance_element_line[iCorr][jCorr]->Eval(FloodZ);
                            }
                            Dt_best += weight_i * Dt_j;
                          }
                          Dt_best = Dt_best / sum_inv_cov;
                          double likeCTR = Dt_best + centralcorrection ;
                          crystal[iCry].likeCTR->Fill(likeCTR);
                          crystal[iCry].vLike.push_back(likeCTR);
                        }
                        else
                        {
                          Float_t sum_inv_cov = 0;
                          Float_t Dt_best = 0;
                          int nLike = crystal[iCry].slice[0].tChannel.size();
                          for(int jCorr = 0; jCorr < nLike ; jCorr++)
                          {
                            int tChannel = crystal[iCry].slice[0].tChannel[jCorr];
                            // find delay as function of doi
                            Float_t delay;
                            if(jCorr == 0)
                            {
                              delay = 0; // no delay for direct detector
                            }
                            else
                            {
                              delay = crystal[iCry].delay_line[jCorr-1]->Eval(FloodZ);
                            }
                            // find Dt of this detector (t - tR) - delay;
                            Float_t Dt_j = (timeStamp[tChannel] - timeStamp[taggingCrystalTimingChannel]) - delay;

                            // now calc the weight_i
                            // keeping j fixed, run on i and sum the inv cov elements
                            Float_t weight_i = 0;
                            for(int iCorr = 0; iCorr < nLike ; iCorr++)
                            {
                              weight_i += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                              sum_inv_cov += crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Eval(FloodZ);
                            }
                            Dt_best += weight_i * Dt_j;
                          }
                          Dt_best = Dt_best / sum_inv_cov;
                          double likeCTR = Dt_best + centralcorrection ;
                          crystal[iCry].likeCTR->Fill(likeCTR);
                          crystal[iCry].vLike.push_back(likeCTR);
                        }
                      }
                    }
                  }
                  // crystal[iCry].allCTR->Fill(averageTimeStamp + zeroCorrection  - timeStamp[taggingCrystalTimingChannel]);
                }
                // hybrid correction
                //
                if(hybridCorrection)
                {
                  //central time stamp
                  // std::cout << "----------------" << std::endl;
                  float weight = 0.0;
                  float meanTimeStamp = 0.0;
                  float sumWeight = 0.0;
                  // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;

                  weight = pow(crystal[iCry].fwhmForPolishedCorrection[0],-2);
                  float t_0 = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[0]];

                  meanTimeStamp += weight * t_0;
                  sumWeight += weight;
                  // std::cout << weight << " " << t_0 << " " << meanTimeStamp << " " << sumWeight << std::endl;

                  // std::cout << t_0 << "\t";

                  for(unsigned int iPoli = 1; iPoli < crystal[iCry].tChannelsForPolishedCorrection.size(); iPoli++)
                  {
                    // std::cout << iPoli << std::endl;
                    // std::cout << "timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] = " << timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] << std::endl;
                    // std::cout << "crystal[iCry].meanForPolishedCorrection[iPoli] = " << crystal[iCry].meanForPolishedCorrection[iPoli] << std::endl;

                    float t_i = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] - crystal[iCry].meanForPolishedCorrection[iPoli];
                    float weight_i = pow(crystal[iCry].fwhmForPolishedCorrection[iPoli],-2);
                    meanTimeStamp += weight_i * t_i;
                    sumWeight += weight_i;
                    // std::cout << weight_i << " " << t_i << " " << meanTimeStamp << " " << sumWeight << std::endl;
                    // std::cout << t_i << "\t";
                  }
                  // std::cout << std::endl;
                  meanTimeStamp = meanTimeStamp / sumWeight;
                  // std::cout << std::endl
                  //           << meanTimeStamp << "\t"
                  //           << std::endl;
                  double hybridCorrCTR = meanTimeStamp - timeStamp[taggingCrystalTimingChannel];
                  //correct by doi

                  hybridCorrCTR = hybridCorrCTR + centralcorrection; //FIXME temp

                  crystal[iCry].hybridCTR->Fill(hybridCorrCTR);
                  crystal[iCry].vhybrid.push_back(hybridCorrCTR);


                }



              }

              if(crystal[iCry].polishedCorrection)
              {
                //central time stamp
                // std::cout << "----------------" << std::endl;
                float weight = 0.0;
                float meanTimeStamp = 0.0;
                float sumWeight = 0.0;
                // std::cout << "crystal[iCry].fwhmForPolishedCorrection[0] = " << crystal[iCry].fwhmForPolishedCorrection[0]<< std::endl;

                weight = pow(crystal[iCry].fwhmForPolishedCorrection[0],-2);
                float t_0 = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[0]];

                meanTimeStamp += weight * t_0;
                sumWeight += weight;
                // std::cout << weight << " " << t_0 << " " << meanTimeStamp << " " << sumWeight << std::endl;

                // std::cout << t_0 << "\t";

                for(unsigned int iPoli = 1; iPoli < crystal[iCry].tChannelsForPolishedCorrection.size(); iPoli++)
                {
                  // std::cout << iPoli << std::endl;
                  // std::cout << "timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] = " << timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] << std::endl;
                  // std::cout << "crystal[iCry].meanForPolishedCorrection[iPoli] = " << crystal[iCry].meanForPolishedCorrection[iPoli] << std::endl;

                  float t_i = timeStamp[crystal[iCry].tChannelsForPolishedCorrection[iPoli]] - crystal[iCry].meanForPolishedCorrection[iPoli];
                  float weight_i = pow(crystal[iCry].fwhmForPolishedCorrection[iPoli],-2);
                  meanTimeStamp += weight_i * t_i;
                  sumWeight += weight_i;
                  // std::cout << weight_i << " " << t_i << " " << meanTimeStamp << " " << sumWeight << std::endl;
                  // std::cout << t_i << "\t";
                }
                // std::cout << std::endl;
                meanTimeStamp = meanTimeStamp / sumWeight;
                // std::cout << std::endl
                //           << meanTimeStamp << "\t"
                //           << std::endl;
                double poliCorrCTR = meanTimeStamp - timeStamp[taggingCrystalTimingChannel];
                crystal[iCry].poliCorrCTR->Fill(poliCorrCTR);
                crystal[iCry].vPoli.push_back(poliCorrCTR);


              }



              // end of temp commented
            }
          }
        }
        // delete Formula;
      }
    }

    counter++;

    int perc = ((100*counter)/nevent); //should strictly have not decimal part, written like this...
    if( (perc % 10) == 0 )
    {
      std::cout << "\r";
      std::cout << perc << "% done... ";
      //std::cout << counter << std::endl;
    }
  }
  std::cout << "Good events = " << goodEvents << std::endl;
  std::sort(crystal.begin(), crystal.end(), compareByNumber);


  TH1F* noCorr = new TH1F("No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* centralCorr = new TH1F("Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* fullCorr = new TH1F("Full Correction","Full Correction",bins,minCTR,maxCTR);
  TH1F* likeCorr = new TH1F("Likelihood Correction","Likelihood Correction",bins,minCTR,maxCTR);
  TH1F* poliCorr = new TH1F("Polished Correction","Polished Correction",bins,minCTR,maxCTR);
  TH1F* hybridCorr = new TH1F("Hybrid Correction","Hybrid Correction",bins,minCTR,maxCTR);

  TH1F* unbinnednoCorr      = new TH1F("Unbinned No Correction","No Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedcentralCorr = new TH1F("Unbinned Central Correction","Central Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedfullCorr    = new TH1F("Unbinned Full Correction","Full Correction",bins,minCTR,maxCTR);
  TH1F* unbinnedpoliCorr    = new TH1F("Unbinned Polished Correction","Polished Correction",bins,minCTR,maxCTR);

  // std::vector<TH1F*> histograms;

  // do summary canvases for checking the fits

  int sqrtCrystals = ceil(sqrt( crystal.size() ) );

  TCanvas *cSumSimple  = new TCanvas("Summary Basic CTR","Summary Basic CTR",1200,1200);
  TCanvas *cSumCentral = new TCanvas("Summary Central CTR","Summary Central CTR",1200,1200);
  TCanvas *cSumAll     = new TCanvas("Summary Full CTR","Summary Full CTR",1200,1200);
  TCanvas *cPoliAll     = new TCanvas("Summary Polished CTR","Summary Polished CTR",1200,1200);
  cSumSimple ->Divide(sqrtCrystals,sqrtCrystals);
  cSumCentral->Divide(sqrtCrystals,sqrtCrystals);
  cSumAll->Divide(sqrtCrystals,sqrtCrystals);
  cPoliAll->Divide(sqrtCrystals,sqrtCrystals);

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  for(unsigned int iCry = 0 ;  iCry < crystal.size() ; iCry++)
  {

    if(likelihood)
    {
      if(crystal[iCry].accepted)
      {
        //create the 2d matrix of tgraphs
        int nLike = crystal[iCry].slice[0].tChannel.size();
        for(int iCorr = 0; iCorr < nLike ; iCorr++)
        {
          for(int jCorr = 0; jCorr < nLike ; jCorr++)
          {
            crystal[iCry].inverse_covariance_element[iCorr][jCorr]->Write();
          }
        }
      }
    }

    std::stringstream sname;

    if(smooth)
    {
      if(crystal[iCry].simpleCTR)
      {
        crystal[iCry].simpleCTR ->Smooth(smooth);
      }
      if(crystal[iCry].centralCTR)
      {
        crystal[iCry].centralCTR ->Smooth(smooth);
      }
      if(crystal[iCry].allCTR)
      {
        crystal[iCry].allCTR ->Smooth(smooth);
      }
      if(crystal[iCry].poliCorrCTR)
      {
        crystal[iCry].poliCorrCTR ->Smooth(smooth);
      }
    }


    Float_t realBasicCTRfwhm,realBasicCTRfwtm ;
    Float_t realCentralCTRfwhm,realCentralCTRfwtm;
    Float_t realAllCTRfwhm,realAllCTRfwtm;
    Float_t poliCorrCTRfwhm,poliCorrCTRfwtm;
    Float_t reallikeCTRfwhm,reallikeCTRfwtm;
    Float_t realhybridCTRfwhm,realhybridCTRfwtm;
    double unbinnedSimpleCTR;
    double unbinnedCentralCTR;
    double unbinnedAllCTR;
    double unbinnedPoliCTR;
    double ret[2];
    double fitRes[3];
    Int_t CTRentries;
    Float_t lightCentral;
    Float_t lightAll;

    // get data on entries and light collected
    // CTR entries
    CTRentries = crystal[iCry].basicCTRhisto->GetEntries();
    crystal[iCry].basicCTRhisto->Write();

    // light central
    TF1 *gaussCentral = new TF1("gaussCentral","gaus");
    crystal[iCry].lightCentralHisto->Fit(gaussCentral,"Q");
    lightCentral = gaussCentral->GetParameter(1);
    crystal[iCry].lightCentralHisto->Write();

    // light central
    TF1 *gaussAll = new TF1("gaussAll","gaus");
    crystal[iCry].lightAllHisto->Fit(gaussAll,"Q");
    lightAll = gaussAll->GetParameter(1);
    crystal[iCry].lightAllHisto->Write();

    // std::cout << "BASIC CTRs --------------------" << std::endl;
    // std::cout << crystal[iCry]
    if(crystal[iCry].simpleCTR)
    {

      crystal[iCry].simpleCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].simpleCTR->SetFillStyle(3001);
      crystal[iCry].simpleCTR->SetFillColor(kGreen);
      crystal[iCry].simpleCTR->SetLineColor(kGreen);
      crystal[iCry].simpleCTR->SetStats(0);
      crystal[iCry].simpleCTR_norm = (TH1F*) crystal[iCry].simpleCTR->Clone();
      if(func == 0)
      {
        extractCTR(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].simpleCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].simpleCTR,0.68,true,tagFwhm);
        }
      }

      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;


      //
      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      textfile  << "No corr    - cry " << crystal[iCry].number << "\t"
                << ret[0]*1e12 << "\t"
                << ret[1]*1e12 << "\t"
                << CTRentries << "\t"
                << lightCentral << "\t"
                << lightAll << "\t"
                << fitRes[0] << "\t"
                << fitRes[1] << "\t"
                << fitRes[2] << "\t"
                << std::endl;

      realBasicCTRfwhm = ret[0]*1e12;
      realBasicCTRfwtm = ret[1]*1e12;
      noCorr->Fill(ret[0]*1e12);

      crystal[iCry].simpleCTR->Write();
      cSumSimple->cd(iCry+1);
      crystal[iCry].simpleCTR->Draw();
      crystal[iCry].simpleCTR_norm->Scale(1.0/crystal[iCry].simpleCTR_norm->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vSimple,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedSimpleCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnednoCorr->Fill(unbinnedSimpleCTR);

        std::cout << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
                  << unbinnedSimpleCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned No corr    - cry " << crystal[iCry].number << "\t"
                  << unbinnedSimpleCTR << "\t"
                  << 0 << std::endl;
      }
    }

    if(crystal[iCry].centralCTR)
    {
      crystal[iCry].centralCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].centralCTR->SetFillStyle(3001);
      crystal[iCry].centralCTR->SetFillColor(kBlue);
      crystal[iCry].centralCTR->SetLineColor(kBlue);
      crystal[iCry].centralCTR->SetStats(0);
      crystal[iCry].centralCTR_norm = (TH1F*) crystal[iCry].centralCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].centralCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].centralCTR,0.68,true,tagFwhm);
        }
      }


      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Central    - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realCentralCTRfwhm = ret[0]*1e12;
      realCentralCTRfwtm = ret[1]*1e12;
      centralCorr->Fill(ret[0]*1e12);

      crystal[iCry].centralCTR->Write();
      cSumCentral->cd(iCry+1);
      crystal[iCry].centralCTR->Draw();
      // crystal[iCry].centralCTR->Scale(1.0/crystal[iCry].centralCTR->GetMaximum());
      crystal[iCry].centralCTR_norm->Scale(1.0/crystal[iCry].centralCTR_norm->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vCentral,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedCentralCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedcentralCorr->Fill(unbinnedCentralCTR);

        std::cout << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
                  << unbinnedCentralCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Central    - cry " << crystal[iCry].number << "\t"
                  << unbinnedCentralCTR << "\t"
                  << 0 << std::endl;
      }
    }

    if(crystal[iCry].allCTR)
    {
      crystal[iCry].allCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].allCTR->SetFillStyle(3001);
      crystal[iCry].allCTR->SetFillColor(kRed);
      crystal[iCry].allCTR->SetLineColor(kRed);
      crystal[iCry].allCTR->SetStats(0);

      crystal[iCry].allCTR_norm = (TH1F*) crystal[iCry].allCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].allCTR,0.68,true,tagFwhm);
        }
      }
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      std::cout << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      realAllCTRfwhm = ret[0]*1e12;
      realAllCTRfwtm = ret[1]*1e12;
      fullCorr->Fill(ret[0]*1e12);
      crystal[iCry].allCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].allCTR->Draw();
      crystal[iCry].allCTR_norm->Scale(1.0/crystal[iCry].allCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vAll,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedfullCorr->Fill(unbinnedAllCTR);

        std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedAllCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedAllCTR << "\t"
                  << 0 << std::endl;
      }
    }





    if(crystal[iCry].likeCTR)
    {
      crystal[iCry].likeCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].likeCTR->SetFillStyle(3001);
      crystal[iCry].likeCTR->SetFillColor(kRed);
      crystal[iCry].likeCTR->SetLineColor(kRed);
      crystal[iCry].likeCTR->SetStats(0);

      crystal[iCry].likeCTR_norm = (TH1F*) crystal[iCry].likeCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].likeCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].likeCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].likeCTR,0.68,true,tagFwhm);
        }
      }



      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Like corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;


      reallikeCTRfwhm = ret[0]*1e12;
      reallikeCTRfwtm = ret[1]*1e12;
      likeCorr->Fill(ret[0]*1e12);
      crystal[iCry].likeCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].likeCTR->Draw();
      crystal[iCry].likeCTR_norm->Scale(1.0/crystal[iCry].likeCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vAll,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   unbinnedfullCorr->Fill(unbinnedAllCTR);
      //
      //   std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      // }
    }



    if(crystal[iCry].hybridCTR)
    {
      crystal[iCry].hybridCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].hybridCTR->SetFillStyle(3001);
      crystal[iCry].hybridCTR->SetFillColor(kRed);
      crystal[iCry].hybridCTR->SetLineColor(kRed);
      crystal[iCry].hybridCTR->SetStats(0);

      crystal[iCry].hybridCTR_norm = (TH1F*) crystal[iCry].hybridCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].allCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].hybridCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].hybridCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].hybridCTR,0.68,true,tagFwhm);
        }
      }



      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Hybr corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Full corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;


      realhybridCTRfwhm = ret[0]*1e12;
      realhybridCTRfwtm = ret[1]*1e12;
      hybridCorr->Fill(ret[0]*1e12);
      crystal[iCry].hybridCTR->Write();
      cSumAll->cd(iCry+1);
      crystal[iCry].hybridCTR->Draw();
      crystal[iCry].hybridCTR_norm->Scale(1.0/crystal[iCry].hybridCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());

      // use unbinned method
      // if(unbinned)
      // {
      //   double mean,meanErr,min,max;
      //   double delta = FindSmallestInterval(mean,
      //                                       meanErr,
      //                                       min,
      //                                       max,
      //                                       crystal[iCry].vAll,
      //                                       0.68,
      //                                       true);
      //   //now pass to fwhm
      //   double fwhm = 2.355 * (delta/2.0);
      //   unbinnedAllCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
      //   unbinnedfullCorr->Fill(unbinnedAllCTR);
      //
      //   std::cout << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      //
      //   textfile  << "Unbinned Full corr. - cry " << crystal[iCry].number << "\t"
      //             << unbinnedAllCTR << "\t"
      //             << 0 << std::endl;
      // }
    }








    if(crystal[iCry].poliCorrCTR)
    {
      crystal[iCry].poliCorrCTR->GetXaxis()->SetTitle("Time [s]");
      crystal[iCry].poliCorrCTR->SetFillStyle(3001);
      crystal[iCry].poliCorrCTR->SetFillColor(kBlack);
      crystal[iCry].poliCorrCTR->SetLineColor(kBlack);
      crystal[iCry].poliCorrCTR->SetStats(0);

      crystal[iCry].poliCorrCTR_norm = (TH1F*) crystal[iCry].poliCorrCTR->Clone();
      // if(func == 0)
      // {
      //   extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      // else
      // {
      //   extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
      // }
      if(func == 0)
      {
        extractCTR(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret,fitRes);
      }
      else
      {
        if(func ==1)
        {
          extractWithGaussAndExp(crystal[iCry].poliCorrCTR,fitPercMin,fitPercMax,divs,tagFwhm,ret);
        }
        else
        { //(double* retValues, TH1F* histo, const float& fraction, const bool& verbosity)
          FindSmallestInterval(ret,crystal[iCry].poliCorrCTR,0.68,true,tagFwhm);
        }
      }


      // std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // std::cout << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;
      //
      // textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      //           << ret[0]*1e12 << "\t"
      //           << ret[1]*1e12 << std::endl;
      // //
      // textfile << "Polished FIT   - cry " << crystal[iCry].number << "\t"
      //           << fitRes[0] << "\t"
      //           << fitRes[1] << "\t"
      //           << fitRes[2] << "\t"
      //           << std::endl;


      //
      std::cout << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;


      std::cout << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      textfile << "# Condition" << "\t"
                << "CTRfwhm"<< "\t"
                << "CTRfwtm"<< "\t"
                << "CTRentries"<< "\t"
                << "lightCentral"<< "\t"
                << "lightAll"<< "\t"
                << "ChiSquare"<< "\t"
                << "NDF"<< "\t"
                << "Prob"<< "\t"
                << std::endl;

      textfile  << "Polished corr. - cry " << crystal[iCry].number << "\t"
      << ret[0]*1e12 << "\t"
      << ret[1]*1e12 << "\t"
      << CTRentries << "\t"
      << lightCentral << "\t"
      << lightAll << "\t"
      << fitRes[0] << "\t"
      << fitRes[1] << "\t"
      << fitRes[2] << "\t"
      << std::endl;

      poliCorrCTRfwhm = ret[0]*1e12;
      poliCorrCTRfwtm = ret[1]*1e12;
      poliCorr->Fill(ret[0]*1e12);
      crystal[iCry].poliCorrCTR->Write();
      cPoliAll->cd(iCry+1);
      crystal[iCry].poliCorrCTR->Draw();
      crystal[iCry].poliCorrCTR_norm->Scale(1.0/crystal[iCry].poliCorrCTR_norm->GetMaximum());
      // crystal[iCry].allCTR->Scale(1.0/crystal[iCry].allCTR->GetMaximum());


      // use unbinned method
      if(unbinned)
      {
        double mean,meanErr,min,max;
        double delta = FindSmallestInterval(mean,
                                            meanErr,
                                            min,
                                            max,
                                            crystal[iCry].vPoli,
                                            0.68,
                                            true);
        //now pass to fwhm
        double fwhm = 2.355 * (delta/2.0);
        unbinnedPoliCTR = 1e12*sqrt(2)*sqrt(pow(fwhm,2)-pow(tagFwhm,2));
        unbinnedpoliCorr->Fill(unbinnedPoliCTR);

        std::cout << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedPoliCTR << "\t"
                  << 0 << std::endl;

        textfile  << "Unbinned Polished corr. - cry " << crystal[iCry].number << "\t"
                  << unbinnedPoliCTR << "\t"
                  << 0 << std::endl;
      }
    }



    sname.str("");

    sname << "Summary - Crystal " << crystal[iCry].number;
    TCanvas* c_summary = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
    c_summary->cd();
    THStack *hs = new THStack("hs","");
    hs->Add(crystal[iCry].simpleCTR_norm);
    hs->Add(crystal[iCry].centralCTR_norm);
    hs->Add(crystal[iCry].allCTR_norm);
    hs->Add(crystal[iCry].likeCTR_norm);
    hs->Add(crystal[iCry].poliCorrCTR_norm);
    hs->Add(crystal[iCry].hybridCTR_norm);


    // std::cout << "Crystal " << crystal[iCry].number << std::endl;
    // std::cout << "CTR FWHM [ps] " << std::endl;
    hs->Draw("hist nostack");
    sname.str("");
    sname << "CTR - Crystal " << crystal[iCry].number << " - width in FWHM";
    hs->SetTitle(sname.str().c_str());
    hs->GetXaxis()->SetTitle("Time [s]");
    hs->GetXaxis()->SetTitleOffset(1);
    hs->GetXaxis()->SetTitleSize(0.045);
    hs->GetXaxis()->SetLabelSize(0.045);
    hs->GetYaxis()->SetLabelSize(0.045);
    TLegend *legend = new TLegend(0.54,0.62,0.89,0.89,"");
    legend->SetFillStyle(0);
    if(crystal[iCry].simpleCTR)
    {
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].simpleCTR,sname.str().c_str(),"f");
      // std::cout << "No correction       = "<< realBasicCTRfwhm   << " ps" << std::endl;
    }
    if(crystal[iCry].centralCTR)
    {
      sname.str("");
      sname << "Central correction = " << realCentralCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].centralCTR,sname.str().c_str(),"f");
      // std::cout << "Central correction  = "<< realCentralCTRfwhm << " ps" << std::endl;
    }
    if(crystal[iCry].allCTR)
    {
      sname.str("");
      sname << "Full correction       = " << realAllCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].allCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].likeCTR)
    {
      sname.str("");
      sname << "Likelihood correction = " << reallikeCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].likeCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].poliCorrCTR)
    {
      sname.str("");
      sname << "Polished correction       = " << poliCorrCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].poliCorrCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }
    if(crystal[iCry].hybridCTR)
    {
      sname.str("");
      sname << "Polished correction       = " << realhybridCTRfwhm << "ps";
      legend->AddEntry(crystal[iCry].hybridCTR,sname.str().c_str(),"f");
      // std::cout << "Full correction     = "<< realAllCTRfwhm     << " ps" << std::endl;
    }


    sname.str("");
    legend->Draw();
    gStyle->SetOptTitle(0);
    TPaveLabel *title = new TPaveLabel(.11,.95,.35,.99,"new title","brndc");
    title->Draw();
    // std::cout << std::endl;

    c_summary->Write();


    TH1F* cloneBasic;
    TH1F* cloneCentral;
    TH1F* cloneAll;
    TH1F* cloneLike;
    TH1F* clonePoli;
    TH1F* cloneHybrid;
    THStack *cloneHs = (THStack*) hs->Clone();
    TLegend *legend1 = new TLegend(0.15,0.69,0.49,0.89,"");
    legend1->SetFillStyle(0);
    sname.str("");
    sname << "Multi - Crystal " << crystal[iCry].number;
    TCanvas* c_multi = new TCanvas(sname.str().c_str(),sname.str().c_str(),1800,1400);
    c_multi->Divide(2,2);

    if(crystal[iCry].simpleCTR_norm)
    {
      cloneBasic   = (TH1F*) crystal[iCry].simpleCTR->Clone();
      c_multi->cd(1);
      sname.str("");
      sname << "No correction        = " << realBasicCTRfwhm << "ps";
      legend1->AddEntry(cloneBasic,sname.str().c_str(),"f");
      cloneBasic->Draw();
      legend1->Draw();
    }
    if(crystal[iCry].centralCTR_norm)
    {
      cloneCentral = (TH1F*) crystal[iCry].centralCTR->Clone();
      c_multi->cd(2);
      TLegend *legend2 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend2->SetFillStyle(0);
      sname.str("");
      sname << "Central correction   = " << realCentralCTRfwhm << "ps";
      legend2->AddEntry(cloneCentral,sname.str().c_str(),"f");
      cloneCentral->Draw();
      legend2->Draw();
    }
    if(crystal[iCry].allCTR_norm)
    {
      cloneAll     = (TH1F*) crystal[iCry].allCTR->Clone();
      c_multi->cd(3);
      TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend3->SetFillStyle(0);
      sname.str("");
      sname << "Full correction      = " << realAllCTRfwhm << "ps";
      legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
      cloneAll->Draw();
      legend3->Draw();
    }
    if(crystal[iCry].likeCTR_norm)
    {
      cloneLike     = (TH1F*) crystal[iCry].likeCTR->Clone();
      c_multi->cd(3);
      TLegend *legend4 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend4->SetFillStyle(0);
      sname.str("");
      sname << "Likelihood correction = " << reallikeCTRfwhm << "ps";
      legend4->AddEntry(cloneLike,sname.str().c_str(),"f");
      cloneLike->Draw();
      legend4->Draw();
    }
    if(crystal[iCry].hybridCTR_norm)
    {
      cloneHybrid     = (TH1F*) crystal[iCry].hybridCTR->Clone();
      c_multi->cd(3);
      TLegend *legend5 = new TLegend(0.15,0.69,0.49,0.89,"");
      legend5->SetFillStyle(0);
      sname.str("");
      sname << "Hybrid correction = " << realhybridCTRfwhm << "ps";
      legend5->AddEntry(cloneHybrid,sname.str().c_str(),"f");
      cloneHybrid->Draw();
      legend5->Draw();
    }


    // if(crystal[iCry].poliCorrCTR_norm)
    // {
    //   clonePoli     = (TH1F*) crystal[iCry].poliCorrCTR->Clone();
    //   c_multi->cd(3);
    //   TLegend *legend3 = new TLegend(0.15,0.69,0.49,0.89,"");
    //   legend3->SetFillStyle(0);
    //   sname.str("");
    //   sname << "Polished correction      = " << poliCorrCTRfwhm << "ps";
    //   legend3->AddEntry(cloneAll,sname.str().c_str(),"f");
    //   cloneAll->Draw();
    //   legend3->Draw();
    // }


    c_multi->cd(4);
    c_multi->cd(4)->SetGrid();
    cloneHs->Draw("hist nostack");
    c_multi->Write();

  }
  noCorr->Write();
  centralCorr->Write();
  fullCorr->Write();
  poliCorr->Write();
  likeCorr->Write();
  hybridCorr->Write();

  unbinnednoCorr->Write();
  unbinnedcentralCorr->Write();
  unbinnedfullCorr->Write();
  unbinnedpoliCorr->Write();

  cSumSimple ->Write();
  cSumCentral->Write();
  cSumAll->Write();
  cPoliAll->Write();


  TNamed CommandNameD("Command",streamCommand.str().c_str());
  CommandNameD.Write();
  // treeFile->Close();

  calibrationFile->Close();
  outputFile->Close();
  textfile.close();
  std::cout << std::endl;
  std::cout << "Histograms saved in   " << outputFileName << std::endl;
  std::cout << "Text summary saved in " << textFileName << std::endl;
  return 0;
}