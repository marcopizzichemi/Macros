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


void thScan(std::string inputFileName = "TTree_", std::string outputFileName = "plots.root" ,int chA = 2, int chB = 9, double sigmas = 2.0, double percDown = 0.04, double percUp = 0.06, double ENres = 0.10)
{
  gStyle->SetOptFit(1111111);
  std::vector<std::string> v;
  read_directory(".", v);
  // extract files with correct prefix
  std::vector<std::string> listInputFiles;
  // std::string inputFileName = "TTree_";
  for(unsigned int i = 0 ; i < v.size() ; i++)
  {
    if(!v[i].compare(0,inputFileName.size(),inputFileName))
    {
      listInputFiles.push_back(v[i]);
    }
  }

  TChain *tree = new TChain("adc"); // read input files
  for(unsigned int i = 0 ; i < listInputFiles.size() ; i++)
  {
    std::cout << listInputFiles[i] << std::endl;
    tree->Add(listInputFiles[i].c_str());
  }

  std::stringstream sname;
  sname << "t_ch" << chA;
  TH1F *t_chA = new TH1F(sname.str().c_str(),sname.str().c_str(),1000,-0.15e-6,0.05e-6);
  sname.str("");

  sname << "c_" << chA;
  TCanvas *c_tA = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
  sname.str("");
  sname << "t" << chA << " >> " << "t_ch" << chA;
  tree->Draw(sname.str().c_str());
  sname.str("");

  sname << "t_ch" << chB;
  TH1F *t_chB = new TH1F(sname.str().c_str(),sname.str().c_str(),1000,-0.15e-6,0.05e-6);
  sname.str("");

  sname << "c_" << chB;
  TCanvas *c_tB = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
  sname.str("");
  sname << "t" << chB << " >> " << "t_ch" << chB;
  tree->Draw(sname.str().c_str());
  sname.str("");

  // TH1F *t_ch9 = new TH1F("t_ch9","t_ch9",1000,-0.15e-6,0.05e-6);
  // TCanvas *c_t9 = new TCanvas("c_t9","c_t9",1200,800);
  // tree->Draw("t9 >> t_ch9");


  sname << "h_ch" << chA;
  TH1F *h_chA = new TH1F(sname.str().c_str(),sname.str().c_str(),1000,0,40000);
  sname.str("");

  sname << "c_h" << chA;
  TCanvas *c_hA = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
  sname.str("");
  sname << "ch" << chA << " >> " << "h_ch" << chA;
  tree->Draw(sname.str().c_str());
  sname.str("");

  sname << "h_ch" << chB;
  TH1F *h_chB = new TH1F(sname.str().c_str(),sname.str().c_str(),1000,0,40000);
  sname.str("");

  sname << "c_h" << chB;
  TCanvas *c_hB = new TCanvas(sname.str().c_str(),sname.str().c_str(),1200,800);
  sname.str("");
  sname << "ch" << chB << " >> " << "h_ch" << chB;
  tree->Draw(sname.str().c_str());
  sname.str("");

  // find cuts
  //ch2
  c_hA->cd();
  TSpectrum *spectrumCH2;
  spectrumCH2 = new TSpectrum(5);
  Int_t peaksN = spectrumCH2->Search(h_chA,1,"",0.5);
  Double_t *PeaksCH2  = spectrumCH2->GetPositionX();
  Double_t *PeaksYCH2 = spectrumCH2->GetPositionY();
  float maxPeak = 0.0;
  int peakID = 0;
  for (int peakCounter = 0 ; peakCounter < peaksN ; peakCounter++ )
  {
    if(PeaksYCH2[peakCounter] > maxPeak)
    {
      maxPeak = PeaksYCH2[peakCounter];
      peakID = peakCounter;
    }
  }
  TF1 *gaussCH2 = new TF1("gaussCH2","gaus");
  gaussCH2->SetParameter(0,PeaksYCH2[peakID]); // heigth
  gaussCH2->SetParameter(1,PeaksCH2[peakID]); //mean
  gaussCH2->SetParameter(2,(ENres*PeaksCH2[peakID])/2.355); // sigma
  h_chA->Fit(gaussCH2,"","",PeaksCH2[peakID] - (PeaksCH2[peakID] * percDown),PeaksCH2[peakID] + (PeaksCH2[peakID] * percUp));

  sname << "ch" << chA << " > " << PeaksCH2[peakID] - (sigmas * gaussCH2->GetParameter(2)) << "&& ch" << chA << " < " << PeaksCH2[peakID] + (sigmas * gaussCH2->GetParameter(2)) ;
  TCut peakch2 = sname.str().c_str();
  sname.str("");

  //ch9
  c_hB->cd();
  TSpectrum *spectrumCH9;
  spectrumCH9 = new TSpectrum(5);
  Int_t peaksN9 = spectrumCH9->Search(h_chB,1,"",0.5);
  Double_t *PeaksCH9  = spectrumCH9->GetPositionX();
  Double_t *PeaksYCH9 = spectrumCH9->GetPositionY();
  float maxPeak9 = 0.0;
  int peakID9 = 0;
  for (int peakCounter = 0 ; peakCounter < peaksN9 ; peakCounter++ )
  {
    if(PeaksYCH9[peakCounter] > maxPeak9)
    {
      maxPeak9 = PeaksYCH9[peakCounter];
      peakID9 = peakCounter;
    }
  }
  TF1 *gaussCH9 = new TF1("gaussCH9","gaus");
  gaussCH9->SetParameter(0,PeaksYCH9[peakID9]); // heigth
  gaussCH9->SetParameter(1,PeaksCH9[peakID9]); //mean
  gaussCH9->SetParameter(2,(ENres*PeaksCH9[peakID9])/2.355); // sigma
  h_chB->Fit(gaussCH9,"","",PeaksCH9[peakID9] - (PeaksCH9[peakID9] * percDown),PeaksCH9[peakID9] + (PeaksCH9[peakID9] * percUp));
  // std::stringstream sname;
  sname << "ch" << chB << " > " << PeaksCH9[peakID9] - (sigmas * gaussCH9->GetParameter(2)) << "&& ch" << chB << " < " << PeaksCH9[peakID9] + (sigmas * gaussCH9->GetParameter(2)) ;
  TCut peakch9 = sname.str().c_str();
  sname.str("");

  // TCut peakch9 = "ch9 > 19000 && ch9 < 20500";

  //avoid zeroes!
  sname << "t" << chA << " != 0 && t" << chB <<" != 0";
  TCut noZeros = sname.str().c_str();
  sname.str("");

  TH1F *ctr = new TH1F("ctr","ctr",500,-1e-9,1e-9);
  TCanvas *c_ht = new TCanvas("c_ht","c_ht",1200,800);
  sname << "t" << chB << " - t" << chA <<" >> ctr";
  tree->Draw(sname.str().c_str(),peakch2+peakch9+noZeros);
  sname.str("");

  TF1 *gaussCTR = new TF1("gaussCTR","gaus");
  ctr->Fit(gaussCTR);
  std::cout << "CTR FWHM = "<< gaussCTR->GetParameter(2)*2.355*1e12 << " ps" << std::endl;

  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  h_chA->Write();
  h_chB->Write();
  t_chA->Write();
  t_chB->Write();
  ctr->Write();
  outputFile->Close();
}
