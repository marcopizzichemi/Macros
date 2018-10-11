struct crystal_t
{
  int num;
  float fwhm,peakNoSat,sumNeighNoSat,sumAllNoSat,peakSat,sumNeighSat,sumAllSat;
};


void timePlots_v2(std::string fileName, int nmppcx = 4, int nmppcy = 4, int ncryx = 2, int ncryy = 2)
{
  std::ifstream inputFile;
  inputFile.open(fileName.c_str(), std::ios::in);
  int nx = nmppcx*ncryx;
  int ny = nmppcy*ncryy;

  int num;
  std::vector<crystal_t> crystal;
  if (inputFile.is_open())
  {
    while (!inputFile.eof())
    {
      int num;
      float fwhm,peakNoSat,sumNeighNoSat,sumAllNoSat,peakSat,sumNeighSat,sumAllSat;
      inputFile >> num >> fwhm >> peakNoSat >> sumNeighNoSat >> sumAllNoSat >> peakSat >> sumNeighSat >> sumAllSat;
      if(!inputFile.eof())
      {
        crystal_t temp_cry;
        temp_cry.num = num;
        temp_cry.fwhm = fwhm;
        temp_cry.peakNoSat = peakNoSat;
        temp_cry.sumNeighNoSat = sumNeighNoSat;
        temp_cry.sumAllNoSat = sumAllNoSat;
        temp_cry.peakSat = peakSat;
        temp_cry.sumNeighSat = sumNeighSat;
        temp_cry.sumAllSat = sumAllSat;

        crystal.push_back(temp_cry);
        // histo_fwhm->Fill(fwhm);
        // int i = num / nx;
        // int j = num % nx;
        // if( (i > (ncryx-1)) && (i < (nx - (ncryx))) && (j > (ncryy-1)) && (j < (ny -(ncryy))) )
        // {
        //   histo_fwhm_central->Fill(fwhm);
        // }
        // histo2d->Fill(i,j,fwhm);
      }
    }
  }

  // calc some limits
  float CTRmin = +INFINITY;
  float CTRmax = -INFINITY;
  float peakNoSatmin = +INFINITY;
  float peakNoSatmax = -INFINITY;

  for(unsigned int iCry = 0 ; iCry < crystal.size(); iCry++)
  {
    if(crystal[iCry].fwhm < CTRmin)
    {
      CTRmin = crystal[iCry].fwhm ;
    }
    if(crystal[iCry].fwhm > CTRmax)
    {
      CTRmax = crystal[iCry].fwhm ;
    }
    if(crystal[iCry].peakNoSat < peakNoSatmin)
    {
      peakNoSatmin = crystal[iCry].peakNoSat ;
    }
    if(crystal[iCry].peakNoSat > peakNoSatmax)
    {
      peakNoSatmax = crystal[iCry].peakNoSat ;
    }
  }

  CTRmin = CTRmin - (CTRmax-CTRmin)/2.0;
  CTRmax = CTRmax + (CTRmax-CTRmin)/2.0;
  peakNoSatmin = peakNoSatmin - (peakNoSatmax-peakNoSatmin)/2.0;
  peakNoSatmax = peakNoSatmax + (peakNoSatmax-peakNoSatmin)/2.0;

  TH2F* histo2d = new TH2F("histo2","CTR FWHM vs. i,j",nx,0,nx,ny,0,ny);
  histo2d->GetXaxis()->SetTitle("i");
  histo2d->GetYaxis()->SetTitle("j");
  histo2d->GetZaxis()->SetTitle("[ps]");
  histo2d->GetZaxis()->SetRangeUser(CTRmin,CTRmax);
  TH1F* histo_fwhm = new TH1F("histo_fwhm","CTR FWHM",50,CTRmin,CTRmax);
  histo_fwhm->GetXaxis()->SetTitle("[ps]");
  TH1F* histo_fwhm_central = new TH1F("histo_fwhm_central","CTR FWHM central channels",50,CTRmin,CTRmax);
  histo_fwhm_central->GetXaxis()->SetTitle("[ps]");

  TH2F* histoCTRvsPixelsFired = new TH2F("histoCTRvsPixelsFired","CTR FWHM vs. pixels fired trigger detector",50,peakNoSatmin,peakNoSatmax,50,CTRmin,CTRmax);

  for(unsigned int iCry = 0 ; iCry < crystal.size(); iCry++)
  {
    histo_fwhm->Fill(crystal[iCry].fwhm);
    int i = crystal[iCry].num / nx;
    int j = crystal[iCry].num % nx;
    if( (i > (ncryx-1)) && (i < (nx - (ncryx))) && (j > (ncryy-1)) && (j < (ny -(ncryy))) )
    {
      histo_fwhm_central->Fill(crystal[iCry].fwhm);
    }
    histo2d->Fill(i,j,crystal[iCry].fwhm);
    histoCTRvsPixelsFired->Fill(crystal[iCry].peakNoSat,crystal[iCry].fwhm);
  }




  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  histo_fwhm->Draw();
  TCanvas *c1_2 = new TCanvas("c1_2","c1_2",1200,800);
  histo_fwhm_central->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  histo2d->Draw("LEGO2");
  TCanvas *c3 = new TCanvas("c3","c3",800,800);
  histoCTRvsPixelsFired->Draw("COLZ");

  //output averages
  std::cout << "Summary " << std::endl;
  std::cout << "CTR FWHM [ps]" << std::endl;
  std::cout << "All crystals " << "(" << histo_fwhm->GetEntries()<< ")" << std::endl;
  std::cout << "Central crystals " << "(" << histo_fwhm_central->GetEntries()<< ")" << std::endl;
  std::cout << histo_fwhm->GetMean() << " +/- " << histo_fwhm->GetRMS() / sqrt(histo_fwhm->GetEntries())  << std::endl;
  std::cout << histo_fwhm_central->GetMean() << " +/- " << histo_fwhm_central->GetRMS() / sqrt(histo_fwhm_central->GetEntries())  << std::endl;

  //save plots to file
  std::string outputFileName = fileName.substr(0,fileName.size()-4);
  outputFileName += ".root";
  TFile *outputFile = new TFile(outputFileName.c_str(),"RECREATE");
  outputFile->cd();
  c1->Write();
  c1_2->Write();
  c2->Write();
  histo_fwhm->Write();
  histo_fwhm_central->Write();
  histo2d->Write();
  histoCTRvsPixelsFired->Write();
  outputFile->Close();
}
