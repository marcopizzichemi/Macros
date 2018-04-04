void timePlots(std::string fileName, int nmppcx = 4, int nmppcy = 4, int ncryx = 2, int ncryy = 2)
{
  std::ifstream inputFile;
  inputFile.open(fileName.c_str(), std::ios::in);
  int nx = nmppcx*ncryx;
  int ny = nmppcy*ncryy;
  TH2F* histo2d = new TH2F("histo2","CTR FWHM vs. i,j",nx,0,nx,ny,0,ny);
  histo2d->GetXaxis()->SetTitle("i");
  histo2d->GetYaxis()->SetTitle("j");
  histo2d->GetZaxis()->SetTitle("[ps]");
  histo2d->GetZaxis()->SetRangeUser(150,450);
  TH1F* histo_fwhm = new TH1F("histo_fwhm","CTR FWHM",50,150,450);
  histo_fwhm->GetXaxis()->SetTitle("[ps]");
  TH1F* histo_fwhm_central = new TH1F("histo_fwhm_central","CTR FWHM central channels",50,150,450);
  histo_fwhm_central->GetXaxis()->SetTitle("[ps]");
  int num;
  double fwhm,fwtm;
  if (inputFile.is_open())
  {
    while (!inputFile.eof())
    {
      inputFile >> num >> fwhm >> fwtm;
      if(!inputFile.eof())
      {
        histo_fwhm->Fill(fwhm);
        int i = num / nx;
        int j = num % nx;
        if( (i > (ncryx-1)) && (i < (nx - (ncryx))) && (j > (ncryy-1)) && (j < (ny -(ncryy))) )
        {
          histo_fwhm_central->Fill(fwhm);
        }
        histo2d->Fill(i,j,fwhm);
      }
    }
  }
  TCanvas *c1 = new TCanvas("c1","c1",1200,800);
  histo_fwhm->Draw();
  TCanvas *c1_2 = new TCanvas("c1_2","c1_2",1200,800);
  histo_fwhm_central->Draw();
  TCanvas *c2 = new TCanvas("c2","c2",800,800);
  histo2d->Draw("LEGO2");

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
  outputFile->Close();
}
