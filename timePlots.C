void timePlots(std::string fileName, int nx = 8, int ny = 8)
{
  std::ifstream inputFile;
  inputFile.open(fileName.c_str(), std::ios::in);

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
        if(i > 1 && i < (nx - 2) && j > 1 && j < (ny -2))
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
}
