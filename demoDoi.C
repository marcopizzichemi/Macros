// macro to extract demonstration of difference of DOI (+ amplitude walk) impact and amplitude walk impact alone
//procedure can be found in amplitudeWalk.cpp source file


void demoDoi(std::string fileName,double min = 0,double max = 3, double length = 15)
{
  gStyle->SetOptStat(0);
  TFile *inputFile = new TFile(fileName.c_str());
  inputFile->cd();

  // with saturation correction
  TH2F* histoADCvsDt = (TH2F*) gDirectory->Get("histoADCvsDt");
  std::stringstream sname;
  sname << "histoSat_" << min << "_" << max;
  TH2F* histo1_5     = (TH2F*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  TH2F* clone_histoADCvsDt = (TH2F*) histoADCvsDt->Clone();
  TH2F* clone_histo1_5     = (TH2F*) histo1_5->Clone();

  clone_histoADCvsDt->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histo1_5    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histoADCvsDt->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histo1_5    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");

  clone_histoADCvsDt->SetTitle("All crystal -> DOI + Amplitude effect");

  sname << "DOI "<< length-min << "mm to "<< length-max<< "mm - Only Amplitude effect";
  clone_histo1_5    ->SetTitle(sname.str().c_str());
  sname.str("");

  clone_histoADCvsDt->Rebin2D(10,10);
  clone_histo1_5->Rebin2D(10,10);

  histoADCvsDt->Rebin2D(30,30);
  histo1_5->Rebin2D(30,30);
  histoADCvsDt->ProfileX();
  histo1_5->ProfileX();

  TProfile *histoADCvsDt_pfx = (TProfile*) gDirectory->Get("histoADCvsDt_pfx");

  sname << "histoSat_" << min << "_" << max << "_pfx";
  TProfile *histo1_5_pfx = (TProfile*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  TProfile* clone_histoADCvsDt_pfx = (TProfile*) histoADCvsDt_pfx->Clone();
  TProfile* clone_histo1_5_pfx     = (TProfile*) histo1_5_pfx->Clone();


  clone_histoADCvsDt_pfx->GetYaxis()->SetRangeUser(-8.8e-9,-7.5e-9);
  clone_histo1_5_pfx    ->GetYaxis()->SetRangeUser(-8.8e-9,-7.5e-9);
  clone_histoADCvsDt_pfx->SetLineColor(kBlue);
  clone_histo1_5_pfx    ->SetLineColor(kRed);
  clone_histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  clone_histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  clone_histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  TLegend *legend1 = new TLegend(0.54,0.62,0.89,0.89,"");
  legend1->SetFillStyle(0);
  legend1->AddEntry(clone_histoADCvsDt_pfx,"DOI + Amplitude","l");
  legend1->AddEntry(clone_histo1_5_pfx,"Only Amplitude","l");

  TLegend *legend2 = new TLegend(0.54,0.62,0.89,0.89,"");
  legend2->SetFillStyle(0);
  legend2->AddEntry(histoADCvsDt_pfx,"DOI + Amplitude","l");
  legend2->AddEntry(histo1_5_pfx,"Only Amplitude","l");

  histoADCvsDt_pfx->GetYaxis()->SetRangeUser(-8.8e-9,-8.3e-9);
  histo1_5_pfx    ->GetYaxis()->SetRangeUser(-8.8e-9,-8.3e-9);
  histoADCvsDt_pfx->SetLineColor(kBlue);
  histo1_5_pfx    ->SetLineColor(kRed);
  histoADCvsDt_pfx->GetXaxis()->SetRangeUser(50000,95000);
  histo1_5_pfx->GetXaxis()->SetRangeUser(50000,95000);
  histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  // fitting
  // TF1 *pol1_doi = new TF1("pol1_doi","pol1");
  // TF1 *pol1_ampl = new TF1("pol1_ampl","pol1");
  // pol1_doi ->SetLineColor(kBlue);
  // pol1_ampl->SetLineColor(kRed);
  // histoADCvsDt_pfx->Fit(pol1_doi);
  // histo1_5_pfx->Fit(pol1_ampl);

  TCanvas* canvas = new TCanvas("canvas","canvas",1600,1300);
  canvas->Divide(2,2);
  canvas->cd(1);
  clone_histoADCvsDt->Draw("COLZ");
  canvas->cd(2);
  clone_histo1_5->Draw("COLZ");
  canvas->cd(3);
  clone_histoADCvsDt_pfx->SetTitle("All Amplitude range");
  clone_histoADCvsDt_pfx->Draw();
  clone_histo1_5_pfx->Draw("same");
  legend1->Draw("same");

  canvas->cd(4);
  histoADCvsDt_pfx->SetTitle("511 KeV Amplitude range");
  histoADCvsDt_pfx->Draw("hist");
  histo1_5_pfx->Draw("same hist");
  legend2->Draw("same");


  //without saturation correction
  TH2F* noSat_histoADCvsDt = (TH2F*) gDirectory->Get("histoADCvsDtNoSat");

  sname << "histoNoSat_" << min << "_" << max;
  TH2F* noSat_histo1_5     = (TH2F*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  TH2F* noSat_clone_histoADCvsDt = (TH2F*) noSat_histoADCvsDt->Clone();
  TH2F* noSat_clone_histo1_5     = (TH2F*) noSat_histo1_5->Clone();

  noSat_clone_histoADCvsDt->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  noSat_clone_histo1_5    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  noSat_clone_histoADCvsDt->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  noSat_clone_histo1_5    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");

  noSat_clone_histoADCvsDt->SetTitle("All crystal -> DOI + Amplitude effect");

  sname << "DOI "<< length-min << "mm to "<< length-max<< "mm - Only Amplitude effect";
  noSat_clone_histo1_5    ->SetTitle(sname.str().c_str());
  sname.str("");
  noSat_clone_histoADCvsDt->Rebin2D(10,10);
  noSat_clone_histo1_5->Rebin2D(10,10);

  noSat_histoADCvsDt->Rebin2D(30,30);
  noSat_histo1_5->Rebin2D(30,30);
  noSat_histoADCvsDt->ProfileX();
  noSat_histo1_5->ProfileX();

  TProfile *noSat_histoADCvsDt_pfx = (TProfile*) gDirectory->Get("histoADCvsDtNoSat_pfx");
  sname << "histoNoSat_" << min << "_" << max << "_pfx";
  TProfile *noSat_histo1_5_pfx = (TProfile*) gDirectory->Get(sname.str().c_str());
  sname.str("");

  TProfile* noSat_clone_histoADCvsDt_pfx = (TProfile*) noSat_histoADCvsDt_pfx->Clone();
  TProfile* noSat_clone_histo1_5_pfx     = (TProfile*) noSat_histo1_5_pfx->Clone();


  noSat_clone_histoADCvsDt_pfx->GetYaxis()->SetRangeUser(-8.8e-9,-7.5e-9);
  noSat_clone_histo1_5_pfx    ->GetYaxis()->SetRangeUser(-8.8e-9,-7.5e-9);
  noSat_clone_histoADCvsDt_pfx->SetLineColor(kBlue);
  noSat_clone_histo1_5_pfx    ->SetLineColor(kRed);
  noSat_clone_histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  noSat_clone_histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  noSat_clone_histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  noSat_clone_histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  TLegend *noSat_legend1 = new TLegend(0.54,0.62,0.89,0.89,"");
  noSat_legend1->SetFillStyle(0);
  noSat_legend1->AddEntry(clone_histoADCvsDt_pfx,"DOI + Amplitude","l");
  noSat_legend1->AddEntry(clone_histo1_5_pfx,"Only Amplitude","l");

  TLegend *noSat_legend2 = new TLegend(0.54,0.62,0.89,0.89,"");
  noSat_legend2->SetFillStyle(0);
  noSat_legend2->AddEntry(histoADCvsDt_pfx,"DOI + Amplitude","l");
  noSat_legend2->AddEntry(histo1_5_pfx,"Only Amplitude","l");

  noSat_histoADCvsDt_pfx->GetYaxis()->SetRangeUser(-8.8e-9,-8.3e-9);
  noSat_histo1_5_pfx    ->GetYaxis()->SetRangeUser(-8.8e-9,-8.3e-9);
  noSat_histoADCvsDt_pfx->SetLineColor(kBlue);
  noSat_histo1_5_pfx    ->SetLineColor(kRed);
  noSat_histoADCvsDt_pfx->GetXaxis()->SetRangeUser(35000,45000);
  noSat_histo1_5_pfx    ->GetXaxis()->SetRangeUser(35000,45000);
  noSat_histoADCvsDt_pfx->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  noSat_histo1_5_pfx    ->GetYaxis()->SetTitle("T_tag - T_mppc [s]");
  noSat_histoADCvsDt_pfx->GetXaxis()->SetTitle("Amplitude [ADC ch]");
  noSat_histo1_5_pfx    ->GetXaxis()->SetTitle("Amplitude [ADC ch]");

  // fitting
  // TF1 *pol1_doi = new TF1("pol1_doi","pol1");
  // TF1 *pol1_ampl = new TF1("pol1_ampl","pol1");
  // pol1_doi ->SetLineColor(kBlue);
  // pol1_ampl->SetLineColor(kRed);
  // histoADCvsDt_pfx->Fit(pol1_doi);
  // histo1_5_pfx->Fit(pol1_ampl);

  TCanvas* noSat_canvas = new TCanvas("canvas2","canvas2",1600,1300);
  noSat_canvas->Divide(2,2);
  noSat_canvas->cd(1);
  noSat_clone_histoADCvsDt->Draw("COLZ");
  noSat_canvas->cd(2);
  noSat_clone_histo1_5->Draw("COLZ");
  noSat_canvas->cd(3);
  noSat_clone_histoADCvsDt_pfx->SetTitle("All Amplitude range");
  noSat_clone_histoADCvsDt_pfx->Draw();
  noSat_clone_histo1_5_pfx->Draw("same");
  noSat_legend1->Draw("same");

  noSat_canvas->cd(4);
  noSat_histoADCvsDt_pfx->SetTitle("511 KeV Amplitude range");
  noSat_histoADCvsDt_pfx->Draw("hist");
  noSat_histo1_5_pfx->Draw("same hist");
  noSat_legend2->Draw("same");



}
