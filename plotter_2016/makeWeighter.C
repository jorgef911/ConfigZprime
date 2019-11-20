void makeWeighter() {
  TFile* Datafile = new TFile("Data.root");
  TFile* MCfile = new TFile("MC.root");
  TFile* extrafile = new TFile("extra.root");

  Datafile->cd("METCut");
  TH1D* Data = (TH1D*)gDirectory->FindObjectAny("NVertices");

  MCfile->cd("METCut");
  TH1D* MC = (TH1D*)gDirectory->FindObjectAny("NVertices");

  double allData = Data->Integral();
  double allMC = MC->Integral();
  double factor = allMC/allData;

  TFile* newfile = new TFile("newWeight.root", "RECREATE");
  newfile->cd();
  TH1D* weight = new TH1D("extra", "extra", 100, 0, 100);
  for(int i = 0; i < weight->GetXaxis()->GetNbins(); ++i) {
    if(MC->GetBinContent(i+1) == 0) continue;
    double tmp = Data->GetBinContent(i+1)*factor/MC->GetBinContent(i+1);
    weight->SetBinContent(i+1, tmp);
    cout <<  (i+1) << " | "  << tmp << " | " << 
      tmp*MC->GetBinContent(i+1) << " | " << ((tmp!=0) ? Data->GetBinContent(i+1)/(tmp*MC->GetBinContent(i+1)) : 0) <<  endl;
  }
  weight->Write();
  newfile->Close();
}
