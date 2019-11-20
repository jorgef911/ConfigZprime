void findweight() {
  TFile* Datafile = new TFile("Data.root");
  Datafile->cd("METCut");
  TH1D* Data = (TH1D*)gDirectory->FindObjectAny("NVertices");
    
  TList* lister = new TList();
  TList* histos = new TList();
  lister->Add(new TFile("files/DYJetsToLL_M-10to50_amcatnloFXFX.root"));
  lister->Add(new TFile("files/DYJetsToLL_M-50_HT0to100_amcatnloFXFX.root"));
  lister->Add(new TFile("files/DYJetsToLL_M-50_HT100to200_madgraphMLM_v0-ext1-v1.root"));
  lister->Add(new TFile("files/DYJetsToLL_M-50_HT200to400_madgraphMLM_v0-ext1-v1.root"));
  lister->Add(new TFile("files/DYJetsToLL_M-50_HT400to600_madgraphMLM_v0-ext1-v1.root"));
  lister->Add(new TFile("files/DYJetsToLL_M-50_HT600toInf_madgraphMLM_v0-ext1-v1.root"));

  double lumi = 15900;
  double values[6] = {18610*0.999226, 6025.2*0.915341, 171.59*0.998792, 52.62 *0.998419, 6.77*0.997889, 2.72*0.996243};
  double integral[6];

  int j = 0;
  TFile *curfile = (TFile*)lister->First();
  while(curfile) {
    curfile->cd("METCut");
    histos->Add((TH1D*)gDirectory->FindObjectAny("NVertices"));
    integral[j] = ((TH1D*)gDirectory->FindObjectAny("Events"))->GetBinContent(1);
    cout << integral[j] << endl;
    j++;
    curfile = (TFile*)lister->After(curfile);
  }

  TH1* curhisto = (TH1*)histos->First();
  j = 0;
  while(curhisto && j < 6) {
    values[j] *= lumi;
    j++;
    curhisto = (TH1*)histos->After(curhisto);
  }

  for(int i=0; i < Data->GetXaxis()->GetNbins(); i++) {
    double first=0, second=0, third=0;
    double abar=0, a2bar=0;
    TH1* curhisto = (TH1*)histos->First();
    j = 0;
    while(curhisto && j < 6) {
      abar += values[j]*curhisto->GetBinContent(i+1)/integral[j];
      a2bar += values[j]*pow(curhisto->GetBinContent(i+1)/integral[j],2);
      j++;
      curhisto = (TH1*)histos->After(curhisto);
    }
    
    if(a2bar <= 0) continue;
    double discrim = pow(abar+a2bar,2)-4*Data->GetBinContent(i+1)*a2bar;
    bool Dodiscrim = (discrim >=0);
    
    if(Dodiscrim) cout << ((abar+a2bar)/(2*a2bar)) << " | " << ((abar+a2bar+discrim)/(2*a2bar)) << " | " << ((abar+a2bar+discrim)/(2*a2bar)) << endl;
    else cout << ((abar+a2bar)/(2*a2bar)) << " | DID NOT PASS" << endl;
  }

  // Datafile->cd("METCut");
  // TH1D* Data = (TH1D*)gDirectory->FindObjectAny("NVertices");

  // MCfile->cd("METCut");
  // TH1D* MC = (TH1D*)gDirectory->FindObjectAny("NVertices");

  // extrafile->cd();
  // TH1D* oldextra = (TH1D*)gDirectory->FindObjectAny("extra");


  // double allData = Data->Integral();
  // double allMC = MC->Integral();
  // double factor = allMC/allData;

  // TFile* newfile = new TFile("newWeight.root", "RECREATE");
  // newfile->cd();
  // TH1D* weight = new TH1D("extra", "extra", 100, 0, 100);
  // for(int i = 0; i < weight->GetXaxis()->GetNbins(); ++i) {
  //   if(MC->GetBinContent(i+1) == 0) continue;
  //   double tmp = Data->GetBinContent(i+1)*factor/MC->GetBinContent(i+1);
  //   weight->SetBinContent(i+1, tmp);
  //   cout <<  (i+1) << " | "  << tmp << " | " << 
  //     tmp*MC->GetBinContent(i+1) << " | " << ((tmp!=0) ? Data->GetBinContent(i+1)/(tmp*MC->GetBinContent(i+1)) : 0) << " | " << oldextra->GetBinContent(i+1) << endl;
  // }
  // weight->Write();
  // newfile->Close();
}
