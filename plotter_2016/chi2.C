int chi2() {
  TFile* Datafile = new TFile("Data.root");
  TFile* MCfile = new TFile("MC.root");

  string path = "NRecoBJet";
  const char* cpath = path.c_str();

  Datafile->cd(cpath);

  TIter nextkey( gDirectory->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {
    //keep only the highest cycle number for each key
    
    TObject *obj = key->ReadObj();
    if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      TH1D* data = (TH1D*)obj;
      MCfile->cd(cpath);
      TH1D* MC = (TH1D*)((TKey*)gDirectory->GetListOfKeys()->FindObject(data->GetName()))->ReadObj();
      MC->Scale(data->Integral()/MC->Integral());
      cout << data->GetName() << " " << MC->GetName() << endl;      
      data->Chi2Test(MC, "P");
      cout << endl;
    }
  }
  return 0;
}
