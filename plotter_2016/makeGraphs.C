void CreateGraph(TDirectory*, TList*);

void makeGraphs() {
  cout << "begin" << endl;
  TFile* newfile = new TFile("compareFake.root", "RECREATE");
  TList* files = new TList();
  TFile* f4 = new TFile("W+Jets.root");
  TFile* f1 = new TFile("Sub.root");
  TFile* f2 = new TFile("Lead.root");
  TFile* f3 = new TFile("Rand.root");
  files->Add(f1);
  files->Add(f2);
  files->Add(f3);
  files->Add(f4);
  cout << "before" << endl;
  gROOT->SetBatch(kTRUE);
  CreateGraph(newfile, files);
  gROOT->SetBatch(kFALSE);
}



 void CreateGraph( TDirectory *target, TList* files) {
   cout << "start" << endl;


  TString path( (char*)strstr( target->GetPath(), ":" ) );
  path.Remove( 0, 2 );

  TFile* openFile = (TFile*)files->First();
  openFile->cd( path );
  TDirectory *current_sourcedir = gDirectory;

  //  Bool_t status = TH1::AddDirectoryStatus();
  //  TH1::AddDirectory(kFALSE);


  // loop over all keys in this directory
  TIter nextkey( current_sourcedir->GetListOfKeys() );
  TKey *key, *oldkey=0;
  while ( (key = (TKey*)nextkey())) {

    //keep only the highest cycle number for each key
    if (oldkey && !strcmp(oldkey->GetName(),key->GetName())) continue;

    openFile->cd( path );

    TObject *obj = key->ReadObj();
    if ( obj->IsA() ==  TH1D::Class() ) {

      TH1D* h1 = (TH1D*)obj;
          
      TList* allFiles = new TList();
      TFile *nextfile = (TFile*)files->After(openFile);
      string title = openFile->GetTitle();
      title = title.substr(0, title.size()-5);
      h1->SetTitle(title.c_str());
      h1->SetLineWidth(3);
      h1->SetLineColor(50);
	  
      allFiles->Add(h1);
      int nfile = 1;
      while ( nextfile ) {
	nextfile->cd( path );
	TKey *key2 = (TKey*)gDirectory->GetListOfKeys()->FindObject(h1->GetName());
	if(key2) {
	  
	  TH1D* h2 = (TH1D*)key2->ReadObj();
	  title = nextfile->GetTitle();
	  title = title.substr(0, title.size()-5);
	  h2->SetTitle(title.c_str());
	  h2->SetLineWidth(3);
	  h2->SetLineColor(50+nfile*10);
	  h2->Scale(h1->Integral()/h2->Integral());
	  allFiles->Add(h2);

	}
	nfile++;
	nextfile = (TFile*)files->After(nextfile);
      }
    
      TLegend* leg = new TLegend(0.73,0.70,0.93,0.90);

      TH1D* tmp = (TH1D*)allFiles->First();
      while( tmp ) {
	leg->AddEntry(tmp, tmp->GetTitle(), "lep");
	tmp = (TH1D*)allFiles->After(tmp);
      }

      target->cd();

      TCanvas *c = new TCanvas(h1->GetName(), h1->GetName());//403,50,600,600);
      double max = 0;
      tmp = (TH1D*)allFiles->First();
      TH1D* first = tmp;
      while(tmp) {
	max = (max < tmp->GetMaximum()) ? tmp->GetMaximum() : max;
	tmp->Draw("sameC");
      	tmp = (TH1D*)allFiles->After(tmp);
      }
      leg->Draw();
      first->GetYaxis()->SetRangeUser(0, max*(1+1./15));


      ///second pad
      c->Write(c->GetName());
      c->Close();
      allFiles->Delete();

    } else if ( obj->IsA()->InheritsFrom( TDirectory::Class() ) ) {
      target->cd();
      TDirectory *newdir = target->mkdir( obj->GetName(), obj->GetTitle() );

      CreateGraph( newdir, files );

    } else if ( obj->IsA()->InheritsFrom( TH1::Class() ) ) {
      continue;

    } else {
         cout << "Unknown object type, name: "
	   << obj->GetName() << " title: " << obj->GetTitle() << endl;
    }
  }

  //  TH1::AddDirectory(kTRUE);
}


