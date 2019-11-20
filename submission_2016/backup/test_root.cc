#include <TFile.h>


int main (int argc, char* argv[]) {
  TFile f;
  for(size_t i=1; i<argc; i++){
    f.Open(argv[i]);
    f.Close();
  }
}

