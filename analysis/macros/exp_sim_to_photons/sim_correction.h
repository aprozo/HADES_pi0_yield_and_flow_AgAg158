
static TH1D *correction=  NULL;

void setupSimEffCorr(TString filename) {
    TFile *fileSimEffCorr= new TFile(filename);
    if (fileSimEffCorr) {
        correction = (TH1D*) fileSimEffCorr->Get("correction");;
        if (!correction) {
            cout << "No histogram with sim eff correction found in file. Exitting..." << endl;
            exit(-1);
        }
    }
    else {
        cout << "No file with with sim eff correction found. Exitting..." << endl;
        exit(-1);
    }
}

Float_t correctEMCenergy(const Float_t & energy, const Int_t & sector, const Int_t & cell) {
     Double_t corr = correction->GetBinContent(sector*200+cell+1);
     Float_t emcEnCorr = energy/corr;
     return emcEnCorr;
}
