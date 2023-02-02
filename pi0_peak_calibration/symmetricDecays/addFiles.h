#ifndef __ADD_FILES__
#define __ADD_FILES__

#include "TTree.h"
#include "TChain.h"
#include "TH1.h"
#include "TFile.h"
#include "TString.h"
#include "TSystem.h"
#include "TFile.h"
#include "TObjString.h"
#include <iostream>
#include <stdio.h>
#include <iomanip>
using namespace std;
void printProgress(Long64_t ct,Long64_t total,Long64_t step, TString comment )
{
        // print the progress of a loop in %, where ct is the loop
        // counter, total the total number of loops and step the
        // step size in precent for which the print should be done.

        if(total <= 0 ) {
                cout<<"Error::printProgress(), total is <=0 !"<<endl;
                return;
        }
        if(ct < 0 ) {
                cout<<"Error::printProgress(), ct is < 0 !"<<endl;
                return;
        }
        if(ct == 0 ) cout<< comment.Data() << flush;
        Long64_t frac   =  (Long64_t) ((( ct     /(Float_t)total) ) * 100);
        Long64_t frac_1 =  (Long64_t) ((((ct - 1)/(Float_t)total) ) * 100);

        if( frac % step == 0 &&  (frac != frac_1) )
        {
                if(frac_1/step > 0 ) cout<<"\b\b\b\b\b"<<flush;
                cout<<setw(3)<<frac<<" %"<<flush;
                if(frac + step >= 100) cout<<endl;
        }
}

Bool_t addFile (TString infile, TChain *fChain)
{
        // add a single root file to the chain.
        // It will be checked if the file exists.s

        TObjArray* elements = fChain->GetListOfFiles();
        Int_t nfiles        = elements->GetEntries();
        TFile * temp;
        if(gSystem->AccessPathName(infile.Data()) == 0) {
                std::cout<<  std::setw(4)<<nfiles<<" add file : "<<infile.Data()<<  std::endl;
                fChain->AddFile(infile.Data());
                temp=TFile::Open(infile.Data(),"read");
                temp->Close();
                return kTRUE;
        } else {
                //Error("addFile()","File = %s not found! Will be skipped .... ",infile.Data());
                return kFALSE;
        }

}

void addMultFiles(TString commaSeparatedList, TChain *fChain)
{
        // add multiple files as comma separated file list : "file1,file2,file3"

        TObjArray* arr = commaSeparatedList.Tokenize(",");
        Int_t n = arr->GetEntries();
        if(n == 0) {
                //Warning("addMultFiles()","No file in input string '%s' found !",commaSeparatedList.Data());
        } else {
                for(Int_t i = 0; i < n; i++) {
                        TString name = ((TObjString*)arr->At(i))->GetString();
                        if(!addFile(name, fChain)) {
                                delete arr;
                                //return kFALSE;
                        }
                }
        }
        delete arr;
        // return kTRUE;
}
#endif
