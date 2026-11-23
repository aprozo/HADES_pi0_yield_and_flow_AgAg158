#include "TString.h"
#include "TROOT.h"
#include <iostream>

#include "photons_to_histograms.C"

int main(int argc, char **argv)
{
    TROOT Analysis("photons_to_histograms","Select, combine photons to diphotons and fill the histograms");
    
    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1
    std::cout<<argc<<" Arguments: "<<std::endl;
    if(argc>1) std::cout<<"arg1 ="<<argv[1]<<std::endl; // input filename
    if(argc>2) std::cout<<"arg2 ="<<argv[2]<<std::endl; // output filename
    if(argc>3) std::cout<<"arg3 ="<<argv[3]<<std::endl; // Number of events
    if(argc>4) std::cout<<"arg4 ="<<argv[4]<<std::endl; // Is simulation?
    
    TString nEvents       ;
    TString isExperiment  ;
    switch (argc)
    {
        case 3:
            return photons_to_histograms(TString(argv[1]), TString(argv[2]), kTRUE) ? 0 : 1;
            break;
        case 4:       // just inputfile name + nEvents
            isExperiment = argv[3];
            return photons_to_histograms(TString(argv[1]), TString(argv[2]), atoi(isExperiment.Data())!=0) ? 0 : 1;
            break;
        case 5:       // inputfile name + nEvents + isExperiment
            isExperiment = argv[3];
            nEvents      = argv[4];
            return photons_to_histograms(TString(argv[1]), TString(argv[2]), atoi(isExperiment.Data())!=0) ? 0 : 1;
            break;
        default:
            cerr<<"ERROR: photons_to_histograms() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="", isExperiment=1, nevents=-1"<<std::endl;
            return 1; // fail
    }
}
