#include "TString.h"
#include "TROOT.h"
#include <iostream>

#include "exp_sim_to_photons.C"

int main(int argc, char **argv)
{
    TROOT Analysis("exp_sim_to_photons","Compiled analysis macro filling the photon tree");
    
    // argc is the number of arguments in char* array argv
    // CAUTION: argv[0] contains the progname
    // argc has to be nargs+1
    std::cout<<argc<<" Arguments: "<<std::endl;
    if(argc>1) std::cout<<"input filename   = "<<argv[1]<<std::endl; // input filename
    if(argc>2) std::cout<<"output filename  = "<<argv[2]<<std::endl; // output filename
    if(argc>3) std::cout<<"sim/exp          = "<<(atoi(argv[3])==0?"SIMULATION":"EXPERIMENT")<<std::endl; // Is experiment?
    if(argc>4) std::cout<<"number of events = "<<argv[3]<<std::endl; // Number of events
    
    
    TString nEvents       ;
    TString isExperiment  ;
    switch (argc)
    {
        case 3:
            return exp_sim_to_photons(TString(argv[1]), TString(argv[2]), kTRUE, -1);
            break;
        case 4:       // just inputfile name + nEvents
            isExperiment = argv[3];
            return exp_sim_to_photons(TString(argv[1]), TString(argv[2]), atoi(isExperiment.Data())!=0, -1);
            break;
        case 5:       // inputfile name + nEvents + isExperiment
            isExperiment = argv[3];
            nEvents      = argv[4];
            return exp_sim_to_photons(TString(argv[1]), TString(argv[2]), atoi(isExperiment.Data())!=0, atoi(nEvents.Data()));
            break;
        default:
            cerr<<"ERROR: exp_sim_to_photons() : WRONG NUMBER OF ARGUMENTS! TString infile="",TString outfile="", isExperiment=1, nevents=-1"<<std::endl;
            return 1; // fail
    }
}

