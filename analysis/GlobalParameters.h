#pragma once
#include "iostream"
#include "TString.h"
#include "TROOT.h"

#define VERBOSITY 0
#define PRINT(s) std::cout <<#s << " " << s << std::endl

class GlobalParameters {
public:

  GlobalParameters(TString parameters_filename)
  {
    FILE * file = fopen(parameters_filename, "r");
    TString s;

    while (s.Gets(file)) {
      if (s.BeginsWith("#")) continue;
      if (s.Contains("RPC_TIME_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("RPC_TIME_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          RPC_TIME_CUT = s.Atof();
          if (VERBOSITY > 0) PRINT(RPC_TIME_CUT);
        }
      }
      if (s.Contains("RPC_XY_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("RPC_XY_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          RPC_XY_CUT = s.Atof();
          if (VERBOSITY > 0) PRINT(RPC_XY_CUT);
        }
      }
      if (s.Contains("ENERGY_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("ENERGY_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          ENERGY_CUT = s.Atof();
          if (VERBOSITY > 0) PRINT(ENERGY_CUT);
        }
      }
      if (s.Contains("BETA_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("BETA_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          BETA_CUT = s.Atof();
          if (VERBOSITY > 0) PRINT(BETA_CUT);
        }
      }
      if (s.Contains("CLUSTER_SIZE_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("CLUSTER_SIZE_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          CLUSTER_SIZE_CUT = s.Atoi();
          if (VERBOSITY > 0) PRINT(CLUSTER_SIZE_CUT);
        }
      }
      if (s.Contains("ANGLE_CUT", TString::kIgnoreCase)) {
        s.ReplaceAll("ANGLE_CUT", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          ANGLE_CUT = s.Atof();
          if (VERBOSITY > 0) PRINT(ANGLE_CUT);
        }
      }
      if (s.Contains("N_MASS_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_MASS_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_MASS_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_MASS_BINS);
        }
      }
      if (s.Contains("N_PT_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_PT_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_PT_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_PT_BINS);
        }
      }
      if (s.Contains("N_Y_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_Y_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_Y_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_Y_BINS);
        }
      }
      if (s.Contains("PATH_TO_URQMD_DATA", TString::kIgnoreCase)) {
        s.ReplaceAll("PATH_TO_URQMD_DATA", "");
        s.ReplaceAll(" ", "");
        // if(s.IsAlpha())
        sprintf(PATH_TO_URQMD_DATA, "%s", (s.Data()));
        if (VERBOSITY > 0) PRINT(PATH_TO_URQMD_DATA);
      }
      if (s.Contains("PATH_TO_CORRECTION_DATA", TString::kIgnoreCase)) {
        s.ReplaceAll("PATH_TO_CORRECTION_DATA", "");
        s.ReplaceAll(" ", "");
        // if(s.IsAlpha())
        sprintf(PATH_TO_CORRECTION_DATA, "%s", (s.Data()));
        if (VERBOSITY > 0) PRINT(PATH_TO_CORRECTION_DATA);
      }
      if (s.Contains("PATH_TO_DST_EXP", TString::kIgnoreCase)) {
        s.ReplaceAll("PATH_TO_DST_EXP", "");
        s.ReplaceAll(" ", "");
        // if(s.IsAlpha())
        sprintf(PATH_TO_DST_EXP, "%s", (s.Data()));
        if (VERBOSITY > 0) PRINT(PATH_TO_DST_EXP);
      }
      if (s.Contains("PATH_TO_DST_SIM", TString::kIgnoreCase)) {
        s.ReplaceAll("PATH_TO_DST_SIM", "");
        s.ReplaceAll(" ", "");
        // if(s.IsAlpha())
        sprintf(PATH_TO_DST_SIM, "%s", (s.Data()));
        if (VERBOSITY > 0) PRINT(PATH_TO_DST_SIM);
      }
      if (s.Contains("N_PT_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_PT_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_PT_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_PT_BINS);
        }
      }
      if (s.Contains("N_Y_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_Y_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_Y_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_Y_BINS);
        }
      }
      if (s.Contains("N_MASS_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_MASS_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsDigit()) {
          N_MASS_BINS = s.Atoi();
          if (VERBOSITY > 0) PRINT(N_MASS_BINS);
        }
      }
      if (s.Contains("SECTOR_0", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_0", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_0 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_0);
        }
      }
      if (s.Contains("SECTOR_1", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_1", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_1 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_1);
        }
      }
      if (s.Contains("SECTOR_2", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_2", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_2 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_2);
        }
      }
      if (s.Contains("SECTOR_3", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_3", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_3 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_3);
        }
      }
      if (s.Contains("SECTOR_4", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_4", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_4 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_4);
        }
      }
      if (s.Contains("SECTOR_5", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_5", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_5 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_5);
        }
      }

      if (s.Contains("SECTOR_5", TString::kIgnoreCase)) {
        s.ReplaceAll("SECTOR_5", "");
        s.ReplaceAll(" ", "");
        if (s.IsBin()) {
          SECTOR_5 = (s.Atoi() != 0);
          if (VERBOSITY > 0) PRINT(SECTOR_5);
        }
      }

      if (s.Contains("MIDRAPIDITY", TString::kIgnoreCase)) {
        s.ReplaceAll("MIDRAPIDITY", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          MIDRAPIDITY = s.Atof();
          if (VERBOSITY > 0) PRINT(MIDRAPIDITY);
        }
      }

      if (s.Contains("N_PHI_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_PHI_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          N_PHI_BINS = s.Atof();
          if (VERBOSITY > 0) PRINT(N_PHI_BINS);
        }
      }

      if (s.Contains("N_PT_FLOW_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_PT_FLOW_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          N_PT_FLOW_BINS = s.Atof();
          if (VERBOSITY > 0) PRINT(N_PT_FLOW_BINS);
        }
      }

      if (s.Contains("N_Y_FLOW_BINS", TString::kIgnoreCase)) {
        s.ReplaceAll("N_Y_FLOW_BINS", "");
        s.ReplaceAll(" ", "");
        if (s.IsFloat()) {
          N_Y_FLOW_BINS = s.Atof();
          if (VERBOSITY > 0) PRINT(N_Y_FLOW_BINS);
        }
      }
    }
    fclose(file);
  }

  void print(void) const
  {
    PRINT(PATH_TO_URQMD_DATA);
    PRINT(PATH_TO_CORRECTION_DATA);
    PRINT(PATH_TO_DST_EXP);
    PRINT(PATH_TO_DST_SIM);
    PRINT(N_MASS_BINS);
    PRINT(N_PT_BINS);
    PRINT(N_Y_BINS);
    PRINT(RPC_TIME_CUT);
    PRINT(RPC_XY_CUT);
    PRINT(ENERGY_CUT);
    PRINT(CLUSTER_SIZE_CUT);
    PRINT(BETA_CUT);
    PRINT(ANGLE_CUT);
    PRINT(N_MASS_BINS);
    PRINT(N_PT_BINS);
    PRINT(N_Y_BINS);
    PRINT(SECTOR_0);
    PRINT(SECTOR_1);
    PRINT(SECTOR_2);
    PRINT(SECTOR_3);
    PRINT(SECTOR_4);
    PRINT(SECTOR_5);
    PRINT(N_PHI_BINS);
    PRINT(N_PT_FLOW_BINS);
    PRINT(N_Y_FLOW_BINS);
  }

  Int_t N_MASS_BINS, N_PT_BINS, N_Y_BINS, N_PHI_BINS, N_Y_FLOW_BINS, N_PT_FLOW_BINS;
  Float_t RPC_TIME_CUT, RPC_XY_CUT, RK_CUT, ENERGY_CUT, BETA_CUT, ANGLE_CUT, MIDRAPIDITY;
  char INPUT[256];
  char PATH_TO_URQMD_DATA[256];
  char PATH_TO_DST_EXP[256];
  char PATH_TO_DST_SIM[256];
  char PATH_TO_CORRECTION_DATA[256];
  Bool_t SECTOR_0, SECTOR_1, SECTOR_2, SECTOR_3, SECTOR_4, SECTOR_5;
  Int_t CLUSTER_SIZE_CUT;
};
