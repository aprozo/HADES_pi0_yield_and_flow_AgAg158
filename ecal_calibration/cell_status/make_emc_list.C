#include <iostream>
#include "TString.h"
#include "exp_setup_cells.h"
using namespace std;

void make_emc_list(){
   FILE *fp2;
  TString outpara="emcCellStatus_all_modules_work.txt";
  if ((fp2 = fopen(outpara, "w")) == NULL)
  { printf("I cannot open new params-file. \n");}       
  fprintf(fp2,"##############################################################################\n");
  fprintf(fp2,"# Class:   HEmcCellStatusPar\n");
  fprintf(fp2,"# Context: EmcCellStatusPar\n");
  fprintf(fp2,"##############################################################################\n");
  fprintf(fp2,"[EmcCellStatusPar]\n");
  fprintf(fp2,"//----------------------------------------------------------------------------\n");
  fprintf(fp2,"cellStatus: Int_t \\\n");

  for(Int_t sector=0; sector<6;sector++){
    for(Int_t cell=0; cell<163;cell++){
      Int_t counter = sector*163+cell+1;      
      fprintf(fp2,"  %u", 1);
      if(counter==5&&cell==163)  {fprintf(fp2,"\n");break;}
      if(counter%6==0) fprintf(fp2,"  \\\n");
    }
  }
  fprintf(fp2,"##############################################################################\n");      
}
