{
  Hades* myHades=new Hades;
  HSpectrometer* spec=gHades->getSetup();
  HRuntimeDb* rtdb=gHades->getRuntimeDb();

  HEmcDetector* emc=new HEmcDetector;
  spec->addDetector(emc);
  Int_t emcMods   [1]  = {1};
  for (Int_t is=0; is<6; is++) {
    spec->getDetector("Emc")->setModules(is,emcMods);
  }
  
  HParAsciiFileIo* inp1=new HParAsciiFileIo;
  TString inputFile="./emcCellStatus_all_modules_work.txt";
  if (!inp1->open(inputFile.Data())) {
    printf("*** File %s not found ****\n",inputFile.Data());
    return;
  }
  rtdb->setFirstInput(inp1);

  HParOra2Io* ora=new HParOra2Io;
  ora->open("emc_oper");
  rtdb->setOutput(ora);
  HOra2Info* oraInfo=ora->getOra2Info();
 
  Int_t runId=1;

  HEmcCellStatusPar* pPar= rtdb->getContainer("EmcCellStatusPar");

  rtdb->initContainers(runId);
  if(pPar->hasChanged()){
    pPar->setAuthor("A. Prozorov");
    pPar->setDescription("EmcCellStatus parameters for ECAL test 5th sector 07.09.20");

    pPar->print();

    pPar->write(ora);
  }

  delete myHades;
}
