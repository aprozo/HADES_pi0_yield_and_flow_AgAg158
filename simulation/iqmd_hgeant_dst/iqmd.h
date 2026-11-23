//////////////////////////////////////////////////////////
// This class has been automatically generated on
// Wed Apr  5 10:22:11 2017 by ROOT version 5.34/34
// from TTree h90/KINE
// found on file: auau1250smw12kp-cp.5idst-5.root
//////////////////////////////////////////////////////////

#ifndef iqmd_h
#define iqmd_h

#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>

// Header file for the classes stored in the TTree if any.

// Fixed size dimensions of array or collections stored in the TTree if any.

class iqmd {
public :
   TTree          *fChain;   //!pointer to the analyzed TTree or TChain
   Int_t           fCurrent; //!current Tree number in a TChain

   // Declaration of leaf types
   Float_t         b;
   Float_t         phr;
   Float_t         time;
   Float_t         apar;
   Float_t         appa;
   Int_t           npart;
   Float_t         px[3000];   //[npart]
   Float_t         py[3000];   //[npart]
   Float_t         pz[3000];   //[npart]
   Float_t         etot[3000];   //[npart]
   Float_t         y[3000];   //[npart]
   Float_t         imod[3000];   //[npart]
   Float_t         mass[3000];   //[npart]
   Float_t         chrg[3000];   //[npart]
   Float_t         phi[3000];   //[npart]
   Float_t         thet[3000];   //[npart]
   Float_t         ncol[3000];   //[npart]
   Float_t         ri[3000];   //[npart]

   // List of branches
   TBranch        *b_b;   //!
   TBranch        *b_phr;   //!
   TBranch        *b_time;   //!
   TBranch        *b_apar;   //!
   TBranch        *b_appa;   //!
   TBranch        *b_phi;   //!
   TBranch        *b_npart;   //!
   TBranch        *b_px;   //!
   TBranch        *b_py;   //!
   TBranch        *b_pz;   //!
   TBranch        *b_etot;   //!
   TBranch        *b_y;   //!
   TBranch        *b_imod;   //!
   TBranch        *b_mass;   //!
   TBranch        *b_chrg;   //!
   TBranch        *b_thet;   //!
   TBranch        *b_ncol;   //!
   TBranch        *b_ri;   //!

   iqmd(TTree *tree=0);
   virtual ~iqmd();
   virtual Int_t    Cut(Long64_t entry);
   virtual Int_t    GetEntry(Long64_t entry);
   virtual Long64_t LoadTree(Long64_t entry);
   virtual void     Init(TTree *tree);
   virtual void     Loop(Int_t nevt);
   virtual Bool_t   Notify();
   virtual void     Show(Long64_t entry = -1);
};

#endif

#ifdef iqmd_cxx
iqmd::iqmd(TTree *tree) : fChain(0) 
{
// if parameter tree is not specified (or zero), connect the file
// used to generate this class and read the Tree.
    if (tree == 0) {

	TChain * chain = new TChain("h90","");

        chain->Add("/lustre/fopi/leifels/hadesAuAu10M/*.root/h90");
        //chain->Add("/lustre/hades/user/sudol/HADES_iQMD_1M/root/auau1250smw12kp-cp.5idst-*.root/h90");
	//chain->Add("/lustre/fopi/leifels/hadesAuAu10M/auau1230smmbt140pi_V9.5idst-*.root/h90");
	tree = chain;
      }

   Init(tree);
}

iqmd::~iqmd()
{
   if (!fChain) return;
   delete fChain->GetCurrentFile();
}

Int_t iqmd::GetEntry(Long64_t entry)
{
// Read contents of entry.
   if (!fChain) return 0;
   return fChain->GetEntry(entry);
}
Long64_t iqmd::LoadTree(Long64_t entry)
{
// Set the environment to read one entry
   if (!fChain) return -5;
   Long64_t centry = fChain->LoadTree(entry);
   if (centry < 0) return centry;
   if (fChain->GetTreeNumber() != fCurrent) {
      fCurrent = fChain->GetTreeNumber();
      Notify();
   }
   return centry;
}

void iqmd::Init(TTree *tree)
{
   // The Init() function is called when the selector needs to initialize
   // a new tree or chain. Typically here the branch addresses and branch
   // pointers of the tree will be set.
   // It is normally not necessary to make changes to the generated
   // code, but the routine can be extended by the user if needed.
   // Init() will be called many times when running on PROOF
   // (once per file to be processed).

   // Set branch addresses and branch pointers
   if (!tree) return;
   fChain = tree;
   fCurrent = -1;
   fChain->SetMakeClass(1);

   fChain->SetBranchAddress("b", &b, &b_b);
   fChain->SetBranchAddress("phr", &phr, &b_phr);
   fChain->SetBranchAddress("time", &time, &b_time);
   fChain->SetBranchAddress("apar", &apar, &b_apar);
   fChain->SetBranchAddress("appa", &appa, &b_appa);
   fChain->SetBranchAddress("npart", &npart, &b_npart);
   fChain->SetBranchAddress("px", px, &b_px);
   fChain->SetBranchAddress("py", py, &b_py);
   fChain->SetBranchAddress("pz", pz, &b_pz);
   fChain->SetBranchAddress("etot", etot, &b_etot);
   fChain->SetBranchAddress("y", y, &b_y);
   fChain->SetBranchAddress("imod", imod, &b_imod);
   fChain->SetBranchAddress("mass", mass, &b_mass);
   fChain->SetBranchAddress("chrg", chrg, &b_chrg);
   fChain->SetBranchAddress("phi", phi, &b_phi);
   fChain->SetBranchAddress("thet", thet, &b_thet);
   fChain->SetBranchAddress("ncol", ncol, &b_ncol);
   fChain->SetBranchAddress("ri", ri, &b_ri);
   Notify();
}

Bool_t iqmd::Notify()
{
   // The Notify() function is called when a new file is opened. This
   // can be either for a new TTree in a TChain or when when a new TTree
   // is started when using PROOF. It is normally not necessary to make changes
   // to the generated code, but the routine can be extended by the
   // user if needed. The return value is currently not used.

   return kTRUE;
}

void iqmd::Show(Long64_t entry)
{
// Print contents of entry.
// If entry is not specified, print current entry
   if (!fChain) return;
   fChain->Show(entry);
}
Int_t iqmd::Cut(Long64_t entry)
{
// This function may be called from Loop.
// returns  1 if entry is accepted.
// returns -1 otherwise.
   return 1;
}
#endif // #ifdef iqmd_cxx
