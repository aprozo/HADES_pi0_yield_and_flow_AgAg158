{  // test thermal source of eta's   (R.H. 28/7/2000)
//gROOT->Reset();
Float_t Eb    = 1.58;     // beam energy in AGeV
Float_t ctr   = 1.0;     //
Float_t bmax  = 7.6;         // Poisson sampling
if (bmax>0.) ctr = 1.;   // use b sampling instead of Poisson sampling

Float_t T1    = 0.049;   // temperature in GeV (for pion 2-component spectra)
Float_t T2    = 0.093;   // temperature in GeV (assume this for thermalized source)
Float_t frac  = 0.98;     // fraction of pion low-T component (from Jehad's fit to QMD)


Float_t blast = 0.0;     // radial expansion velocity
Float_t A2    = 0.0;     // polar distribution
Float_t A4    = 0.0;
Float_t v1    = 0.0;     // side flow
Float_t v2    = 0.0;     // elliptic flow

for (Int_t i=15;i<25;i++){
PFireball *source=new PFireball("pi0",Eb,T1,T2,frac,blast,A2,A4,v1,v2);

//source->SetW(0.1);  // set a weight
//source->Print();
PParticle *pi0=new PParticle("pi0");
PParticle *s1[]={source,pi0};
PParticle *g1=new PParticle("g");
PParticle *g2=new PParticle("g");
PParticle *s2[]={pi0,g1,g2};

PChannel  *c1=new PChannel(s1,1,1);
PChannel  *c2=new PChannel(s2);
PChannel  *cc[]={c1,c2};

//PChannel  *cc[]={c1};

PReaction *r=new PReaction(cc,Form("pi0_10mil_evts_%i",i),2,0,0,0,1);
r->Print();
r->setHGeant(0);   // set to 1, if PLUTO run from HGeant prompt

r->loop(10000000);
}
}
