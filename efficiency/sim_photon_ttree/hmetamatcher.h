#ifndef HMetaMatcher_h
#define HMetaMatcher_h




#include "hades.h"
#include "hspectrometer.h"
#include "hruntimedb.h"
#include "hcategory.h"
#include "hcategorymanager.h"
#include "hlocation.h"
#include "hgeomvector.h"
#include "hgeomvolume.h"
#include "hgeomcompositevolume.h"
#include "hgeomtransform.h"
#include "hreconstructor.h"

#include "hrpcgeompar.h"
#include "hrpcgeomcellpar.h"
#include "hemcgeompar.h"
#include "htofgeompar.h"
#include "hspecgeompar.h"

#include "hparticlecand.h"
#include "hrpccluster.h"
#include "hemccluster.h"
#include "htofcluster.h"

#include "rpcdef.h"
#include "emcdef.h"
#include "tofdef.h"

#include "TMath.h"

// from hydra
class HCategory;
class HGeomVector;
class HMetaMatchPar;
class HRpcGeomPar;
class HEmcGeomPar;
class HTofGeomPar;

// from ROOT
class     TVector3;

class HMetaMatcher : public HReconstructor {

protected:

    HCategory*       fCatEmcCluster;  // pointer to the Emc cluster category
    HCategory*       fCatTofHit;      // pointer to the Tof hit category
    HCategory*       fCatTofCluster;  // pointer to the Tof cluster category
    HCategory*       fCatRpcCluster;  // pointer to the Rpc cluster category.

    HRpcGeomPar*     fRpcGeometry;    // Rpc geometry
    HRpcGeomCellPar* pGeomCellPar;
    HEmcGeomPar*     fEmcGeometry;    // Emc geometry
    HTofGeomPar*     fTofGeometry;    // TOF geometry
    HSpecGeomPar*    pSpecGeomPar;

    HGeomTransform   labSecTrans[6];       // [sec]
    HGeomTransform   modSecTransRpc[6];    // [sec]
    HGeomTransform   modSecTransEmc[6];    // [sec]
    HGeomTransform   modSecTransTof[6][8]; // [sec][mod]

    TVector3         normVecRpc[6][3];   // [sec][plane col{0,2,4}, {1,3,5}, mid]  normal vector on the RPC module in the sector coordinate system
    TVector3         normVecEmc[6];      // [sec] normal vector on the Emc module in the sector coordinate system
    TVector3         normVecTof[6][8];   // [sec][mod] normal vector on each Tof module in the sector coordinate system

    TVector3         centerVecRpc[6][3]; // [sec][plane col{0,2,4}, {1,3,5}, mid] center vector on the RPC module in the sector coordinate system
    TVector3         centerVecEmc[6];    // [sec] center vector on the Emc module in the sector coordinate system
    TVector3         centerVecTof[6][8]; // [sec][mod] center vector on each Tof module in the sector coordinate system

    HLocation fLoc;
    Float_t DPlanesRpc;

    HRpcCluster rpccluster;
    HTofCluster tofcluster;
    HEmcCluster emccluster;

    HGeomVector  trackMetaMod;
    HGeomVector  trackMetaSec;
    HGeomVector  trackMetaLab;

    Bool_t fDebug;
    void   calcSegPoints(TVector3& p1,TVector3& p2,Double_t phiseg,Double_t thetaseg,Double_t roseg,Double_t zseg);
    Bool_t findIntersectionLinePlane(TVector3 &pointIntersect,
				     const TVector3 &pos, const TVector3 &dir,
				     const TVector3 &planeCenter, const  TVector3 &planeNormal);
    Bool_t traceToMeta(const HParticleCand* cand,TVector3 &metaHit, const TVector3 &metaCenter, const TVector3 &metaNorm);

public:

    HMetaMatcher(const Text_t *name = "HMetaMatcher",const Text_t *title ="HMetaMatcher");

    virtual ~HMetaMatcher() { ; }
    void     setDebug(Bool_t debug=kTRUE) {fDebug = debug;}

    HRpcCluster*  matchWithRpc(HParticleCand* cand);
    const HTofCluster*  matchWithTof(HParticleCand* cand);
    const HEmcCluster*  matchWithEmc(HParticleCand* cand);

    const HGeomVector&  getTrackMetaMod() { return trackMetaMod;}  // only valid direct following matchWithXXX() call!
    const HGeomVector&  getTrackMetaSec() { return trackMetaSec;}  // only valid direct following matchWithXXX() call!
    const HGeomVector&  getTrackMetaLab() { return trackMetaLab;}  // only valid direct following matchWithXXX() call!

    Bool_t  init           ();
    Bool_t  reinit         ();
    Int_t   execute()  { return 0; }
    Bool_t  finalize() { return kTRUE; }
};




HMetaMatcher::HMetaMatcher(const Text_t *name,const Text_t *title)
	: HReconstructor(name,title) {

    fCatRpcCluster  = NULL;
    fCatEmcCluster  = NULL;
    fCatTofCluster  = NULL;
    fCatTofHit      = NULL;

    fRpcGeometry    = NULL;
    fEmcGeometry    = NULL;
    fTofGeometry    = NULL;
    pSpecGeomPar    = NULL;

    DPlanesRpc = 0;
    fDebug = kFALSE;
}

void HMetaMatcher::calcSegPoints(TVector3& p1,TVector3& p2,Double_t phiseg,Double_t thetaseg,Double_t roseg,Double_t zseg)
{
    // in sector sys : calc p1 p2 points (not MDC Hit!) on MdcSeg

    Double_t Xseg[2],Yseg[2],Zseg[2];

    Xseg[0]=roseg*cos(phiseg+TMath::Pi()/2);
    Yseg[0]=roseg*sin(phiseg+TMath::Pi()/2);
    Zseg[0]=zseg;
    Xseg[1]=Xseg[0]+cos(phiseg)*sin(thetaseg);
    Yseg[1]=Yseg[0]+sin(phiseg)*sin(thetaseg);
    Zseg[1]=Zseg[0]+cos(thetaseg);

    p1.SetXYZ(Xseg[0],Yseg[0],Zseg[0]);
    p2.SetXYZ(Xseg[1],Yseg[1],Zseg[1]);
}

Bool_t HMetaMatcher::findIntersectionLinePlane(TVector3 &pointIntersect,
					       const TVector3 &pos, const TVector3 &dir,
					       const TVector3 &planeCenter, const  TVector3 &planeNormal) {
    // Finds the intersection point of a straight line with any plane.
    // pointIntersect: the intersection point (return value).
    // pos:            a point on the straight line.
    // dir:            direction of the straight line.
    // planeCenter:    a point on the plane.
    // planeNormal:    normal vector of the plane.

    Double_t denom = planeNormal.Dot(dir);
    if (denom != 0.0) {
	Double_t t = ((planeCenter.x() - pos.x()) * planeNormal.x() +
		      (planeCenter.y() - pos.y()) * planeNormal.y() +
		      (planeCenter.z() - pos.z()) * planeNormal.z()) / denom;

	pointIntersect = pos + (t * dir);
	return kTRUE;
    } else {
	cout<<"Warning \"findIntersection()\" : , No intersection point found : (plane || track)"<<endl;
        return kFALSE;
    }
    return kFALSE;
}


Bool_t HMetaMatcher::traceToMeta(const HParticleCand* cand,TVector3 &metaHit, const TVector3 &metaCenter, const TVector3 &metaNorm)
{
    // calculate intersection of outer MdcSeg (RK) straightline with META plane
    // in sector sys

    TVector3 p1;
    TVector3 p2;
    TVector3 dir;

    Double_t phiLab   = cand->getPhi2();    // outer MdcSeg (RK) 0-360 deg  lab sys

    Double_t phi = fmod(phiLab,60.F) + 60;  // sec sys

    Double_t theta = cand->getTheta2();  // deg
    Double_t r     = cand->getR2();      // mm
    Double_t z     = cand->getZ2();      // mm

    calcSegPoints(p1,p2,phi*TMath::DegToRad(),theta*TMath::DegToRad(),r,z);  // 2 points from segment in sec sys

    dir = p2 - p1;

    return  findIntersectionLinePlane(metaHit,p1,dir,metaCenter,metaNorm); // in sec system
}


Bool_t HMetaMatcher::init()
{
    // Get categories and coordinate transformations.
    // Call this after HADES has been initialized.

    if(!gHades) {
        Error("init()", "HADES not initialized.");
        return kFALSE;
    }

    HRuntimeDb *rtdb=gHades->getRuntimeDb();

    if(!rtdb) {
        Error("init()", "HADES runtime database not found.");
        return kFALSE;
    }

    fCatRpcCluster  = HCategoryManager::getCategory(catRpcCluster);
    fCatEmcCluster  = HCategoryManager::getCategory(catEmcCluster);
    fCatTofCluster  = HCategoryManager::getCategory(catTofCluster);
    fCatTofHit      = HCategoryManager::getCategory(catTofHit);

    if(!fCatEmcCluster) Warning("init()", "Category catEmcCluster not found.");
    if(!fCatRpcCluster) Warning("init()", "Category catRpcCluster not found.");
    if(!fCatTofCluster) Warning("init()", "Category catTofCluster not found.");
    if(!fCatTofHit)     Warning("init()", "Category catTofHit not found.");


    fRpcGeometry = (HRpcGeomPar*)    rtdb->getContainer("RpcGeomPar");
    pGeomCellPar = (HRpcGeomCellPar*)rtdb->getContainer("RpcGeomCellPar");
    pSpecGeomPar = (HSpecGeomPar*)   rtdb->getContainer("SpecGeomPar");
    fEmcGeometry = (HEmcGeomPar *)   rtdb->getContainer("EmcGeomPar");
    fTofGeometry = (HTofGeomPar *)   rtdb->getContainer("TofGeomPar");

    return kTRUE;
}

Bool_t HMetaMatcher::reinit() {

    if(pGeomCellPar){
	DPlanesRpc  = pGeomCellPar->getDPlanes();
    }

    // Geometry transformation from lab to sector coordinate system.
    for(Int_t i = 0; i < 6; i++) {
        labSecTrans[i] = pSpecGeomPar->getSector(i)->getTransform();
    }

    for(Int_t iSec = 0; iSec < 6; iSec++) {
	if(fRpcGeometry) {
	    HModGeomPar *module = fRpcGeometry->getModule(iSec,0);
	    if(module) {
		// Transformation module->lab.
                HGeomTransform modTrans(module->getLabTransform());
                // Combine with transformation lab->sector
                // to get the transformation from module to sector.
                modTrans.transTo(labSecTrans[iSec]);

                modSecTransRpc[iSec] = modTrans;

                // Centre point plane 0
                HGeomVector r0_mod0(0.0, 0.0, -0.5*DPlanesRpc);  // plane col 0,2,4
                // Normal vector.
                HGeomVector rz_mod(0.0, 0.0, 1.0);
                HGeomVector nRpc = modTrans.transFrom(rz_mod) - modTrans.transFrom(r0_mod0);
                normVecRpc[iSec][0].SetXYZ(nRpc.getX(), nRpc.getY(), nRpc.getZ());
                HGeomVector cRpc = modTrans.transFrom(r0_mod0);
                centerVecRpc[iSec][0].SetXYZ(cRpc.getX(), cRpc.getY(), cRpc.getZ());

                // Centre point plane 1
                HGeomVector r0_mod1(0.0, 0.0, 0.5*DPlanesRpc);   // plane col 1,3,5

		nRpc = modTrans.transFrom(rz_mod) - modTrans.transFrom(r0_mod1);
                normVecRpc[iSec][1].SetXYZ(nRpc.getX(), nRpc.getY(), nRpc.getZ());
                cRpc = modTrans.transFrom(r0_mod1);
                centerVecRpc[iSec][1].SetXYZ(cRpc.getX(), cRpc.getY(), cRpc.getZ());

                HGeomVector r0_mod2(0.0, 0.0, 0);                 // plane between cols for cluster type 2

		nRpc = modTrans.transFrom(rz_mod) - modTrans.transFrom(r0_mod2);
                normVecRpc[iSec][2].SetXYZ(nRpc.getX(), nRpc.getY(), nRpc.getZ());
                cRpc = modTrans.transFrom(r0_mod2);
                centerVecRpc[iSec][2].SetXYZ(cRpc.getX(), cRpc.getY(), cRpc.getZ());

	    }
        }

	if (fEmcGeometry) {
	    HModGeomPar *pmodgeom = fEmcGeometry->getModule(iSec);
	    HGeomTransform modTrans(pmodgeom->getLabTransform());

	    modTrans.transTo(labSecTrans[iSec]);

	    modSecTransEmc[iSec] = modTrans;

	    HGeomVector r0_mod(0.0, 0.0, 0.);
	    HGeomVector rz_mod(0.0, 0.0, 1.);
	    HGeomVector nEmc = modTrans.transFrom(rz_mod) - modTrans.transFrom(r0_mod);
	    normVecEmc[iSec].SetXYZ(nEmc.getX(), nEmc.getY(), nEmc.getZ()) ;
	    HGeomVector cEmc = modTrans.transFrom(r0_mod);
	    centerVecEmc[iSec].SetXYZ(cEmc.getX(), cEmc.getY(), cEmc.getZ()) ;
	}

	if (fTofGeometry) {
	    for(Int_t iTofMod = 0; iTofMod < 8; iTofMod++) {
		HModGeomPar *module = fTofGeometry->getModule(iSec,iTofMod);
		if(module) {
		    HGeomTransform modTrans(module->getLabTransform());
		    modTrans.transTo(labSecTrans[iSec]);

		    modSecTransTof[iSec][iTofMod] = modTrans;

                    HGeomVector r0_mod(0.0, 0.0, 0.0);
                    HGeomVector rz_mod(0.0, 0.0, 1.0);
                    HGeomVector nTof = modTrans.transFrom(rz_mod) - modTrans.transFrom(r0_mod);
                    normVecTof[iSec][iTofMod].SetXYZ(nTof.getX(), nTof.getY(), nTof.getZ());
                    HGeomVector cTof = modTrans.transFrom(r0_mod);
                    centerVecTof[iSec][iTofMod].SetXYZ(cTof.getX(), cTof.getY(), cTof.getZ());
                }
            }
        }
    }
    return kTRUE;
}


HRpcCluster* HMetaMatcher::matchWithRpc(HParticleCand* cand) {


    if(!cand) return 0;
    if(!cand->isRpcClstUsed()) return 0;


    if(fDebug){ cout <<"matchRpc--------------------------------"<<endl;}
    Int_t sec   = cand->getSector();
    Int_t mod0  = cand->getMetaModule(0);
    Int_t mod1  = cand->getMetaModule(1);
    Int_t cell0 = cand->getMetaCell(0);
    Int_t cell1 = cand->getMetaCell(1);

    if(mod0==-1) return 0;

    Short_t type = 0;
    if(mod0>-1) type=1;
    if(mod1>-1) type=2;

    Float_t dx = cand->getRkMetaDx();
    Float_t dy = cand->getRkMetaDy();

    if(fDebug){
	HRpcCluster* pRpc    = 0 ;

	if(fCatRpcCluster) {
	    pRpc = (HRpcCluster*) fCatRpcCluster ->getObject(cand->getRpcInd());
	}

	if(pRpc){

	    Float_t x,y,z;
	    pRpc->getXYZLab(x,y,z);
	    cout<<"rpc sec "<<pRpc->getSector()<<" col "<<pRpc->getColumn1()<<" cell "<<setw(2)<<pRpc->getCell1()<<" type "<<type
		<<" ,pos mod "<<setw(10)<<pRpc->getXMod()
		<<" "<<setw(10)<<pRpc->getYMod()
		<<" "<<setw(10)<<pRpc->getZMod()
		<<" ,pos sec "<<setw(10)<<pRpc->getXSec()
		<<" "<<setw(10)<<pRpc->getYSec()
		<<" "<<setw(10)<<pRpc->getZSec()
		<<" ,pos lab "<<setw(10)<<x
		<<" "<<setw(10)<<y
		<<" "<<setw(10)<<z
		<<" ,theta "<<setw(10)<<pRpc->getTheta()
		<<" phi "<<setw(10)<<pRpc->getPhi()
		<<" dx "<<dx
		<<" dy "<<dy
                <<" type "<<pRpc->getClusterType()
		<<endl;
	}
    }

    TVector3 metaHit;
    Int_t ind = mod0%2; // finding the correct plane
    if(mod1>-1) ind=2;  // both cells, plane in middle

    traceToMeta(cand,metaHit,centerVecRpc[sec][ind],normVecRpc[sec][ind]);
    HGeomVector vsec(metaHit.X(),metaHit.Y(),metaHit.Z());    // hit point mdc on META sec sys

    trackMetaSec = vsec;

    HGeomTransform& TransSec = pSpecGeomPar->getSector(sec)->getTransform();
    HGeomVector vlab  = TransSec.transFrom(vsec); // lab sys

    trackMetaLab = vlab;

    HGeomTransform& modToLab=fRpcGeometry->getModule(sec,0)->getLabTransform();
    HGeomVector vmod = modToLab.transTo(vlab);  // mod sys

    trackMetaMod = vmod;

    vmod.setX(vmod.X()+dx);   // reconstruct the meta hit
    vmod.setY(vmod.Y()+dy);   // reconstruct the meta hit

    vlab = modToLab.transFrom(vmod); // lab sys
    vsec = TransSec.transTo(vlab);   // sec sys

    Double_t theta = TMath::RadToDeg() * TMath::ATan2(TMath::Sqrt(vlab.X()*vlab.X()+vlab.Y()*vlab.Y()),vlab.Z());
    Double_t phi   = TMath::RadToDeg() * TMath::ATan2(vlab.Y(),vlab.X());
    if (phi < 0.) phi += 360.;

    if(fDebug){
	cout<<"mdc sec "<<sec<<" col "<<mod0<<" cell "<<setw(2)<<cell0<<" type "<<type
	    <<" ,pos mod "<<setw(10)<<vmod.X()
	    <<" "<<setw(10)<<vmod.Y()
	    <<" "<<setw(10)<<vmod.Z()
	    <<" ,pos sec "<<setw(10)<<vsec.X()
	    <<" "<<setw(10)<<vsec.Y()
	    <<" "<<setw(10)<<vsec.Z()
	    <<" ,pos lab "<<setw(10)<<vlab.X()
	    <<" "<<setw(10)<<vlab.Y()
	    <<" "<<setw(10)<<vlab.Z()
	    <<" ,theta "<<setw(10)<<theta
	    <<" phi "<<setw(10)<<phi
            <<" type "<<type
	    <<endl;
    }
    rpccluster.setDetID1(sec,mod0,cell0);
    rpccluster.setDetID2(sec,mod1,cell1);
    rpccluster.setTheta(theta);
    rpccluster.setPhi(phi);
    rpccluster.setXYZLab(vlab.X(),vlab.Y(),vlab.Z());
    rpccluster.setRSEC(vsec.X(),vsec.Y(),vsec.Z());
    rpccluster.setCluster(cand->getTof(),1.,vmod.X(),vmod.Y(),vmod.Z());
    rpccluster.setClusterType(type);

    return &rpccluster;

}

const HTofCluster* HMetaMatcher::matchWithTof(HParticleCand* cand) {


    if(!cand) return 0;
    if(! (cand->isTofClstUsed() || cand->isTofHitUsed()) ) return 0;


    if(fDebug){ cout <<"matchTof--------------------------------"<<endl;}
    Int_t sec   = cand->getSector();
    Int_t mod0  = cand->getMetaModule(0);
    Int_t mod1  = cand->getMetaModule(1);
    Int_t cell0 = cand->getMetaCell(0);

    if(mod0==-1) return 0;

    Short_t type = 0;
    if(mod0>-1) type=1;
    if(mod1>-1) type=2;

    Float_t dx = cand->getRkMetaDx();
    Float_t dy = cand->getRkMetaDy();

    if(fDebug){
	HTofCluster* pTof    = 0 ;

	if(cand->isTofClstUsed() && fCatTofCluster) {
	    pTof = (HTofCluster*) fCatTofCluster ->getObject(cand->getTofClstInd());
	}
	if(cand->isTofHitUsed() && fCatTofHit) {
	    pTof = (HTofCluster*) fCatTofHit ->getObject(cand->getTofHitInd());
	}

	if(pTof){

	    Float_t x,y,z;
	    Float_t ph,th;
            pTof->getTheta(th);
            pTof->getPhi(ph);
	    pTof->getXYZLab(x,y,z);
	    cout<<"tof sec "<<(Int_t)pTof->getSector()<<" mod "<<(Int_t)pTof->getModule()<<" cell "<<setw(2)<<(Int_t)pTof->getCell()<<" type "<<type
		<<" ,pos mod "<<setw(10)<<pTof->getXpos()
		<<" ,pos lab "<<setw(10)<<x
		<<" "<<setw(10)<<y
		<<" "<<setw(10)<<z
		<<" ,theta "<<setw(10)<<th
		<<" phi "<<setw(10)<<ph
		<<" dx "<<dx
		<<" dy "<<dy
		<<endl;
	}
    }

    TVector3 metaHit;

    traceToMeta(cand,metaHit,centerVecTof[sec][mod0],normVecTof[sec][mod0]);
    HGeomVector vsec(metaHit.X(),metaHit.Y(),metaHit.Z());    // hit point mdc on META sec sys

    trackMetaSec = vsec;

    HGeomTransform& TransSec = pSpecGeomPar->getSector(sec)->getTransform();
    HGeomVector vlab  = TransSec.transFrom(vsec); // lab sys

    trackMetaLab = vlab;

    HGeomTransform& modToLab=fTofGeometry->getModule(sec,mod0)->getLabTransform();
    HGeomVector vmod = modToLab.transTo(vlab);  // mod sys

    trackMetaMod = vmod;

    vmod.setX(vmod.X()+dx);   // reconstruct the meta hit
    vmod.setY(vmod.Y()+dy);   // reconstruct the meta hit

    vlab = modToLab.transFrom(vmod); // lab sys
    vsec = TransSec.transTo(vlab);   // sec sys

    Double_t theta = TMath::RadToDeg() * TMath::ATan2(TMath::Sqrt(vlab.X()*vlab.X()+vlab.Y()*vlab.Y()),vlab.Z());
    Double_t phi   = TMath::RadToDeg() * TMath::ATan2(vlab.Y(),vlab.X());
    if (phi < 0.) phi += 360.;

    if(fDebug){
	cout<<"mdc sec "<<sec<<" mod "<<mod0<<" cell "<<setw(2)<<cell0<<" type "<<type
	    <<" ,pos mod "<<setw(10)<<vmod.X()
	    <<" ,pos lab "<<setw(10)<<vlab.X()
	    <<" "<<setw(10)<<vlab.Y()
	    <<" "<<setw(10)<<vlab.Z()
	    <<" ,theta "<<setw(10)<<theta
	    <<" phi "<<setw(10)<<phi
	    <<endl;
    }
    tofcluster.setSector(sec);
    tofcluster.setModule(mod0);
    tofcluster.setCell(cell0);
    tofcluster.setTheta(theta);
    tofcluster.setPhi(phi);
    tofcluster.setXYZLab(vlab.X(),vlab.Y(),vlab.Z());
    tofcluster.setXpos(vmod.X());
    tofcluster.setTof(cand->getTof());
    tofcluster.setEdep(cand->getTofdEdx());
    tofcluster.setClusterSize(type);

    return &tofcluster;

}

const HEmcCluster* HMetaMatcher::matchWithEmc(HParticleCand* cand)
{
    // results not exact at the moment! No dx,dy available for Emc.


    if(!cand) return 0;
    if(cand->getEmcInd()==-1) return 0;


    if(fDebug){ cout <<"matchEmc--------------------------------"<<endl;}
    Int_t sec   = cand->getSector();

    Short_t type = 0;

    Float_t dx = cand->getMetaMatchRadiusEmc();    // currently no dx,dy available!!!!!
    Float_t dy = cand->getMetaMatchRadiusEmc();

    if(fDebug){
	HEmcCluster* pEmc    = 0 ;

	if(cand->getEmcInd()>-1 && fCatEmcCluster) {
	    pEmc = (HEmcCluster*) fCatEmcCluster ->getObject(cand->getEmcInd());
	}

	if(pEmc){

	    Float_t x,y,z;
	    pEmc->getXYZLab(x,y,z);
	    cout<<"emc sec "<<(Int_t)pEmc->getSector()<<" cell "<<setw(3)<<(Int_t)pEmc->getCell()<<" type "<<(Int_t)pEmc->getNCells()
		<<" ,pos mod "<<setw(10)<<pEmc->getXMod()
		<<" "<<setw(10)<<pEmc->getYMod()
		<<" ,pos lab "<<setw(10)<<x
		<<" "<<setw(10)<<y
		<<" "<<setw(10)<<z
		<<" ,theta "<<setw(10)<<pEmc->getTheta()
		<<" phi "<<setw(10)<<pEmc->getPhi()
		<<" dx "<<dx
		<<" dy "<<dy
		<<endl;
	}
    }

    TVector3 metaHit;

    traceToMeta(cand,metaHit,centerVecEmc[sec],normVecEmc[sec]);
    HGeomVector vsec(metaHit.X(),metaHit.Y(),metaHit.Z());    // hit point mdc on META sec sys

    trackMetaSec = vsec;

    HGeomTransform& TransSec = pSpecGeomPar->getSector(sec)->getTransform();
    HGeomVector vlab  = TransSec.transFrom(vsec); // lab sys

    trackMetaLab = vlab;

    HGeomTransform& modToLab=fEmcGeometry->getModule(sec,0)->getLabTransform();
    HGeomVector vmod = modToLab.transTo(vlab);  // mod sys

    trackMetaMod = vmod;

    vmod.setX(vmod.X()+dx);   // reconstruct the meta hit
    vmod.setY(vmod.Y()+dy);   // reconstruct the meta hit

    vlab = modToLab.transFrom(vmod); // lab sys
    vsec = TransSec.transTo(vlab);   // sec sys

    Double_t theta = TMath::RadToDeg() * TMath::ATan2(TMath::Sqrt(vlab.X()*vlab.X()+vlab.Y()*vlab.Y()),vlab.Z());
    Double_t phi   = TMath::RadToDeg() * TMath::ATan2(vlab.Y(),vlab.X());
    if (phi < 0.) phi += 360.;

    if(fDebug){
	cout<<"mdc sec "<<sec<<" cell "<<setw(3)<<0<<" type "<<type
	    <<" ,pos mod "<<setw(10)<<vmod.X()
	    <<" "<<setw(10)<<vmod.Y()
	    <<" ,pos lab "<<setw(10)<<vlab.X()
	    <<" "<<setw(10)<<vlab.Y()
	    <<" "<<setw(10)<<vlab.Z()
	    <<" ,theta "<<setw(10)<<theta
	    <<" phi "<<setw(10)<<phi
	    <<endl;
    }
    emccluster.setSector(sec);
    emccluster.setTheta(theta);
    emccluster.setPhi(phi);
    emccluster.setXYZLab(vlab.X(),vlab.Y(),vlab.Z());
    emccluster.setXYMod(vmod.X(),vmod.Y());
    emccluster.setTime(cand->getEmcTime());
    emccluster.setEnergy(cand->getEmcEnergy());

    return &emccluster;

}
#endif // HMetaMatcher_h
