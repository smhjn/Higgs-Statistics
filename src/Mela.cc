/** \file
 *
 *  See header file for documentation.
 *
 *  $Date: 2012/09/25 01:39:24 $
 *  $Revision: 1.13 $
 */

#include <Mela.h>
//#include <DataFormats/GeometryVector/interface/Pi.h>
//#include <FWCore/ParameterSet/interface/FileInPath.h>

#include "computeAngles.h"
#include "AngularPdfFactory.h"
#include "RooqqZZ_JHU_ZgammaZZ.h"
#include "RooTsallis.h"
#include "RooTsallisExp.h"
#include "RooRapidityBkg.h"
#include "RooRapiditySig.h"

#include <RooMsgService.h>
#include <TFile.h>
#include <TH1F.h>
#include <TH2F.h>
#include <TH3F.h>
#include <vector>

#include <string>

using namespace RooFit;

Mela::Mela(bool usePowhegTemplate){ 

  usePowhegTemplate_=usePowhegTemplate;

  edm::FileInPath fip("ZZMatrixElement/MELA/data/my8DTemplateNotNorm.root");
  string fullPath = fip.fullPath();

  // Original code from KDProducer::beginJob
  f = new TFile(fullPath.c_str(),"READ");
  h_mzz= (TH1F*)(f->Get("h_mzz"));
  h_mzzm1m2= (TH3F*)(f->Get("h_mzzm1m2"));
  h_mzzcosthetastar= (TH2F*)(f->Get("h_mzzcosthetastar"));
  h_mzzcostheta1= (TH2F*)(f->Get("h_mzzcostheta1"));
  h_mzzcostheta2= (TH2F*)(f->Get("h_mzzcostheta2"));
  h_mzzphi1= (TH2F*)(f->Get("h_mzzphi1")); // This is phistar1
  h_mzzphi= (TH2F*)(f->Get("h_mzzphi"));

  RooMsgService::instance().getStream(1).removeTopic(NumIntegration);

}

Mela::~Mela(){ 
  delete f;
}

double Mela::sigPdfNorm(double mzz){

  double p0 =     -16.3523;
  double p1 =     0.545736;
  double p2 =  -0.00228891;
  double p3 =  5.06627e-06;
  double p4 = -6.12505e-09;
  double p5 =  3.83126e-12;
  double p6 =  -9.7065e-16;

  double poly = p0+
    p1*mzz+
    p2*mzz*mzz+
    p3*mzz*mzz*mzz+
    p4*mzz*mzz*mzz*mzz+
    p5*mzz*mzz*mzz*mzz*mzz+
    p6*mzz*mzz*mzz*mzz*mzz*mzz;

  return exp(poly);

}

double Mela::bkgPdfNorm(double mzz){

  mzz++;
    
  double p0;
  double p1;
  double p2;
  double p3;
  double p4;
  double p5;
  double p6;
  if(mzz<182){
 
    p0 =  1.36714e-09;
    p1 = -6.38091e-11;
    p2 =  1.22195e-12;
    p3 = -1.23122e-14;
    p4 =  6.89716e-17;
    p5 = -2.03987e-19;
    p6 =  2.49154e-22;

  }else if(mzz<210){

    p0 =  -4.83519e-08;
    p1 =   9.61952e-10;
    p2 =   -7.1795e-12;
    p3 =   2.38262e-14;
    p4 =  -2.96631e-17;
    p5 = 0.0; 
    p6 = 0.0;

  }else{

    p0 =  1.15125e-11  ;
    p1 =  1.21307e-14  ;
    p2 =  -3.1029e-16  ;
    p3 =  1.01661e-18  ;
    p4 = -1.48631e-21  ;
    p5 =   1.0428e-24  ;
    p6 = -2.85659e-28  ; 

  }
  return p0+
    p1*mzz+
    p2*mzz*mzz+
    p3*mzz*mzz*mzz+
    p4*mzz*mzz*mzz*mzz+
    p5*mzz*mzz*mzz*mzz*mzz+
    p6*mzz*mzz*mzz*mzz*mzz*mzz;

}

void Mela::computeKD(TLorentzVector Z1_lept1, int Z1_lept1Id,
		     TLorentzVector Z1_lept2, int Z1_lept2Id,
		     TLorentzVector Z2_lept1, int Z2_lept1Id,
		     TLorentzVector Z2_lept2, int Z2_lept2Id,
		     // return variables:
		     float& costhetastar,
		     float& costheta1, 
		     float& costheta2,
		     float& phi,
		     float& phistar1,
		     float& kd, 
		     float& psig,
		     float& pbkg,
		     bool withPt,
		     bool withY,
		     int LHCsqrts) {
  
  //compute angles  
  float m1=(Z1_lept1 + Z1_lept2).M();
  float m2=(Z2_lept1 + Z2_lept2).M();

  TLorentzVector ZZ = (Z1_lept1 + Z1_lept2 + Z2_lept1 + Z2_lept2);
  float mzz = ZZ.M();
  float pt  = ZZ.Pt();
  float Y   = ZZ.Rapidity(); // Fixme: should probably protect against NaN?

  mela::computeAngles(Z1_lept1, Z1_lept1Id, Z1_lept2, Z1_lept2Id, 
		      Z2_lept1, Z2_lept1Id, Z2_lept2, Z2_lept2Id,
		      costhetastar,costheta1,costheta2,phi,phistar1);

  //compute kd
  pair<float,float> P = likelihoodDiscriminant(mzz,m1,m2,costhetastar,costheta1,costheta2,phi,phistar1,
					       LHCsqrts, withPt, pt, withY, Y);
  psig=P.first;
  pbkg=P.second;
  kd = psig/(psig+pbkg);
  
}



void Mela::computeKD(float mzz, float mZ1, float mZ2, 
		     float costhetastar,
		     float costheta1, 
		     float costheta2,
		     float phi,
		     float phistar1,
		     float& kd, 
		     float& psig,
		     float& pbkg,
		     bool withPt,
		     float pt4l,
		     bool withY,
		     float Y4l,
		     int LHCsqrts) {

  // Skip candidates where KD is irrelevant.
  if (mzz<100.){
    kd = 0;
    psig = 0;
    pbkg = 0;
    return;
  }

  pair<float,float> P = likelihoodDiscriminant(mzz, mZ1, mZ2, 
					       costhetastar, 
					       costheta1, 
					       costheta2, 
					       phi, 
					       phistar1,
					       LHCsqrts,
					       withPt,pt4l,
					       withY, Y4l);
  psig = P.first;
  pbkg = P.second;
  kd = psig/(psig + pbkg);

}



vector<float> Mela::my8DTemplate(bool normalized,float mZZ, float m1, float m2, float costhetastar, float costheta1, float costheta2, float phi, float phistar1){

  //multiply the P values
  float n = h_mzz->GetBinContent(h_mzz->FindBin(mZZ));
  float Pmzzm1m2 = h_mzzm1m2->GetBinContent(h_mzzm1m2->FindBin(mZZ,m1,m2));

  // - - - - - - - - - - - - - - - whitbeck
  // if bin has no events: add 1
  // safety feature to prevent KD = 1 as a
  // result of low statistics

  if(Pmzzm1m2==0){
    Pmzzm1m2++;
    }
  // - - - - - - - - - - - - - - - 

  float Pmzzcosthetastar = h_mzzcosthetastar->GetBinContent(h_mzzcosthetastar->FindBin(mZZ,costhetastar));
  float Pmzzcostheta2 = h_mzzcostheta2->GetBinContent(h_mzzcostheta2->FindBin(mZZ,costheta2));
  float Pmzzcostheta1 = h_mzzcostheta1->GetBinContent(h_mzzcostheta1->FindBin(mZZ,costheta1));
  float Pmzzphi1 = h_mzzphi1->GetBinContent(h_mzzphi1->FindBin(mZZ,phistar1));
  float Pmzzphi = h_mzzphi->GetBinContent(h_mzzphi->FindBin(mZZ,phi));

  //normalization
  float binwidth_mzzm1m2 = h_mzzm1m2->GetYaxis()->GetBinWidth(1) * h_mzzm1m2->GetZaxis()->GetBinWidth(1);
  float binwidth_mzzcosthetastar = h_mzzcosthetastar->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzcostheta1 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzcostheta2 = h_mzzcostheta1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzphi1 = h_mzzphi1->GetYaxis()->GetBinWidth(1);
  float binwidth_mzzphi = h_mzzphi->GetYaxis()->GetBinWidth(1);

  float Pmzzm1m2_norm = Pmzzm1m2/(n*binwidth_mzzm1m2); 
  float Pmzzcosthetastar_norm = Pmzzcosthetastar/(n*binwidth_mzzcosthetastar);
  float Pmzzcostheta1_norm = Pmzzcostheta1/(n*binwidth_mzzcostheta1);
  float Pmzzcostheta2_norm = Pmzzcostheta2/(n*binwidth_mzzcostheta2);
  float Pmzzphi1_norm = Pmzzphi1/(n*binwidth_mzzphi1);
  float Pmzzphi_norm = Pmzzphi/(n*binwidth_mzzphi);

  vector <float> P;
  P.push_back(Pmzzm1m2);
  P.push_back(Pmzzcosthetastar);
  P.push_back(Pmzzcostheta1);
  P.push_back(Pmzzcostheta2);
  P.push_back(Pmzzphi);
  P.push_back(Pmzzphi1);

  vector <float> P_norm;
  P_norm.push_back(Pmzzm1m2_norm);
  P_norm.push_back(Pmzzcosthetastar_norm);
  P_norm.push_back(Pmzzcostheta1_norm);
  P_norm.push_back(Pmzzcostheta2_norm);
  P_norm.push_back(Pmzzphi_norm);
  P_norm.push_back(Pmzzphi1_norm);

  /*delete h_mzz;
  delete h_mzzm1m2;
  delete h_mzzcosthetastar;
  delete h_mzzcostheta1;
  delete h_mzzcostheta2;
  delete h_mzzphi1;
  delete h_mzzphi;*/
  
  if(normalized)
    return P_norm;
  else
    return P;
}

pair<float,float> Mela::likelihoodDiscriminant (float mZZ, float m1, float m2, float costhetastar, float costheta1, float costheta2, float phi, float phistar1,
						int LHCsqrts, 
						bool withPt, float pt, 
						bool withY, float y){

  RooRealVar* z1mass_rrv = new RooRealVar("z1mass","m_{Z1}",0.,180.);
  RooRealVar* z2mass_rrv = new RooRealVar("z2mass","m_{Z2}",0.,120.); 
  RooRealVar* costhetastar_rrv = new RooRealVar("costhetastar","cos#theta^{*}",-1.,1.);  
  RooRealVar* costheta1_rrv = new RooRealVar("costheta1","cos#theta_{1}",-1.,1.);  
  RooRealVar* costheta2_rrv = new RooRealVar("costheta2","cos#theta_{2}",-1.,1.);
  RooRealVar* phi_rrv= new RooRealVar("phi","#Phi",-3.1415,3.1415);
  RooRealVar* phi1_rrv= new RooRealVar("phi1","#Phi_{1}",-3.1415,3.1415);
  RooRealVar* pt_rrv= new RooRealVar("pt","p_{T}^{4l}",0.,1000);
  RooRealVar* y_rrv= new RooRealVar("y","Y^{4l}",-4.0,4.0);
  RooRealVar* sqrtS_rrv= new RooRealVar("sqrtS","#sqrt{s}",1000,14000);
  RooRealVar* mzz_rrv= new RooRealVar("mzz","mZZ",80.,1000.);
  RooRealVar* upFrac_rrv= new RooRealVar("upFrac","upFrac",.5);
  upFrac_rrv->setConstant(kTRUE);

  AngularPdfFactory *SMHiggs = new AngularPdfFactory(z1mass_rrv,z2mass_rrv,costheta1_rrv,costheta2_rrv,phi_rrv,mzz_rrv);
  RooqqZZ_JHU_ZgammaZZ* SMZZ = new RooqqZZ_JHU_ZgammaZZ("SMZZ","SMZZ",*z1mass_rrv,*z2mass_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*costhetastar_rrv,*phi1_rrv,*mzz_rrv,*upFrac_rrv);
  SMHiggs->makeSMHiggs();
  SMHiggs->makeParamsConst(true);

  // shapes for pt and y ============

  static const int NptparamsS = 17;
  static const int NptparamsB = 11;

  string rrvnamesB[NptparamsB] = {"m","n0","n1","n2","ndue","bb0","bb1","bb2","T0","T1","T2"};
  RooRealVar *ptparamsB[NptparamsB];
  RooArgSet* allparamsB = new RooArgSet();
  for (int i = 0; i < NptparamsB; i++) {
    ptparamsB[i] = new RooRealVar(rrvnamesB[i].c_str(),rrvnamesB[i].c_str(),-10000.,10000.);
    allparamsB->add(*ptparamsB[i]);
  }

  string rrvnamesS[NptparamsS] = {"ms","ns0","ns1","ns2","ndues","bbs0","bbs1","bbs2","Ts0","Ts1","Ts2","bbdues0","bbdues1","bbdues2","fexps0","fexps1","fexps2"};
  RooRealVar *ptparamsS[NptparamsS];
  RooArgSet* allparamsS = new RooArgSet();
  for (int i = 0; i < NptparamsS; i++) {
    ptparamsS[i] = new RooRealVar(rrvnamesS[i].c_str(),rrvnamesS[i].c_str(),-10000.,10000.);
    allparamsS->add(*ptparamsS[i]);
  }
 
  char fileName[200];
  sprintf(fileName,"ZZMatrixElement/MELA/data/allParamsSig_%dTeV.txt",LHCsqrts);
  edm::FileInPath TsallisParams_Sig(fileName);
  string fullPath = TsallisParams_Sig.fullPath();

  RooTsallisExp* sigPt = new RooTsallisExp("sigPt","sigPt",*pt_rrv,*mzz_rrv,
					   *ptparamsS[0],*ptparamsS[1],*ptparamsS[2],
					   *ptparamsS[3],*ptparamsS[4],*ptparamsS[5],
					   *ptparamsS[6],*ptparamsS[7],*ptparamsS[8],
					   *ptparamsS[9],*ptparamsS[10],*ptparamsS[11],
					   *ptparamsS[12],*ptparamsS[13],*ptparamsS[14],
					   *ptparamsS[15],*ptparamsS[16]);

  allparamsS->readFromFile(fullPath.c_str(),0);

  sprintf(fileName,"ZZMatrixElement/MELA/data/allParamsBkg_%dTeV.txt",LHCsqrts);
  edm::FileInPath TsallisParam_Bkg(fileName);
  fullPath = TsallisParam_Bkg.fullPath();

  RooTsallis* bkgPt = new RooTsallis("bkgPt","bkgPt",*pt_rrv,*mzz_rrv,
				     *ptparamsB[0],*ptparamsB[1],*ptparamsB[2],
				     *ptparamsB[3],*ptparamsB[4],*ptparamsB[5],
				     *ptparamsB[6],*ptparamsB[7],*ptparamsB[8],
				     *ptparamsB[9],*ptparamsB[10]);
  allparamsB->readFromFile(fullPath.c_str(),0);

  RooRapiditySig* sigY = new RooRapiditySig("sigY", "sigY", *y_rrv, *mzz_rrv, *sqrtS_rrv);
  RooRapidityBkg* bkgY = new RooRapidityBkg("bkgY", "bkgY", *y_rrv, *mzz_rrv, *sqrtS_rrv);

  // ================================

  checkZorder(m1,m2,costhetastar,costheta1,costheta2,phi,phistar1);
  
  z1mass_rrv->setVal(m1);  
  z2mass_rrv->setVal(m2);
  costhetastar_rrv->setVal(costhetastar);
  costheta1_rrv->setVal(costheta1);
  costheta2_rrv->setVal(costheta2);
  phi_rrv->setVal(phi);
  phi1_rrv->setVal(phistar1);
  mzz_rrv->setVal(mZZ);
  if(withPt)
    pt_rrv->setVal(pt);
  if(withY)
    y_rrv->setVal(y);

  vector <float> P=my8DTemplate(1, mZZ,  m1,  m2,  costhetastar,  costheta1,  costheta2,  phi,  phistar1);
  
  float Pbackg=-99;
  float Psig=-99; 

  if(usePowhegTemplate_){
    // using template background calculation
    if(mZZ>100 && mZZ<180){
      Pbackg = P[0]*P[1]*P[2]*P[3]*P[4]*P[5]*5.0;
      Psig=SMHiggs->getVal(mZZ);
    }
    if(mZZ>180&&mZZ<=2*91.188){
      z1mass_rrv->setVal(mZZ/2.-1e-9);
      z2mass_rrv->setVal(mZZ/2.-1e-9);
      Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
      Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
    }
    if(mZZ>2*91.188){
      z1mass_rrv->setVal(91.188);
      z2mass_rrv->setVal(91.188);
      Pbackg = SMZZ->getVal()/(SMZZ->createIntegral(RooArgSet(*costhetastar_rrv,*costheta1_rrv,*costheta2_rrv,*phi_rrv,*phi1_rrv))->getVal())*10.0;
      Psig = SMHiggs->PDF->getVal()/(SMHiggs->PDF->createIntegral(RooArgSet(*costheta1_rrv,*costheta2_rrv,*phi_rrv))->getVal());
    }
  }else{

    // using analytic background calculation
    Pbackg = SMZZ->getVal()*1e-4/bkgPdfNorm(mZZ); 
    Psig = SMHiggs->PDF->getVal()/sigPdfNorm(mZZ);

  }

  if (withPt) {
    Pbackg *= bkgPt->getVal()/(bkgPt->createIntegral(RooArgSet(*pt_rrv))->getVal());
    Psig *= sigPt->getVal()/(sigPt->createIntegral(RooArgSet(*pt_rrv))->getVal());
  }
  if(withY) {
    Pbackg *= bkgY->getVal()/(bkgY->createIntegral(RooArgSet(*y_rrv))->getVal());
    Psig *= sigY->getVal()/(sigY->createIntegral(RooArgSet(*y_rrv))->getVal());
  }

  // - - - - - - - - - - - - - - - - - - - - - Whitbeck 
  // check whether P[i] is zero and print warning
  // message if so

  string varName[6]={"m1/m2","costhetastar","costheta1","coshteta2","phi","phistar1"};
  for(int iVar=0; iVar<6; iVar++){

    if(P[iVar]==0 && (m1+m2)<mZZ && m2>4 && mZZ>80 && mZZ<180)
	cout << " uh oh... Probability of " << varName[iVar] << " is zero." << endl;
  }
  // - - - - - - - - - - - - - - - - - - - - - 

  // the rrv's and pdf's should all be made members so that objects don't have to be created and deleted every call ~ AJW

  delete z1mass_rrv; 
  delete z2mass_rrv; 
  delete costhetastar_rrv;
  delete costheta1_rrv;
  delete costheta2_rrv;
  delete phi_rrv;
  delete phi1_rrv;
  delete y_rrv;
  delete pt_rrv;
  delete sqrtS_rrv;
  delete mzz_rrv; 
  delete SMZZ;
  delete SMHiggs;
  delete allparamsB;
  delete allparamsS;

  for(int i=0; i<NptparamsB; i++){
    delete ptparamsB[i]; 
  }
  for(int i=0; i<NptparamsS; i++){
    delete ptparamsS[i]; 
  }
  
  delete sigPt;
  delete bkgPt;
  delete sigY;
  delete bkgY;

  if(Psig<0 || Pbackg<0){
    cout<<"Mela::likelihoodDiscriminant() Error: KD not defined for this mzz (maybe mZZ<100 ?)"<<endl;
    cout << "=========================" << endl;
    cout << "psig: " << Psig << endl;
    cout << "pbkg: " << Pbackg << endl;
    cout << " - - - - - - - - - - - - " << endl;
    cout << "mzz: " << mZZ << endl;
    cout << "m1: " << m1 << endl;
    cout << "m2: " << m2 << endl;
    cout << "costheta1: " << costheta1 << endl;
    cout << "costheta2: " << costheta2 << endl;
    cout << "costhetastar: " << costhetastar << endl;
    cout << "phi: " << phi << endl;
    cout << "phi1: " << phistar1 << endl;
  }

  return make_pair(Psig,Pbackg);

}


// Re-order masses and angles as needed by likelihoodDiscriminant. 
// This follows a different convention than the usual Z1/Z2 definition!
void Mela::checkZorder(float& z1mass, float& z2mass,
		       float& costhetastar, float& costheta1,
		       float& costheta2, float& phi, 
		       float& phistar1){

  float tempZ1mass=z1mass;
  float tempZ2mass=z2mass;
  float tempH1=costheta1;
  float tempH2=costheta2;
  float tempHs=costhetastar;
  float tempPhi1=phistar1;
  float tempPhi=phi;

  if(z2mass>z1mass){
    //cout<<"inverted"<<endl;
    z1mass=tempZ2mass;
    z2mass=tempZ1mass;
    costhetastar=-tempHs;
    costheta1=tempH2;
    costheta2=tempH1;
    phi=tempPhi;
    phistar1=-tempPhi1-tempPhi;
    if(phistar1>3.1415)
      phistar1=phistar1-2*Geom::pi();
    if(phistar1<-3.1415)
      phistar1=phistar1+2*Geom::pi();

  }else
    return;

}
