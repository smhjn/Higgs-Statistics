//-----------------------------------------------------------------------------
// File:        weight.C
// Description: Do some preprocessing to the files, reweight for factors like
//				Data driven background estimates.  You can also apply a few simple
//				Cuts if needed
//-----------------------------------------------------------------------------
// Created:     08-Dec-2012 Joe Bochenek
//-----------------------------------------------------------------------------


#include <vector>
#include <iostream>
#include "Slurper.h"
#include "TLorentzVector.h"
#include "TPaveText.h"
#include "TH1.h"
#include "TCanvas.h"
#include <iomanip>
#include "MELA.h"


using namespace std;


int applycuts(std::map<std::string, double>& values, int channel){
	int notselected = 0;
//	cout << values["issamesign"] << endl;


//	if((int(values["issamesign"]) ))  notselected++;


	if(!( values["lept1_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept2_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept3_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept4_pfx"] < 0.4 ))  notselected++;
	
/*
	if(!( values["mass4l"] >110 ))  notselected++;
	if(!( values["mass4l"] <150 ))  notselected++;
	if(!( values["Z2mass"] >12 ))  notselected++;

	if(!( values["lept1_sip"] < 4 ))  notselected++;
	if(!( values["lept2_sip"] < 4 ))  notselected++;
	if(!( values["lept3_sip"] < 4 ))  notselected++;
	if(!( values["lept4_sip"] < 4 ))  notselected++;

	vector<double> pts;
	pts.push_back(values["lept1_pt"]);
	pts.push_back(values["lept2_pt"]);
	pts.push_back(values["lept3_pt"]);
	pts.push_back(values["lept4_pt"]);

	vector<double> etas;
	etas.push_back(values["lept1_eta"]);
	etas.push_back(values["lept2_eta"]);
	etas.push_back(values["lept3_eta"]);
	etas.push_back(values["lept4_eta"]);

	vector<double> phis;
	phis.push_back(values["lept1_phi"]);
	phis.push_back(values["lept2_phi"]);
	phis.push_back(values["lept3_phi"]);
	phis.push_back(values["lept4_phi"]);
	
	vector<double> mvas;
	mvas.push_back(values["lept1_mvaid"]);
	mvas.push_back(values["lept2_mvaid"]);
	mvas.push_back(values["lept3_mvaid"]);
	mvas.push_back(values["lept4_mvaid"]);
	
	
	
	// pt > 5 for all electrons
//	for(int i = 0; i < 4; ++i) if(pts[i] < 5) notselected++;


	// all mll candidates > 4GeV
	for(int i = 0; i < 4; ++i){
	for(int j = i + 1; j < 4; ++j){
	TLorentzVector pl1, pl2, p2l;
	pl1.SetPtEtaPhiM(pts[i], etas[i], phis[i], 0);
	pl2.SetPtEtaPhiM(pts[j], etas[j], phis[j], 0);
	p2l = pl1 + pl2;
	double mll = p2l.M();	
//	cout << etas[i] << "\t" << etas[j] << "\tmll: "  << mll << endl;
	if(!(mll>4)) notselected++;
	}
	}
	
	
	// Do the MVA electron selection for 4mu (AN 141)
	if(channel == 3){
	for(int i = 0; i < 2; ++i){
	if(pts[i] < 5) notselected++;
	if(pts[i] < 10){
	if((etas[i] < 0.8) && (mvas[i] < 0.47)) notselected++;
	if((etas[i] > 0.8) && (etas[i] < 1.479) && (mvas[i] < 0.004)) notselected++;
	if((etas[i] > 1.479) && (mvas[i] < 0.295)) notselected++;
	} else {
	if((etas[i] < 0.8) && (mvas[i] < 0.5)) notselected++;
	if((etas[i] > 0.8) && (etas[i] < 1.479) && (mvas[i] < 0.12)) notselected++;
	if((etas[i] > 1.479) && (mvas[i] < 0.6)) notselected++;	
	}
	}
	}

	// Do the MVA electron selection for 4mu (AN 141)
	if(channel == 2){
	for(int i = 0; i < 4; ++i){
	if(pts[i] < 5) notselected++;
	if(pts[i] < 10){
	if((etas[i] < 0.8) && (mvas[i] < 0.47)) notselected++;
	if((etas[i] > 0.8) && (etas[i] < 1.479) && (mvas[i] < 0.004)) notselected++;
	if((etas[i] > 1.479) && (mvas[i] < 0.295)) notselected++;
	} else {
	if((etas[i] < 0.8) && (mvas[i] < 0.5)) notselected++;
	if((etas[i] > 0.8) && (etas[i] < 1.479) && (mvas[i] < 0.12)) notselected++;
	if((etas[i] > 1.479) && (mvas[i] < 0.6)) notselected++;	
	}
	}
	}
	


	cout << 
	"\tlept1_pfx: " << values["lept1_pfx"] << 
	"\tlept2_pfx: " << values["lept2_pfx"] << 
	"\tlept3_pfx: " << values["lept3_pfx"] << 
	"\tlept4_pfx: " << values["lept4_pfx"] << 
	"\tlept1_pfx: " << values["lept1_sip"] << 
	"\tlept2_pfx: " << values["lept2_sip"] << 
	"\tlept3_pfx: " << values["lept3_sip"] << 
	"\tlept4_pfx: " << values["lept4_sip"] << 
	endl;
	

	if( RECOELE_PT[i] > 5. &&  RECOELE_PT[i] < 10. ){
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > .47 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
		&& RECOELE_mvaNonTrigV0[i] > .004 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > .295 ) BDT_ok = 1 ;
	}
	else { 
		if( fabs(RECOELE_scl_Eta[i]) < .8 && RECOELE_mvaNonTrigV0[i] > .5 ) BDT_ok = 1 ;
		if( ( fabs(RECOELE_scl_Eta[i]) >= .8 && fabs(RECOELE_scl_Eta[i]) <= 1.479 )
		&& RECOELE_mvaNonTrigV0[i] > .12 ) BDT_ok = 1 ;
		if( fabs(RECOELE_scl_Eta[i]) > 1.479 && RECOELE_mvaNonTrigV0[i] > .6 ) BDT_ok = 1 ;
	}
*/

	return notselected;

}

//-----------------------------------------------------------------------------
// MAIN MACRO
//-----------------------------------------------------------------------------


double calculateMELA(std::map<std::string, double>& values){

		TLorentzVector h;
		TLorentzVector z1;
		TLorentzVector z2;
		TLorentzVector l1;
		TLorentzVector l2;
		TLorentzVector l3;
		TLorentzVector l4;
		TLorentzVector l1_temp;
		TLorentzVector l2_temp;


		l1.SetPtEtaPhiM(values["lept1_pt"], values["lept1_eta"], values["lept1_phi"], 0.000511);
		l2.SetPtEtaPhiM(values["lept2_pt"], values["lept2_eta"], values["lept2_phi"], 0.000511);
		l3.SetPtEtaPhiM(values["lept3_pt"], values["lept3_eta"], values["lept3_phi"], 0.000511);
		l4.SetPtEtaPhiM(values["lept4_pt"], values["lept4_eta"], values["lept4_phi"], 0.000511);

		h = l1 + l2 + l3 + l4;
		z2 = l1 + l2;
		z1 = l3 + l4;
		

		int switched = 0;
		if(  fabs(z1.M() - values["Z1mass"]) < 0.01 || fabs(z2.M() - values["Z2mass"]) < 0.01 ) 
		{
//			cout << "not switching: " << z1.M() - values["Z1mass"] << endl;
		}
		else
		{
//			cout << "switching: " << z1.M() - values["Z1mass"] << endl;
			switched = 1;
			l1_temp = l1;
			l2_temp = l2;
		 	l1 = l3;
		 	l2 = l4;
			l3 = l1_temp;
			l4 = l2_temp;
		 }
		 

		
		double 		costheta1;
		double 		costheta2;
		double 		phi;
		double 		costhetastar;
		double 		phistar1;
		double 		phistar2;
		double 		phistar12;
		double 		phi1;
		double 		phi2;

		calculateAngles(h, z1, l1, l2, z2, l3, l4,  costheta1, costheta2, phi, costhetastar, phistar1, phistar2, phistar12, phi1, phi2);
		std::pair<double, double> result = likelihoodDiscriminant(h.M(), z1.M(), z2.M() ,  costhetastar,  costheta1,  costheta2,  phi,  phi1);


/*		
		TVector3 theZ1X_p3 = TVector3( z1.X(), z1.Y(), z1.Z() );
		TVector3 Z1dot = z1.Dot(l1);
		cout << "calculated values: " << endl;
		cout << "calculated values: " << endl;
		cout << "phi1: " << theZ1X_p3.Phi() << endl;
		cout << "phi2: " << z2.Phi() << endl;
		

		cout << "costheta1: " << costheta1 << endl;
		cout << "costheta2: " << costheta2 << endl;
		cout << "phi: " << phi << endl;
		cout << "phi1: " << phi1 << endl;
		cout << "costhetastar: " << costhetastar << endl;


		if(switched){

		cout << "h: " << values["mass4l"] << endl;

		cout << "Z1: " << values["Z1mass"] << endl;
		cout << "Z1: " << z1.M() << endl;
		
		cout << "Z2: " << values["Z2mass"] << endl;
		cout << "Z2: " << z2.M() << endl;
		
		cout << "l1: " << l1.Pt() << endl;
		cout << "l2: " << l2.Pt() << endl;
		cout << "l3: " << l3.Pt() << endl;
		cout << "l4: " << l4.Pt() << endl;

		cout << "l1: " << values["lept1_pt"] << endl;
		cout << "l2: " << values["lept2_pt"] << endl;
		cout << "l3: " << values["lept3_pt"] << endl;
		cout << "l4: " << values["lept4_pt"] << endl;
		cout << "first: " << result.first << endl;
		cout << "second: " << result.second << endl;
		}
	*/

		return  result.first/(result.first+result.second);;

}

int main()
{

	
	const int maxbkg = 19;
	const int maxdat = 16;
	const int maxsig = 16;
	string channel = "1";
	int channel_ = 1;
		
	string dirSOURCE = "/home/jbochenek/data/HZZ4l_2012_data/lustre/cms/store/user/jpb/data/";
	string dirIN = "/home/jbochenek/data/HZZ4l_2012_data/hzz_trial/";

	TString sigfiles[maxsig] = {
	"GluGluToHToZZTo4L_M-115_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-117_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-119_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-120_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-121_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-123_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-124_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-125_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-126_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-127_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-128_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-130_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-135_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-140_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-145_8TeV-powheg-pythia6_4mu_bnn.txt",
	"GluGluToHToZZTo4L_M-150_8TeV-powheg-pythia6_4mu_bnn.txt"
	};

	int sigNames[maxsig] = {
	115,
	117,
	119,
	120,
	121,
	123,
	124,
	125,
	126,
	127,
	128,
	130,
	135,
	140,
	145,
	150,
	};
	
	TString datfiles[maxdat] = {
	"DoubleElectron_Run2012A_Jun08ReReco_noLowPileUp_noGoldetJsonJun22_noMay23_4mu_bnn.txt",
	"DoubleElectron_Run2012A_May23ReReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"DoubleElectron_Run2012A_PromptReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"DoubleElectron_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"DoubleElectron_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_uptoJun24_190456-196531_noLowPU_4mu_bnn.txt",
	"DoubleMu_Run2012A_Jun08ReReco_noLowPileUp_noGoldetJsonJun22_noMay23_4mu_bnn.txt",
	"DoubleMu_Run2012A_May23ReReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"DoubleMu_Run2012A_PromptReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"DoubleMu_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_1_4mu_bnn.txt",
	"DoubleMu_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_2_4mu_bnn.txt",
	"DoubleMu_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_uptoJun24_190456-196531_noLowPU_4mu_bnn.txt",
	"MuEG_Run2012A_PromptReco_Jun22_190456-196509_noLowPU_4mu_bnn.txt",
	"MuEG_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_1_4mu_bnn.txt",
	"MuEG_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_2_4mu_bnn.txt",
	"MuEG_Run2012B_PromptReco_Jun22_190456-196509_noLowPU_uptoJun24_190456-196531_noLowPU_4mu_bnn.txt",
	};

	TString bkgfiles[maxbkg] = {
	"DYJetsToLL_M-50_TuneZ2Star_8TeV-madgraph-tarball_4mu_bnn.txt",
	"Tbar_s-channel_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"Tbar_t-channel_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"Tbar_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"T_s-channel_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"T_t-channel_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"T_tW-channel-DR_TuneZ2star_8TeV-powheg-tauola_4mu_bnn.txt",
	"TTJets_TuneZ2star_8TeV-madgraph-tauola_4mu_bnn.txt",
	"WZJetsTo3LNu_TuneZ2_8TeV-madgraph-tauola_4mu_bnn.txt",
//	"WZTo3LNu_TuneZ2star_8TeV_pythia6_tauola_4mu_bnn.txt",
	"GluGluToZZTo2L2L_TuneZ2star_8TeV-gg2zz-pythia6_4mu_bnn.txt",
	"GluGluToZZTo4L_8TeV-gg2zz-pythia6_4mu_bnn.txt",
	"WWJetsTo2L2Nu_TuneZ2star_8TeV-madgraph-tauola_4mu_bnn.txt",
	"WWTo2L2Nu_TuneZ2star_8TeV_pythia6_tauola_4mu_bnn.txt",
	"ZZTo2mu2tau_8TeV-powheg-pythia6_4mu_bnn.txt",
	"ZZTo2e2tau_8TeV-powheg-pythia6_4mu_bnn.txt",
	"ZZTo4tau_8TeV-powheg-pythia6_4mu_bnn.txt",
	"ZZTo4mu_8TeV-powheg-pythia6_4mu_bnn.txt",
	"ZZTo2e2mu_8TeV-powheg-pythia6_4mu_bnn.txt",
	"ZZTo4e_8TeV-powheg-pythia6_4mu_bnn.txt",
//	"ZZTo4L_TuneZ2star_8TeV_pythia6_tauola_4mu_bnn.txt",
	};

       
	TString bkgNames[maxbkg] = {
	"Z+Jets",
	"top",
	"top",
	"top",
	"top",
	"top",
	"top",
	"top",
	"WZ",
	"ZZ",
	"ZZ",
	"WW",
	"WW",
	"ZZ",
	"ZZ",
	"ZZ",
	"ZZ",
	"ZZ",
	"ZZ",
	};

	int bkgIndex[maxbkg] = {
	1,
	2,
	2,
	2,
	2,
	2,
	2,
	2,
	3,
	4,
	4,
	5,
	5,
	4,
	4,
	4,
	4,
	4,
	4
	};
	
	TString sampleNames_sig = "sig_ggH";
	
	TString sampleNames[maxbkg] = {
	"bkg_zjets",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"NOSYST",
	"bkg_ggzz",
	"bkg_ggzz"
	"NOSYST",
	"NOSYST",
	"bkg_qqzz",
	"bkg_qqzz",
	"bkg_qqzz",
	"bkg_qqzz",
	"bkg_qqzz",
	"bkg_qqzz",
	};


	TString log_file = "weight_log.dat";
	fstream log;
	log.open(log_file, ios_base::trunc | ios_base::out);
	if (!log || !log.good())
	{
	std::cout << "could not open file! " << log_file << endl;
	return -1;
	}

	//---------------------------------------------------------------------------
	// Set up TMVA NN
	//---------------------------------------------------------------------------

	vector<string> vars;
	string line;
	ifstream varsfile ("scripts/vars.txt");
	if (varsfile.is_open())
	{
	cout << "Variables: " << endl;
	while ( varsfile.good() )
	{
	getline (varsfile,line);
	cout << line << endl;
	vars.push_back(line);
	}
	varsfile.close();
	}

    else 
    {
	cout << "Unable to open mvavars.txt file"; 
	return 0;
	}

	//---------------------------------------------------------------------------
	// OUTPUT
	//---------------------------------------------------------------------------

	TString sig_cut_file =  dirIN + "sig_" + channel + ".dat";
	fstream sig_cuts;
	sig_cuts.open(sig_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << sig_cut_file << endl;
	return -1;
	}
//	sig_cuts << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString bkg_cut_file = dirIN + "bkg_" + channel + ".dat";
	fstream bkg_cuts;
	bkg_cuts.open(bkg_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << bkg_cut_file << endl;
	return -1;
	}
//	bkg_cuts << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString dat_cut_file = dirIN + "data_" + channel + ".dat";
	fstream dat_cuts;
	dat_cuts.open(dat_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << dat_cut_file << endl;
	return -1;
	}

	// Print Headers

	TString dataString = dirSOURCE + sigfiles[0];
	Slurper file(dataString.Data());
//	vector<string> vars = file.nget();
	sig_cuts << "weight\t"; 
	bkg_cuts << "weight\t";
	dat_cuts << "weight\t";
	for(int i= 0; i < vars.size(); i++ ){
		sig_cuts << vars[i] << "\t";
		bkg_cuts << vars[i] << "\t";
		dat_cuts << vars[i] << "\t";
	}
	sig_cuts << "mela\tsample\tindex\tchannel\tend" << endl; 
	bkg_cuts << "mela\tsample\tindex\tchannel\tend" << endl;
	dat_cuts << "mela\tchannel\tend" << endl;


	//---------------------------------------------------------------------------
	// TMVA OUTPUT
	//---------------------------------------------------------------------------

	TString sig_cut_file_tmva =  dirIN + "sig_tmva.dat";
	fstream sig_cuts_tmva;
	sig_cuts_tmva.open(sig_cut_file_tmva, ios_base::trunc | ios_base::out);
	if (!sig_cuts_tmva || !sig_cuts_tmva.good())
	{
	std::cout << "could not open file! " << sig_cut_file_tmva << endl;
	return -1;
	}


	TString bkg_cut_file_tmva = dirIN + "bkg_tmva.dat";
	fstream bkg_cuts_tmva;
	bkg_cuts_tmva.open(bkg_cut_file_tmva, ios_base::trunc | ios_base::out);
	if (!bkg_cuts_tmva || !bkg_cuts_tmva.good())
	{
	std::cout << "could not open file! " << bkg_cut_file_tmva << endl;
	return -1;
	}

	// Print Headers
	sig_cuts_tmva << "weight/F:"; 
	bkg_cuts_tmva << "weight/F:";
	sig_cuts_tmva << "index/I:"; 
	bkg_cuts_tmva << "index/I:";
	for(int i= 0 ; i < vars.size(); i++ ){
		sig_cuts_tmva << vars[i] << "/F:";
		bkg_cuts_tmva << vars[i] << "/F:";
	}

	sig_cuts_tmva << endl;
	bkg_cuts_tmva << endl;

	
	//---------------------------------------------------------------------------
	// BNN OUTPUT
	//---------------------------------------------------------------------------

	TString sig_bnn_train_file =  dirIN + "sig_bnn_train" + channel + ".dat";
	fstream sig_bnn_train;
	sig_bnn_train.open(sig_bnn_train_file, ios_base::trunc | ios_base::out);
	if (!sig_bnn_train || !sig_bnn_train.good())
	{
	std::cout << "could not open file! " << sig_bnn_train_file << endl;
	return -1;
	}
//	sig_bnn_train << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString bkg_bnn_train_file = dirIN + "bkg_bnn_train" + channel + ".dat";
	fstream bkg_bnn_train;
	bkg_bnn_train.open(bkg_bnn_train_file, ios_base::trunc | ios_base::out);
	if (!sig_bnn_train || !sig_bnn_train.good())
	{
	std::cout << "could not open file! " << bkg_bnn_train_file << endl;
	return -1;
	}

	// Print Headers

	sig_bnn_train << "weight\t"; 
	bkg_bnn_train << "weight\t";
	for(int i= 0 ; i < vars.size(); i++ ){
		sig_bnn_train << vars[i] << "\t";
		bkg_bnn_train << vars[i] << "\t";
	}
	sig_bnn_train << "mela\tindex\tchannel\tsample" << endl; 
	bkg_bnn_train << "mela\tindex\tchannel\tsample" << endl;

	TString sig_bnn_test_file =  dirIN + "sig_bnn_test" + channel + ".dat";
	fstream sig_bnn_test;
	sig_bnn_test.open(sig_bnn_test_file, ios_base::trunc | ios_base::out);
	if (!sig_bnn_test || !sig_bnn_test.good())
	{
	std::cout << "could not open file! " << sig_bnn_test_file << endl;
	return -1;
	}
//	sig_bnn_test << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString bkg_bnn_test_file = dirIN + "bkg_bnn_test" + channel + ".dat";
	fstream bkg_bnn_test;
	bkg_bnn_test.open(bkg_bnn_test_file, ios_base::trunc | ios_base::out);
	if (!sig_bnn_test || !sig_bnn_test.good())
	{
	std::cout << "could not open file! " << bkg_bnn_test_file << endl;
	return -1;
	}
//	bkg_bnn_test << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	// Print Headers

	sig_bnn_test << "weight\t"; 
	bkg_bnn_test << "weight\t";
	for(int i= 0 ; i < vars.size(); i++ ){
		sig_bnn_test << vars[i] << "\t";
		bkg_bnn_test << vars[i] << "\t";
	}
	sig_bnn_test << "mela\tindex\tchannel\tsample" << endl; 
	bkg_bnn_test << "mela\tindex\tchannel\tsample" << endl;








	//---------------------------------------------------------------------------
	// BKG
	//---------------------------------------------------------------------------
	// Apply cuts and reweight and put all bkgnals in to one file
	for(int j = 0; j < maxbkg; j++){
	int bkgout = 0;
	int bkgin = 0;
	double bkgin_w = 0;
	double bkgout_w = 0;

	TString dataString = dirSOURCE + bkgfiles[j];
	cout << "\n Processing file " << j << " " << dataString << endl;
	Slurper bkg(dataString.Data());
	while ( bkg.read() )
	{
	double weight = bkg.get("weight");
	std::map<std::string, double> data;
	bkg.mget(data);
	bkgin++;
	bkgin_w+=weight;
	if (applycuts(data, channel_)) continue;
//	cout << "\n Processing file " << j << " " << dataString << endl;

	double melaval = calculateMELA(data);

	bkg_cuts << weight << "\t";
	for(int i=0; i < vars.size(); i++) bkg_cuts << bkg.get(vars[i]) << "\t";
	bkg_cuts << melaval << "\t" << bkgIndex[j] << "\t" << channel << endl;

	bkg_cuts_tmva << weight << "\t"<<  j << "\t" << sampleNames[j] << "\t";
	for(int i=0; i < vars.size(); i++) bkg_cuts_tmva << bkg.get(vars[i]) << "\t";
	bkg_cuts_tmva << endl;

	if(bkgIndex[j] == 4){
	bkgout++;
	if(bkgout%2){
	bkg_bnn_train << weight << "\t";
	for(int i=0; i < vars.size(); i++) bkg_bnn_train << bkg.get(vars[i]) << "\t";
	bkg_bnn_train << melaval << "\t" << bkgIndex[j] << "\t" << channel << endl;
	} else {
	bkg_bnn_test << weight << "\t";
	for(int i=0; i < vars.size(); i++) bkg_bnn_test << bkg.get(vars[i]) << "\t";
	bkg_bnn_test << melaval << "\t" << bkgIndex[j] << "\t" << channel << endl;
	}
	}

	bkgout_w+=weight;
	}	
	bkg.close();
	cout << "bkg (" << bkgNames[j] << ")\t " << bkgout << "/" << bkgin << "\t" << bkgout_w << "/" << bkgin_w << endl;
	log << "bkg (" << bkgNames[j] << ")\t " << bkgout << "/" << bkgin << "\t" << bkgout_w << "/" << bkgin_w << endl;
	}
	cout << endl;



	//---------------------------------------------------------------------------
	// DATA
	//---------------------------------------------------------------------------
	// Apply cuts and reweight and put all details in to one file
	for(int j = 0; j < maxdat; j++){
	int datout = 0;
	int datin = 0;
	double datin_w = 0;
	double datout_w = 0;

	TString dataString = dirSOURCE + datfiles[j];
	cout << "\n Processing file " << dataString << endl;
	Slurper dat(dataString.Data());
	while ( dat.read() )
	{
	double weight = 1;
	std::map<std::string, double> data;
	dat.mget(data);
//	for ( std::map<std::string, double>::const_iterator iter = data.begin();iter != data.end(); ++iter )	cout << iter->first << '\t' << iter->second << '\t';
//	vector<string> vals = file.nget();
	datin++;
	datin_w+=weight;
	if (applycuts(data, channel_)) continue;
	datout++;	

	double melaval = calculateMELA(data);

	dat_cuts << weight << "\t" << "\t";
	for(int i=0; i < vars.size(); i++) dat_cuts << dat.get(vars[i]) << "\t";
	dat_cuts << "\t" << melaval << "\t" << channel;
	dat_cuts << endl;
	datout_w+=weight;
	}	
	dat.close();
	cout << "dat\t " << datout << "/" << datin << "\t" << datout_w << "/" << datin_w << endl;
	log << "dat\t " << datout << "/" << datin << "\t" << datout_w << "/" << datin_w << endl;
	}
	cout << endl;


	//---------------------------------------------------------------------------
	// SIGNAL
	//---------------------------------------------------------------------------
	// Apply cuts and reweight and put all signals in to one file


	TH1F* hsig_mela[maxsig];
	
	
	for(int j = 0; j < maxsig; j++){
	char title[256];
	sprintf(title, "hsig_mela_%d", j);
	
	hsig_mela[j] = new TH1F(title,title, 25, 0, 1);
	hsig_mela[j]->SetDirectory(0);
	hsig_mela[j]->Sumw2();

	int sigout = 0;
	int sigin = 0;
	double sigin_w = 0;
	double sigout_w = 0;
	double average = 0.;
	int numave = 0;

	TString dataString = dirSOURCE + sigfiles[j];
	cout << "\n Processing file " << dataString << endl;
	Slurper sig(dataString.Data());
	while ( sig.read() )
	{
	double weight = sig.get("weight");
	std::map<std::string, double> data;
	sig.mget(data);
	if(!sigin) cout << sig.get("mass4l") << endl;
//	for ( std::map<std::string, double>::const_iterator iter = data.begin();iter != data.end(); ++iter )	cout << iter->first << '\t' << iter->second << '\t';
//	vector<string> vals = file.nget();
	sigin++;
	sigin_w+=weight;
	if (applycuts(data, channel_)) continue;
	sigout++;	


	double melaval = calculateMELA(data);
	if(melaval>0) average +=melaval;	

	hsig_mela[j]->Fill(melaval);

//	cout << sig.get("mass4l") << "\t" <<  melaval << "\t" << average / numave << endl;
	
	
	sig_cuts << weight << "\t";
	for(int i=0; i < vars.size(); i++) sig_cuts << sig.get(vars[i]) << "\t";
	sig_cuts << melaval << "\t" <<  j << "\t" << channel << endl;

	if(sigout%2){
	sig_bnn_train << weight << "\t";
	for(int i=0; i < vars.size(); i++) sig_bnn_train << sig.get(vars[i]) << "\t";
	sig_bnn_train << melaval << "\t" <<  j << "\t" << channel << endl;
	} else {
	sig_bnn_test << weight << "\t";
	for(int i=0; i < vars.size(); i++) sig_bnn_test << sig.get(vars[i]) << "\t";
	sig_bnn_test << melaval << "\t" <<  j << "\t" << channel << endl;
	}
	
	sig_cuts_tmva << weight << "\t" <<  j  << "\t" << sampleNames_sig << "\t";
	for(int i=0; i < vars.size(); i++) sig_cuts_tmva << sig.get(vars[i]) << "\t";
	sig_cuts_tmva << endl;
	
	sigout_w+=weight;
	}	
	sig.close();
	cout << "Signal (" << sigNames[j] << ")\t " << sigout << "/" << sigin << "\t" << sigout_w << "/" << sigin_w << endl;
	log << "Signal (" << sigNames[j] << ")\t " << sigout << "/" << sigin << "\t" << sigout_w << "/" << sigin_w << endl;
	char filename[256];
	sprintf(filename, "output/hsig_mela_%d.gif", j);

	TCanvas *c1 = new TCanvas(title, title, 800, 960);
	c1->cd();
	hsig_mela[j]->Draw("hist");
	c1->SaveAs(filename);
	cout << "average" << average/sigout << endl;
	}
	cout << endl;



	
	sig_cuts.close();
	bkg_cuts.close();
	dat_cuts.close();
	sig_cuts_tmva.close();
	bkg_cuts_tmva.close();

}
