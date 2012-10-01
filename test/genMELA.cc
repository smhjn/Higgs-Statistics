//-----------------------------------------------------------------------------
// File:        geninput.C
// Description: Build input histograms for use by likelihood script
//-----------------------------------------------------------------------------
// Created:     11-April-2011  Joe Bochenek
//-----------------------------------------------------------------------------

#include <time.h>
#include <sstream>
#include <vector>
#include <fstream>
#include <algorithm>
#include <TStyle.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <THStack.h>
#include <TMath.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include "pseudoMELA.h"
#include "MELA.h"

using namespace std;

#include "samples.h"
#include "Slurper.h"
#include "Systematics.h"
#include <ctype.h>
#include <stdio.h>

// apply cuts to the output events
int applycuts(std::map<std::string, double>& values){
	int notselected = 0;
	
	if(!( values["lept1_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept2_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept3_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept4_pfx"] < 0.4 ))  notselected++;
	
	// Cut to reduce pileup effect on MET
//	if(!((values["worst_iso_X"] + values["second_worst_iso_X"]) < 0.35 ) )  notselected++;
//	if(!( values["worst_vertex"] < 4))  notselected++;
//	if(!((values["mZ"] > 50) && (values["mZ"] < 120)))  notselected++;
//	if(!((values["mZstar"] > 12) && (values["mZstar"] < 120)))  notselected++;
//	if(!((values["mass4l"] > 112) && (values["mass4l"] < 162)))  notselected++;

	return notselected;
}




pair<double,double> calculateMELA(std::map<std::string, double>& values){

		TLorentzVector h;
		TLorentzVector z1;
		TLorentzVector z2;
		TLorentzVector l1;
		TLorentzVector l2;
		TLorentzVector l3;
		TLorentzVector l4;
		TLorentzVector z1_temp;


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
			z1_temp = z1;
		 	z1 = z2;
			z2 = z1_temp;
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

//	    pseudoMELA pseudoMELAcalculator;	    
 //       double pmelaval = pseudoMELAcalculator.eval(h.M(), z1.M(), z2.M(), costhetastar, costheta1, costheta2, phi, phi1);
double pmelaval = 0.;
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

		double melaval =  result.first/(result.first+result.second);;

		return make_pair(melaval,pmelaval);
}


int main(int argc, char* argv[]){

	clock_t start, end;
	start = clock();


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
	
	
	string dirSOURCE = "/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/";
	string dirIN = "/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/";
	
	//---------------------------------------------------------------------------
	// OUTPUT
	//---------------------------------------------------------------------------

	TString sig_cut_file =  dirIN + "sig_mela.dat";
	fstream sig_cuts;
	sig_cuts.open(sig_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << sig_cut_file << endl;
	return -1;
	}
//	sig_cuts << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString bkg_cut_file = dirIN + "bkg_mela.dat";
	fstream bkg_cuts;
	bkg_cuts.open(bkg_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << bkg_cut_file << endl;
	return -1;
	}
//	bkg_cuts << "weight/F:dphi2l/F:mt/F:m2l/F:deltar2l/F:pfmet/F:pt2/F:pt1/F:pt2l/F" << endl;


	TString dat_cut_file = dirIN + "data_mela.dat";
	fstream dat_cuts;
	dat_cuts.open(dat_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << dat_cut_file << endl;
	return -1;
	}

	TString zx_cut_file = dirIN + "zx_mela.dat";
	fstream zx_cuts;
	zx_cuts.open(zx_cut_file, ios_base::trunc | ios_base::out);
	if (!sig_cuts || !sig_cuts.good())
	{
	std::cout << "could not open file! " << zx_cut_file << endl;
	return -1;
	}


	// Print Headers

	TString dataString = dirSOURCE + sigfiles[0];
	Slurper file(dataString.Data());
//	vector<string> vars = file.nget();
	sig_cuts << "weight\t"; 
	bkg_cuts << "weight\t";
	dat_cuts << "weight\t";
	zx_cuts << "weight\t";

	for(int i= 0; i < vars.size(); i++ ){
		sig_cuts << vars[i] << "\t";
		bkg_cuts << vars[i] << "\t";
		dat_cuts << vars[i] << "\t";
		zx_cuts << vars[i] << "\t";

	}
	sig_cuts << "mela\tpmela\tsample\tindex\tchannel\tend" << endl; 
	bkg_cuts << "mela\tpmela\tsample\tindex\tchannel\tend" << endl;
	dat_cuts << "mela\tpmela\tchannel\tend" << endl;
	zx_cuts  << "mela\tpmela\tchannel\tend" << endl;
	
	const int maxNum = 3;
	
	TString filelist[maxNum] = {
		"data.dat",
		"bkg_bnn_test_uw.dat",
		"sig_bnn_train_uw.dat"
	};
	
	TString names[maxNum] = {
		"hdata",
		"hbkg",
		"hsig"
	};

	double weights[maxNum] = {	
		1.0,
		1.0,
		36.5721777888,
	};


	int ninputs = vars.size();
	float inputs1[ninputs];
	vector<double> data(ninputs);

	


	//---------------------------------------------------------------------------
	// Event Loop - DATA
	//---------------------------------------------------------------------------

	TString dataString3 = dirSOURCE + filelist[0];
	Slurper dat(dataString3.Data());

	int datcount = 0;
		
	// Loop over events
	while(dat.read())
	{
		if((dat.get("issamesign") )) continue;

		int channel = int(dat.get("channel"));
		datcount++;
		double weight = dat.get("weight");
		std::map<std::string, double> data;
		dat.mget(data);

		std::pair<double, double> melapair  = calculateMELA(data);
		double melaval = melapair.first;
		double pmelaval = melapair.second;

		dat_cuts << weight << "\t" << "\t";
		for(int i=0; i < vars.size(); i++) dat_cuts << dat.get(vars[i]) << "\t";
		dat_cuts << "\t" << melaval << "\t"<< pmelaval << "\t" << channel;
		dat_cuts << endl;
	}
		
	cout <<"data\t" << datcount <<  endl;

	dat.close();
	


	TString dataString1 = dirSOURCE + filelist[0];
	TFile f_fake("scripts/fake_reweight.root");
	f_fake.cd();
	TH2F *cratio_Run2012_e=(TH2F*) gDirectory->Get("cratio_Run2012_e");
	TH2F *cratio_Run2012_mu=(TH2F*) gDirectory->Get("cratio_Run2012_mu");


	//---------------------------------------------------------------------------
	// Event Loop - DATA opposite sign Z+X estimate
	//---------------------------------------------------------------------------

	TString dataString5 = dirSOURCE + filelist[0];
	Slurper dat2(dataString5.Data());

	int datcount2 = 0;
		
	// Loop over events
	while(dat2.read())
	{
		if(!(dat2.get("issamesign") )) continue;

		int channel = int(dat2.get("channel"));
		datcount2++;
		double weight = dat2.get("weight");
		std::map<std::string, double> data;
		dat2.mget(data);

		double reweight = 1.;
		if(channel==1){
		reweight *= cratio_Run2012_mu->Interpolate(float(dat2.get("lept3_eta")), float(dat2.get("lept3_pt"))) * cratio_Run2012_mu->Interpolate(float(dat2.get("lept4_eta")), float(dat2.get("lept4_pt")));
		}
		if(channel==2){
		reweight *= cratio_Run2012_e->Interpolate(float(dat2.get("lept3_eta")), float(dat2.get("lept3_pt"))) * cratio_Run2012_e->Interpolate(float(dat2.get("lept4_eta")), float(dat2.get("lept4_pt")));
		}		
		if(channel==3){
		reweight *= cratio_Run2012_mu->Interpolate(float(dat2.get("lept3_eta")), float(dat2.get("lept3_pt"))) * cratio_Run2012_mu->Interpolate(float(dat2.get("lept4_eta")), float(dat2.get("lept4_pt")));
		}	
		
		cout << "\trewegith: " << reweight << endl;


		std::pair<double, double> melapair  = calculateMELA(data);
		double melaval = melapair.first;
		double pmelaval = melapair.second;

		zx_cuts << reweight << "\t" << "\t";
		for(int i=0; i < vars.size(); i++) zx_cuts << dat2.get(vars[i]) << "\t";
		zx_cuts << "\t" << melaval << "\t" << pmelaval << "\t" << channel;
		zx_cuts << endl;
	}
	dat2.close();



	cout <<"z+x\t" << datcount2 <<  endl;



	//---------------------------------------------------------------------------
	// Event Loop - BACKGROUND
	//---------------------------------------------------------------------------

	TString dataString2 = dirIN + filelist[1];
	Slurper bkg(dataString2.Data());

	cout << dataString2 << endl;

	int bkgcount = 0;
		
	// Loop over events
	while(bkg.read())
	{
		int channel = int(bkg.get("channel"));
		int index = int(bkg.get("index"));

		bkgcount++;
		double weight = bkg.get("weight");
		std::map<std::string, double> data;
		bkg.mget(data);


		std::pair<double, double> melapair  = calculateMELA(data);
		double melaval = melapair.first;
		double pmelaval = melapair.second;
			
		double reweight = weights[1]/10000;
	
		bkg_cuts << reweight << "\t";
		for(int i=0; i < vars.size(); i++) bkg_cuts << bkg.get(vars[i]) << "\t";
		bkg_cuts << melaval << "\t" << pmelaval << "\t"  << index << "\t" << channel << endl;
	}
		
	bkg.close();
	
	cout <<"bkg\t" << bkgcount <<  endl;



	//---------------------------------------------------------------------------
	// Event Loop - SIGNAL
	//---------------------------------------------------------------------------

	TString dataString4 = dirIN + filelist[2];
	Slurper sig(dataString4.Data());

	cout << dataString4 << endl;
	
	int sigcount = 0;
		
	// Loop over events
	while(sig.read())
	{
		int channel = int(sig.get("channel"));
		int index = int(sig.get("index"));

		sigcount++;
		double weight = sig.get("weight");
		std::map<std::string, double> data;
		sig.mget(data);


		std::pair<double, double> melapair  = calculateMELA(data);
		double melaval = melapair.first;
		double pmelaval = melapair.second;
				
		double reweight = weights[1]/10000;
		
		sig_cuts << reweight << "\t";
		for(int i=0; i < vars.size(); i++) sig_cuts << sig.get(vars[i]) << "\t";
		sig_cuts << melaval << "\t" << pmelaval << "\t"  <<  index << "\t" << channel << endl;
	}	
	sig.close();
	
	cout <<"sig\t" << sigcount <<  endl;



	end = clock();
	
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";
}
