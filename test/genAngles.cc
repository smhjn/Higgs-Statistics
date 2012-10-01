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

	    pseudoMELA pseudoMELAcalculator;	    
        double pmelaval = pseudoMELAcalculator.eval(h.M(), z1.M(), z2.M(), costhetastar, costheta1, costheta2, phi, phi1);
//		double pmelaval = 0.;
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
//		melaval = 0.
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
	
	
	string dirSOURCE = "/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/";
	string dirIN = "/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/";

	//---------------------------------------------------------------------------
	// OUTPUT
	//---------------------------------------------------------------------------




	TString ps_file =  dirIN + "zx_test.dat";
	fstream sig_ps;
	sig_ps.open(ps_file, ios_base::trunc | ios_base::out);
	if (!sig_ps || !sig_ps.good())
	{
	std::cout << "could not open file! " << ps_file << endl;
	return -1;
	}
	

	sig_ps << "weight\tmass4l\tZ1mass\tZ2mass\tcosthetastar\tcostheta1\tcostheta2\tphi\tphi1\tpmela\tchannel" << endl;





	//---------------------------------------------------------------------------
	// Event Loop - SCALAR
	//---------------------------------------------------------------------------

	TString dataString1 = dirSOURCE + "zx_mela.dat";
	Slurper sig2(dataString1.Data());

	cout << dataString1 << endl;
	
	double sigcount = 0;
		
	// Loop over events
	while(sig2.read())
	{
	
		TLorentzVector h;
		TLorentzVector z1;
		TLorentzVector z2;
		TLorentzVector l1;
		TLorentzVector l2;
		TLorentzVector l3;
		TLorentzVector l4;
		TLorentzVector z1_temp;


		l1.SetPtEtaPhiM(sig2.get("lept1_pt"), sig2.get("lept1_eta"), sig2.get("lept1_phi"), 0.000511);
		l2.SetPtEtaPhiM(sig2.get("lept2_pt"), sig2.get("lept2_eta"), sig2.get("lept2_phi"), 0.000511);
		l3.SetPtEtaPhiM(sig2.get("lept3_pt"), sig2.get("lept3_eta"), sig2.get("lept3_phi"), 0.000511);
		l4.SetPtEtaPhiM(sig2.get("lept4_pt"), sig2.get("lept4_eta"), sig2.get("lept4_phi"), 0.000511);

		h = l1 + l2 + l3 + l4;
		z2 = l1 + l2;
		z1 = l3 + l4;
		

		int switched = 0;
		if(  fabs(z1.M() - sig2.get("Z1mass")) < 0.01 || fabs(z2.M() - sig2.get("Z2mass")) < 0.01 ) 
		{
//			cout << "not switching: " << z1.M() - sig2.get("Z1mass") << endl;
		}
		else
		{
//			cout << "switching: " << z1.M() - sig2.get("Z1mass") << endl;
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

		int channel = int(sig2.get("channel"));
		int index = int(sig2.get("index"));
		int weight = int(sig2.get("weight"));

		sigcount++;
		
		std::map<std::string, double> data;
		sig2.mget(data);

		std::pair<double, double> melapair  = calculateMELA(data);
		double melaval = melapair.first;
		double pmelaval = melapair.second;

		sig_ps << weight << "\t" <<  h.M() << "\t" << z1.M() << "\t" << z2.M() << "\t" << costhetastar << "\t" << costheta1 << "\t" << costheta2 << "\t" << phi << "\t" << phi1 << "\t" << pmelaval << "\t" << switched << endl;
	}	
	sig2.close();
	
	cout <<"sig\t" << sigcount <<  endl;



	end = clock();
	
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";
}


