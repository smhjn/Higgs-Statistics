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
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <THStack.h>
#include <TMath.h>
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TPaveText.h"
#include <iomanip>
#include <TRandom.h>
#include <TGaxis.h>
#include <TGraph2D.h>

using namespace std;

#include "samples.h"
#include "Slurper.h"
#include "Systematics.h"
#include <ctype.h>
#include <stdio.h>

// apply cuts to the output events
int applycuts(std::map<std::string, double>& values, int channel){
	int notselected = 0;
//	cout << values["issamesign"] << endl;
	if((int(values["issamesign"]) ))  notselected++;


	if(!( values["lept1_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept2_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept3_pfx"] < 0.4 ))  notselected++;
	if(!( values["lept4_pfx"] < 0.4 ))  notselected++;
	

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
	
	
	// Do the MVA electron selection for 2e2mu (AN 141)
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

	// Do the MVA electron selection for 4e (AN 141)
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
	

/*
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

int main(int argc, char* argv[]){

	clock_t start, end;
	start = clock();

	if (argc < 3) {
			std::cerr << "Usage: " << argv[0] << "-r RANDOM_NUMBER_SEED (integer)" << std::endl;
			return 1;
	}

	int random_seed = 3;

	if (std::string(argv[1]) == "-r") {
			if (2 < argc) { // Make sure we aren't at the end of argv!


				if (int seed = atoi(argv[2])) { // Make sure we aren't at the end of argv!
					cout << "Random seed set to " <<  seed << endl;
					random_seed = int(seed);
					cout << "argv: " << argv[0] << "\t" << argv[1] << "\t" << seed << endl; 
				} else {
				  std::cerr << "-r option requires an integer." << std::endl;	
				  std::cerr << "Usage: " << argv[0] << "-r RANDOM_NUMBER_SEED (integer)" << std::endl;
				  return 1;			
				}
			} else { // Uh-oh, there was no argument
				  std::cerr << "-r option requires one argument." << std::endl;
        		  std::cerr << "Usage: " << argv[0] << "-r RANDOM_NUMBER_SEED (integer)" << std::endl;
				  return 1;
			}  
	} else {
			std::cerr << "Usage: " << argv[0] << "-r RANDOM_NUMBER_SEED (integer)" << std::endl;
			return 1;		
	}

	std::cout << argv[1] << std::endl;
        
	//TString html_name = "/afs/cern.ch/user/j/jpb/www/studies/HWW_analysis/plots/index.html";  
	//ofstream event_list;
	//event_list.open(html_name, ios_base::trunc);
	
	TFile f2("output/hist_mela.root", "RECREATE");
	const int thischannel = 0;

	//---------------------------------------------------------------------------
	// Set up TMVA NN
	//---------------------------------------------------------------------------

	int ninputs1 = 0, ninputs2 = 0, ninputs3 = 0;
	vector<string> var1;
	string line;
	ifstream vars1 ("scripts/mvavars1.txt");
	if (vars1.is_open())
	{
	cout << "Variables: " << endl;
	while ( vars1.good() )
	{
	getline (vars1,line);
	cout << line << endl;
	var1.push_back(line);
	}
	vars1.close();
	 ninputs1 = var1.size();
	}

    else 
    {
    	cout << "Unable to open mvavars.txt file"; 
		return 0;
	}


	vector<string> var2;
	ifstream vars2 ("scripts/mvavars2.txt");
	if (vars2.is_open())
	{
	cout << "Variables: " << endl;
	while ( vars2.good() )
	{
	getline (vars2,line);
	cout << line << endl;
	var2.push_back(line);
	}
	vars2.close();
	 ninputs2 = var2.size();
	}

    else 
    {
    	cout << "Unable to open mvavars.txt file"; 
		return 0;
	}

	vector<string> var3;
	ifstream vars3 ("scripts/mvavars3.txt");
	if (vars3.is_open())
	{
	cout << "Variables: " << endl;
	while ( vars3.good() )
	{
	getline (vars3,line);
	cout << line << endl;
	var3.push_back(line);
	}
	vars3.close();
	 ninputs3 = var3.size();
	}

    else 
    {
    	cout << "Unable to open mvavars.txt file"; 
		return 0;
	}
	

//    mvavars.pop_back();

	//---------------------------------------------------------------------------
	// Initialize
	//---------------------------------------------------------------------------

	double mvacut[3] =
	{
	0.512,
	0.428,
	0.836
	};

	const int maxNum = 3;
	
	TString filelist[maxNum] = {
		"data.dat",
		"bkg.dat",
		"sig.dat"
	};
	
	TString names[maxNum] = {
		"hdata",
		"hbkg",
		"hsig"
	};


	int ninputs = var2.size();
	float inputs1[ninputs];
	vector<double> data(ninputs);

	
	

	//---------------------------------------------------------------------------
	// BOOK HISTOGRAMS
	//---------------------------------------------------------------------------

//	TH2F* histos_bkg[maxNum];
	TH1F* hmlp_s[maxsig_plots];
	TH1F* histos_sig_bins[maxsig_plots];

	TH1F* hmlp_b;

	TH2F* histos_bkg;
	TH2F* histos_bkg2;
	TH2F* histos_bkg3;
	TH2F* histos_bkg4;

	TH2F* histos_dat;
	TH2F* histos_sig[maxsig_plots];
	TH2F* histos_bkgs[maxsig_plots];


	TH1F* hmlp_weighted[maxNum];

	hmlp_b =new TH1F("hbkg", "hbkg", xbins, xmin, xmax);
	hmlp_b->SetMarkerColor(1);
	hmlp_b->SetLineColor(1);
	hmlp_b->Sumw2();

	histos_bkg = new TH2F("bkg2d", "bkg2d",  ybins, ymin, ymax, xbins, xmin, xmax);
	histos_bkg->SetDirectory(0);
	histos_bkg->Sumw2();
	
	histos_dat = new TH2F("dat2d", "dat2d", ybins, ymin, ymax, xbins, xmin, xmax);
	histos_dat->SetDirectory(0);
	histos_dat->Sumw2();

	histos_bkg2 = new TH2F("zjets2d", "zjets2d", ybins, ymin, ymax, xbins, xmin, xmax);
	histos_bkg2->SetDirectory(0);
	histos_bkg2->Sumw2();

	histos_bkg3 = new TH2F("zjets2d_sel", "zjets2d_sel", ybins, ymin, ymax, xbins, xmin, xmax);
	histos_bkg3->SetDirectory(0);
	histos_bkg3->Sumw2();

	histos_bkg4 = new TH2F("all_redecible_bkg", "all_redecible_bkg", ybins, ymin, ymax, xbins, xmin, xmax);
	histos_bkg4->SetDirectory(0);
	histos_bkg4->Sumw2();


	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_" << sigNames_mc[k]; 
		TString name = _name.str();

		histos_sig[k] = new TH2F(name, name, ybins, ymin, ymax, xbins, xmin, xmax);
		histos_sig[k]->Sumw2();
	}

	for(Int_t k=0; k < maxbkg; k++) 
	{
		std::stringstream _name;
		_name << "histos_bkg_" << k; 
		TString name = _name.str();

		histos_bkgs[k] = new TH2F(name, name, ybins, ymin, ymax, xbins, xmin, xmax);
		histos_bkgs[k]->Sumw2();
	}

	TH2F* histos_sig_morph[maxsig_plots];
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_morph_" << sigNames_mc[k]; 
		TString name = _name.str();

		histos_sig_morph[k] = new TH2F(name, name, ybins_morph, ymin_morph, ymax_morph, xbins_morph, xmin, xmax);
		histos_sig_morph[k]->Sumw2();
	}


	// Set up systematics
	// Systematics getsystematics(random_seed, mvavars);

	int k = 0;




	//---------------------------------------------------------------------------
	// Event Loop - DATA
	//---------------------------------------------------------------------------

	TString dataString3 = dirIN + filelist[0];
	Slurper dat(dataString3.Data());
	int datcount = 0;
		
	// Loop over events
	while(dat.read())
	{
		int channel = int(dat.get("channel"));
		if(thischannel!=0) if(!(channel==thischannel)) continue;
		if((dat.get("issamesign") )) continue;
		datcount++;
		double weight = dat.get("weight");

		double y = dat.get("mela");	

//		for(int i=0; i < ninputs; i++) inputs[i] = dat.get(mvavars[i]);
//		double y2 = BNNValue->GetMvaValue(inputs);	

		cout << "channel: " << channel << "\t mva:" << y << "\t mass: " << dat.get("mass4l") << endl;

//		for(int i=0; i < ninputs2; i++) cout << var2[i] << ": " << data[i] << "\t";
		std::map<std::string, double> data;
		dat.mget(data);
		if (applycuts(data, channel)) continue;

//		if(y > 0.6	) histos_dat->Fill(y, dat.get("mass4l"), weight);
		histos_dat->Fill(dat.get("mass4l"), y,  weight);
		cout << "mva: " << y << endl;

	}
		
	histos_dat->Write();
	dat.close();
	
	cout <<"data\t" << histos_dat->Integral() << "/" << datcount <<  endl;

//	exit(0);
	//---------------------------------------------------------------------------
	// Event Loop - Z+Jets Bkg
	//---------------------------------------------------------------------------

	// use same sign data events
//	TString dataString1 = dirIN + "ZJets.txt";
	TString dataString1 = dirIN + filelist[0];

	Slurper bkg2(dataString1.Data());
	int bkg2count = 0;
		
	// Loop over events
	while(bkg2.read())
	{
		int channel = int(bkg2.get("channel"));
		if(!(bkg2.get("issamesign") )) continue;
		
//		if(thischannel!=0) if(!(channel==thischannel)) continue;
		bkg2count++;
		double weight = bkg2.get("weight");
		double y = bkg2.get("mela");	
		cout << "channel: " << channel << "\t mva:" << y << "\t mass: " << bkg2.get("mass4l") << endl;

		histos_bkg2->Fill(bkg2.get("mass4l"), y,  1.0);

		std::map<std::string, double> data;
		bkg2.mget(data);
		if (applycuts(data, channel)) continue;

		cout << "selected" <<endl;

		histos_bkg3->Fill(bkg2.get("mass4l"), y,  1.0);
//		if(y > mvacut) histos_bkg2->Fill(y, bkg2.get("mass4l"), weight);
//		if(y > mvacut) histos_bkg3->Fill(y, bkg2.get("mass4l"), weight);
	}
		
	histos_bkg2->Write();
	histos_bkg3->Write();

	bkg2.close();
	
	cout <<"data\t" << histos_dat->Integral() << "/" << datcount <<  endl;
	
	

	//---------------------------------------------------------------------------
	// Event Loop - BACKGROUND
	//---------------------------------------------------------------------------

	TString dataString = dirIN + filelist[1];
	Slurper bkg(dataString.Data());

	cout << "Processing file: " << dataString << endl;
	int thisbkg = 4;
	int l = 0;
	double weighted = 0;
	double squareweights = 0;
	int count = 0;
	// Loop over events
	while(bkg.read())
	{
		int channel = int(bkg.get("channel"));
		int index = int(bkg.get("index"));
		if((bkg.get("issamesign") )) continue;

		if(thischannel!=0) if(!(channel==thischannel)) continue;

		// Only use ZZ events for now

		if(filelist[1] == "sig.dat") if(bkg.get("index") != k) continue;
		l++;
		count ++;

		// Calculate changes due to shift in systematics
		string sample = sampleNames[int(bkg.get("index"))];
//		double weightratio = getsystematics.doline(data, channel, sample);
		double weightratio = 1.0;
		double weight = bkg.get("weight")*weightratio;

		weighted += weight;
		squareweights += weight*weight;
		// NN
		
		double y = bkg.get("mela");	

		
		if(!(count%10000)) cout << count  << "\t mass4l: " << bkg.get("mass4l") <<"\t channel: " << bkg.get("channel") << "\t sample: "  <<  bkg.get("sample") << "\t index: "  <<  bkg.get("index") << "\t y:" << y << endl;

//		if(y > mvacut) histos_bkg->Fill(y, bkg.get("mass4l"), weight);
//		cout << "channel: " << channel << "\t mva:" << y << "\t mass: " << bkg.get("mass4l") << endl;
		std::map<std::string, double> data;
		bkg.mget(data);
		if (applycuts(data, channel)) continue;

		histos_bkgs[index-1]->Fill( bkg.get("mass4l"), y, weight);

		if(!(index==4 || index==6))	histos_bkg4->Fill(bkg.get("mass4l"), y,  weight);

		if(!(index==thisbkg)) continue;
		histos_bkg->Fill( bkg.get("mass4l"), y, weight);
//		if(y > mvacut) histos_bkg->Fill(y, bkg.get("mass4l"), weight);
	}
		
	histos_bkg->Write();
	histos_bkg4->Write();

	bkg.close();

		
	cout <<"background\t" << histos_bkg->Integral() << "/" << weighted <<  endl;
	
//	exit(0);




	//---------------------------------------------------------------------------
	// Event Loop - SIGNAL
	//---------------------------------------------------------------------------

	float signal_weights[maxsig_plots];


	cout << "Making signal histograms" << endl;

	TString dataString2 = dirIN + filelist[2];
	Slurper sig(dataString2.Data());
	
	l = 0;
	count = 0;

	// Loop over events
	while(sig.read())
	{
		int channel = int(sig.get("channel"));
		if(thischannel!=0) if(!(channel==thischannel)) continue;
		if((sig.get("issamesign") )) continue;

		int index = sig.get("index");
		count ++;

		if(!(count%10000)) cout << count << endl;
		// Calculate changes due to shift in systematics
		
		string sample = sampleNames_sig;
//		double weightratio = getsystematics.doline(data, channel, sample);
		double weightratio = 1.0;
		double weight = sig.get("weight") * weightratio;
		
		if(index > maxsig_plots -1) continue;
		l++;
		signal_weights[index] += weight;

		double y;
		y = sig.get("mela");	

		if(!(count%10000)) cout << count  << "\t" << sigNames_mc[index] << "\t" << channel << "\t" << y << endl;
		
//		if(y > mvacut) histos_sig[index]->Fill(y, sig.get("mass4l"), weight);
//		if(y > mvacut) histos_sig_morph[index]->Fill(y, sig.get("mass4l"), weight);

		std::map<std::string, double> data;
		sig.mget(data);
		if (applycuts(data, channel)) continue;

		histos_sig_morph[index]->Fill( sig.get("mass4l"), y, weight);
		histos_sig[index]->Fill( sig.get("mass4l"),y, weight);		
		
		
	}


	cout <<"signal\t" << histos_sig[2]->Integral() << "/" << signal_weights[2] <<  endl;


	cout <<"signal\t" << histos_sig[2]->Integral() << "/" << signal_weights[2] <<  endl;

	// Calculate normalization to get the uncertainties right
	// Loop over files (sig, dat)
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		histos_sig[k]->Write();
		histos_sig_morph[k]->Write();

	}
	
	for(Int_t k=0; k < maxbkg; k++) 
	{
		histos_bkgs[k]->Write();
	}
	
	sig.close();
	f2.Close();

	end = clock();
	
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";
}
