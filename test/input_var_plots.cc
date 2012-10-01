//-----------------------------------------------------------------------------
// File:        input_var_plots.C
// Description: Data vs MC plots 
//-----------------------------------------------------------------------------
// Created:     20-June-2011  Joe Bochenek
//-----------------------------------------------------------------------------

#include <vector>
#include <fstream>
#include <sstream>
#include <algorithm>
#include <TH2.h>
#include <TStyle.h>
#include <TCanvas.h>
#include <TObject.h>
#include <TLorentzVector.h>
#include <THStack.h>
#include <TMath.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "Slurper.h"
#include "TPaveText.h"
#include <iomanip>
#include "samples.h"

using namespace std;

// apply cuts to the output events
int applycuts(std::map<std::string, double>& values, int channel){
	int notselected = 0;
//	cout << values["issamesign"] << endl;

	if(!(int(values["issamesign"]) ))  notselected++;

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

//-----------------------------------------------------------------------------
// MAIN MACRO
//-----------------------------------------------------------------------------
int main()
{
	int thischannel = 3;

	TString label = "_oppsign_pfxsel";
	TString imgDir = "output/";
	
	// Make the directory to store teh images.
	system("mkdir " + imgDir);
	system("mkdir " + imgDir + "inputplots/");
	
	TString description = "Input variables for H->WW->2l2nu after WW preselection<br><br><hr>";
	
	// Make an  html file to display plots.
	TString html_file = imgDir + "index.html";
	fstream html;
	html.open(html_file, ios_base::trunc | ios_base::out);
	if (!html || !html.good())
	{
	std::cout << "could not open file! " << html_file << endl;
	return -1;
	}
	html << "<html>";
	html << description << "\n\n";

	
	style();
	
	
	vector<string> var;
	var.push_back("nJets");
	var.push_back("sumET");
	var.push_back("nVertexes");
	var.push_back("pt2l");
	
	
	//---------------------------------------------------------------------------
	// Plot results
	//---------------------------------------------------------------------------
	
	
	vector<string> var1;
	var1.push_back("lept1_pt");
	var1.push_back("lept1_eta");
	var1.push_back("lept1_phi");
	var1.push_back("lept1_charge");
	var1.push_back("lept1_pfx");
	var1.push_back("lept1_sip");
	var1.push_back("lept1_mvaid");
	var1.push_back("lept2_pt");
	var1.push_back("lept2_eta");
	var1.push_back("lept2_phi");
	var1.push_back("lept2_charge");
	var1.push_back("lept2_pfx");
	var1.push_back("lept2_sip");
	var1.push_back("lept2_mvaid");
	var1.push_back("lept3_pt");
	var1.push_back("lept3_eta");
	var1.push_back("lept3_phi");
	var1.push_back("lept3_charge");
	var1.push_back("lept3_pfx");
	var1.push_back("lept3_sip");
	var1.push_back("lept3_mvaid");
	var1.push_back("lept4_pt");
	var1.push_back("lept4_eta");
	var1.push_back("lept4_phi");
	var1.push_back("lept4_charge");
	var1.push_back("lept4_pfx");
	var1.push_back("lept4_sip");
	var1.push_back("lept4_mvaid");
	var1.push_back("Z1mass");
	var1.push_back("Z2mass");
	var1.push_back("mass4l");

	
	
	string plotsdir = "plots1030/";
  
  
  
  
  
	// -------------------------------------------------------------------------------------------------------------------------------------
	// BOOK HISTOGRAMS
	// -------------------------------------------------------------------------------------------------------------------------------------
	
	// Data Histograms
	TH1F* counts_dat;  
	counts_dat = new TH1F("background", "", 1,0,1);
	counts_dat->SetBit(TH1::kCanRebin);
	counts_dat->SetStats(0);
	
	
	TH1F* dileptontype = new TH1F("dileptontype", "dileptontype", 4, 0, 4);   
	
	// Signal Histograms
	TH1F* counts_sig;  
	counts_sig = new TH1F("signal", "", 1,0,1);
	counts_sig->SetBit(TH1::kCanRebin);
	counts_sig->SetStats(0);
	
	
	// Total Histograms
	TH1F* counts_tot;  
	counts_tot = new TH1F("total bkg", "", 1,0,1);
	counts_tot->SetBit(TH1::kCanRebin);
	counts_tot->SetStats(0);
	
	// Signal Histograms After Cuts
	
	TH1F* histos_sig[N];
	for(Int_t i=0; i<N; i++)
	{
	TString title = var1[i] + "_sig";
	histos_sig[i] = new TH1F(title.Data(),title.Data(), binning[i], histo_min[i], histo_max[i]);
	histos_sig[i]->SetDirectory(0);
	histos_sig[i]->Sumw2();
	}
	
	// Background Histograms After Cuts
	TH1F* histos_bkg[N][maxbkg];
	
	
	TH1F* counts_bkg[maxbkg];  
	for(Int_t j=0; j<maxbkg; j++)
	{
	counts_bkg[j] = new TH1F(bkgNames[j], bkgNames[j], 1,0,1);
	counts_bkg[j]->SetBit(TH1::kCanRebin);
	counts_bkg[j]->SetStats(0);
	for(Int_t i=0; i<N; i++)
	{
	TString title = var1[i] + "_" + bkgNames[j] + "_dat2";
	histos_bkg[i][j]=new TH1F(title.Data(),title.Data(), binning[i], histo_min[i], histo_max[i]);
	histos_bkg[i][j]->SetDirectory(0);
	histos_bkg[i][j]->Sumw2();
	}
	}
	
	
	// Data Histograms After Cuts
	
	TH1F* histos_dat[N];
	for(Int_t i=0; i<N; i++)
	{
	TString title = var1[i] + "_dat1";
	histos_dat[i] = new TH1F(title.Data(),title.Data(), binning[i], histo_min[i], histo_max[i]);
	histos_dat[i]->SetDirectory(0);
	histos_dat[i]->Sumw2();
	}
	
	// Total Histograms After  
	
	TH1F* histos_total_bkg[N];
	for(Int_t i=0; i<N; i++)
	{
	TString title = var1[i] + "_total3";
	histos_total_bkg[i] = new TH1F(title.Data(),title.Data(), binning[i], histo_min[i], histo_max[i]);
	histos_total_bkg[i]->SetDirectory(0);
	histos_total_bkg[i]->Sumw2();
	}
	
	

	// Start processing the files and filling histograms

	//---------------------------------------------------------------------------
	// DATA
	//---------------------------------------------------------------------------
	int totdat = 0;
	int datsel = 0;

	TString datafile2 = dirIN +	"data.dat";
	//	cout << "\n Processing file " << datafile2 << endl;
	Slurper dat(datafile2.Data());
	while ( dat.read() )
	{
		int channel = int(dat.get("channel"));
//		if(!(channel==thischannel)) continue;
		totdat++;

		
		std::map<std::string, double> data;
		dat.mget(data);
//		if (applycuts(data, channel)) continue;
		datsel++;
		cout<< "\tmela:" << dat.get("mela")  << "\t" << dat.get("issamesign") << "\t" << dat.get("channel") << "\t" << dat.get("mass4l") << "\t" << dat.get("Z1mass") << "\t" << dat.get("Z2mass")<< "\t" << dat.get("lept4_pt")<< "\t" << dat.get("lept3_pt")<< "\t" << dat.get("lept2_pt")<< "\t" << dat.get("lept1_pt") << "\t" << dat.get("lept1_mvaid") << "\t" << dat.get("lept2_mvaid") << "\t" << dat.get("lept3_mvaid") << "\t" << dat.get("lept4_mvaid") << endl;

		double weight;
		weight = 1;
		for(Int_t i=0; i<N; i++)
		{
		histos_dat[i]->Fill(dat.get(var1[i]));
		}

	}
	dat.close();	
	cout << "Total Data: " << datsel<< "\t/\t" << totdat << endl;


//	exit(0);

	Double_t reweight[maxbkg] = {
	0.544732*0.925624,
	1.,
	1.,
	1.,
	1.
	};

	//---------------------------------------------------------------------------
	// BACKGROUND MC
	//---------------------------------------------------------------------------
	
	TString dataString = dirIN + "bkg.dat";

	cout << "\nProcessing file " << dataString << endl;

	Slurper bkg(dataString.Data());
	double counts_in[maxbkg] = { 0 };
	double counts_out[maxbkg] = { 0 };
	double counts_in_w[maxbkg] = { 0 };
	double counts_out_w[maxbkg] = { 0 };
	double totbkg = 0;
	int bkgsel = 0;
	int count = 0;

	// Loop over events
	while(bkg.read())
	{

// 			Counts and weights and crap like that
			int channel = int(bkg.get("channel"));
//			if(!(channel==thischannel)) continue;
//			cout << "HELLO!" <<  bkg.get("issamesign") << endl;
			
			int bkg_index = int(bkg.get("sample"));
			double weight = bkg.get("weight") * reweight[bkg_index - 1];
			count++;
			if(!(count%10000)) cout << count << "\tchannel: " << bkg_index << endl;
			counts_in[bkg_index-1] ++;
			counts_in_w[bkg_index-1] += weight;
			totbkg+= weight;

			std::map<std::string, double> data;
			bkg.mget(data);
//			if (applycuts(data, channel)) continue;
			counts_out_w[bkg_index-1] += weight;
			
			
			bkgsel++;			
			// Fill plots
			for(Int_t i=0; i<N; i++)
			{		
//				 cout << i << "/" << N << "\t" << bkg_index << "/" << maxbkg << "\t" << bkg.get(var1[i]) << "\t" << var1[i] << endl;
				 histos_bkg[i][bkg_index-1]->Fill(bkg.get(var1[i]), weight); 
				 histos_total_bkg[i]->Fill(bkg.get(var1[i]), weight);
			}
			

	}
	
	
	for(int j = 0; j < maxbkg; j++){
		cout << j << "\t" << bkgNames[j]  << "\t " << counts_out[j] << " / " << counts_in[j]   << "\t " << counts_out_w[j] << " / " << counts_in_w[j] << endl;
		bkg.close();
	}

	cout << endl;
	cout << "Background: " << totbkg << endl;
	cout << "Background: " << bkgsel << " / " << count << endl;

	cout << endl;

	// Choose only one signal file to print
	int onesignal = 2; 

	
	//---------------------------------------------------------------------------
	// SIGNAL
	//---------------------------------------------------------------------------

	TString sigFile = dirIN + "sig_bnn_train_uw.dat";
	cout << "\n Processing file " << sigFile << endl;
	Slurper sig(sigFile.Data());

	int totalsig = 0;
	double selsig = 0;
	 count = 0;

	for(int j = 0; j < maxsig; j++){
	while ( sig.read() )
	{
		count++;
		if(!(count%100000)) cout << count << endl;

		int sig_index = int(sig.get("index"));
		if(!(sig_index == onesignal)) continue;
		
		int channel = int(sig.get("channel"));
//		if(!(channel==thischannel)) continue;
		
		totalsig++;

		std::map<std::string, double> data;
		sig.mget(data);
//		if (applycuts(data, channel)) continue;


		double weight;
		weight = sig.get("weight");
		selsig+=weight;
		
		for(Int_t i=0; i<N; i++)
		{
			 histos_sig[i]->Fill(sig.get(var1[i]), weight);
		}

	}
	}
	
	sig.close();
	cout << "Signal (" << sigNames[onesignal] << ")\t " << totalsig << "\tweighted: " << selsig << endl;
	cout << endl;
	


	// Print table of yields
	ofstream log;
	log.open("output/preselection_yields.txt", ios_base::app);
	
	log << "---------------------------------------------------------------------------" << endl;
	
	for(int j = 0; j < maxbkg; j++){
	log << bkgNames[j] << "\t\t";
	}
	
	log << "Higgs\t\t";
	log << "Totbkg\t\t";
	log << "Data" << endl;



  
// ----------------------------------------------------------------------------------------
//  INPUT VARIABLE PLOTS
// ----------------------------------------------------------------------------------------

	for(int i=0; i<N; i++)
	{	
		stringstream title_string;
		title_string << varnames[i];
		TString title = title_string.str();
		
		cout << " ----------------- " << title << "---------------" << endl;
		
		double total1=histos_dat[i]->Integral();
		double total2=histos_bkg[i][0]->Integral();
		double total3=histos_total_bkg[i]->Integral();
		double weightratio2 = total2/total1;
		double weightratio1 = total1/total3;
		cout << "Total Data: " << total1 << "Total Background: " << total3 << "\tRatio MC/Data: " << weightratio1 << endl;		
		stringstream imgdir_string;
		imgdir_string << imgDir << "/inputplots/";    
		TString imgdir = imgdir_string.str();
		TString nomeOUT = imgdir + var1[i] + label + ".png";
		TString imageName = "/inputplots/" + var1[i] + label + ".png";

     	html << "<a href=\"" + imageName + "\"><img img width=\"300\" height=\"360\" src=\"" + imageName + "\">\n";
		if(!((i+1)%3)) html << "<br>\n";


		// ------- Kolmogorov and Chi Test ------------
		Double_t koltest2 = histos_dat[i]->KolmogorovTest(histos_total_bkg[i], "");
		Double_t chi2test2 = histos_dat[i]->Chi2Test(histos_total_bkg[i], "UW");
		
		stringstream kol2, kol3, chi2, chi3;
		kol2 << setprecision(4)<< fixed;
		chi2 << setprecision(4)<< fixed;
		kol3 << setprecision(4) << fixed;
		chi3 << setprecision(4) << fixed;
		
		kol2 << "KS-prob = " << koltest2;
		chi2 << "#chi^{2}-prob: = " << chi2test2;
		
		TString kol_2 = kol2.str();
		TString chi_2 = chi2.str();
		TString kol_3 = kol3.str();
		TString chi_3 = chi3.str();
		
		//	cout <<  "Var " << var1[i] << "\t K Test: " << koltest << ", dist: " << koltestdist << "\t Chi: " << chi2test << "; value: " << chi2test2 << "\t Weighting: " << (double) totals/totalb << endl;
		
		TPaveText *pt3 = new TPaveText(0.45,0.85,0.8,0.93, "bordersize=0 NDC"); // NDC sets coords
		pt3->SetFillStyle(4000); // text is black on white
		pt3->SetTextSize(0.03); 
		pt3->SetTextAlign(12);
		pt3->SetShadowColor(0);
		pt3->AddText(kol_3);
		pt3->AddText(chi_3);
		pt3->SetBorderSize(0);
		
		
		TPaveText *pt2 = new TPaveText(0.15,0.81,0.5,0.89, "bordersize=0 NDC"); // NDC sets coords
		pt2->SetFillStyle(4000); // text is black on white
		pt2->SetTextSize(0.03); 
		pt2->SetTextAlign(12);
		pt2->SetShadowColor(0);
		pt2->AddText(kol_2);
		pt2->AddText(chi_2);
		pt2->SetBorderSize(0);
		
		
		// Make Canvas
		TCanvas *c1 = new TCanvas(title, title, 800, 960);
		c1->cd();
		
		TPad* pad1 = new TPad(title, title, 0.0, 0.3, 1.0, 1.0);
		if(logscale[i]) pad1->SetLogy();
		pad1->Draw();
		c1->Update();
		
		TPad* pad2 = new TPad(title, title, 0.0,0.0, 1.0, 0.3);
		pad2->SetGridy();
		pad2->Draw();
		c1->Update();
		
		THStack *htotal1=new THStack(title,title);
		
		TPaveText *ttitle = new TPaveText(0.0,0.93,0.45,1.0, "bordersize=5 NDC"); // NDC sets coords
		ttitle->SetFillColor(0); // text is black on white
		ttitle->SetTextSize(0.05); 
		ttitle->SetTextAlign(12);
		ttitle->SetShadowColor(0);
		ttitle->AddText(title);
		
		TLegend *leg1 = new TLegend(0.7,0.8,1.0, 1.0);
		leg1->SetFillColor(4000);
		leg1->SetShadowColor(0);
		
		histos_dat[i]->SetLineColor(kRed);
		histos_dat[i]->SetLineWidth(2);
		
		histos_dat[i]->SetMarkerStyle(21);
		histos_dat[i]->SetMarkerSize(1);
		histos_dat[i]->SetMarkerColor(1);
		
		leg1->SetNColumns(2);
		leg1->AddEntry(histos_dat[i],"Data (4.71 fb^{-1})","LEP");

		
		histos_dat[i]->SetLineColor(1);
		histos_dat[i]->SetLineWidth(2.0);
		
		histos_total_bkg[i]->SetFillColor(18);
		histos_total_bkg[i]->SetLineColor(18);
		
		
		for(Int_t j=0; j<maxbkg; j++)
		{
		histos_bkg[i][j]->GetXaxis()->SetTitle(title);
		histos_bkg[i][j]->GetYaxis()->SetTitle("Events");
		histos_bkg[i][j]->SetFillColor(colors[j]);
		histos_bkg[i][j]->SetLineColor(1);
		histos_bkg[i][j]->SetLineWidth(0.5);
		htotal1->Add(histos_bkg[i][j], "");
		leg1->AddEntry(histos_bkg[i][j],bkgNames[j],"F");
		}
		

//		leg1->AddEntry(histos_sig[i],signalTag[onesignal],"L");
		histos_sig[i]->SetLineColor(2);
		histos_sig[i]->SetLineWidth(3.0);
		
		
		// ----------------- Upper Plot (Variables) -----------------
		
		pad1->cd();
		histos_dat[i]->GetXaxis()->SetTitle(title);
		leg1->Draw();
//		htotal1->SetMaximum(7);
		htotal1->Draw("HIST");
		htotal1->GetXaxis()->SetTitle(title);
		htotal1->GetYaxis()->SetTitle("Events");

		htotal1->SetTitle(title);

//		htotal1->GetHistogram()->GetYaxis()->SetTitle("Events");

		gPad->Modified();

		gPad->RedrawAxis();

		histos_dat[i]->Draw("LEP same");
//		histos_sig[i]->Draw("HIST same");

		gPad->RedrawAxis();
	
		
		leg1->Draw();
		pad1->Draw();
		pt2->Draw();
//		ttitle->Draw();
		c1->Modified();
		c1->Update();
		
		
		// ----------------- Lower Plot (Ratio Plot) -----------------
		
		TH1F* roofit_plot_pt2l_2 = new TH1F("Pt2l1oPt2l2_2","Pt2l1oPt2l2_2", binning[i], histo_min[i], histo_max[i]);
		for( int k = 1; k < binning[i] + 1; k++){
			double dat1count = histos_dat[i]->GetBinContent(k);
			double dat2count = histos_total_bkg[i]->GetBinContent(k);
		if(dat2count > 0 && dat1count > 0) {
		double ratio = (dat1count)/(dat2count);
		roofit_plot_pt2l_2->SetBinContent(k, ratio);
		double error = sqrt(dat1count/(dat2count*dat2count) + (dat1count*dat1count)/(dat2count*dat2count*dat2count));
		if(error < 2) roofit_plot_pt2l_2->SetBinError(k, error);
		}
		}	
		pad2->cd();
		//	roofit_plot_pt2l_2->GetYaxis()->SetNdivisions(6);
		roofit_plot_pt2l_2->SetMaximum(3.5);
		roofit_plot_pt2l_2->SetMinimum(-1.5);
		roofit_plot_pt2l_2->Draw(); 

		c1->Modified();
		c1->Update();
		
		c1->cd();
		c1->SaveAs(nomeOUT);

		
		/*
		// ------- QQ ------------
		
		int nq = 20;
		Double_t xq[nq];  // position where to compute the quantiles in [0,1]
		Double_t yq1[nq];  // array to contain the quantiles
		Double_t yq2[nq];  // array to contain the quantiles
		
		c1->cd(2);
		for (Int_t i=0;i<nq;i++) xq[i] = Float_t(i+1)/nq;
		histos_dat[i]->GetQuantiles(nq, yq1, xq);
		histos_bkg[i]->GetQuantiles(nq, yq2, xq);
		
		TGraph *gr = new TGraph(nq,yq1,yq2);
		gr->GetXaxis()->SetTitle("MC quantiles");
		gr->GetYaxis()->SetTitle("Data quantiles");
		gr->Draw("alp");
		c1->Update();
		*/	
		}  

		double m2lsum_other = 0;
		for(int i=2; i<maxbkg; i++)
		{
			m2lsum_other+=histos_bkg[8][i]->Integral(0,16);		
		}
		double m2lsum_ZZ = histos_bkg[8][0]->Integral(0,16) + histos_bkg[8][1]->Integral(0,16);		
		double m2lsum_data = histos_dat[8]->Integral(0,16);
		double factorZZ_R =  (m2lsum_data - m2lsum_other) / m2lsum_ZZ;
		cout << "ZZ scale factor: " << factorZZ_R << endl;



// ----------------------------------------------------------------------------------------
//  Plots Showing Counts
// ----------------------------------------------------------------------------------------

	// Make Canvas
	TCanvas *c2 = new TCanvas("counts", "counts", 800, 960);
	c2->cd();
	c2->SetLogy();
	TLegend *legcounts = new TLegend(0.75,0.80,0.95, 1.0);
	legcounts->SetFillColor(kWhite);
	legcounts->SetShadowColor(0);

    counts_sig->SetLineWidth(1.0);
    counts_sig->SetLineColor(1);
    counts_sig->SetFillColor(4000);
    counts_dat->SetLineColor(1);
    counts_dat->SetFillColor(4000);
    counts_dat->SetLineWidth(2.0);


   	THStack *hcounts=new THStack("counts","counts");

    for(Int_t j=0; j<maxbkg; j++)
    {
		counts_bkg[j]->SetLineWidth(1.0);
		counts_bkg[j]->SetLineColor(j+2);
		counts_bkg[j]->SetFillColor(j+2);
		hcounts->Add(counts_bkg[j], "E");
		legcounts->AddEntry(counts_bkg[j],bkgNames[j],"F");
	}

  	legcounts->AddEntry(counts_sig,"signal","L");
  	legcounts->AddEntry(counts_dat,"data","L");

    hcounts->Add(counts_tot, ""); 	

    counts_tot->SetFillColor(18);
	counts_tot->SetLineColor(18);

    counts_dat->Draw("LEP");
	hcounts->Draw("HIST same");
    counts_sig->Draw("L same");


	stringstream imgdir_string;
    imgdir_string << imgDir << "/inputplots/";    
    TString imgdir = imgdir_string.str();
	TString nomeOUT = imgdir + "counts.png";

	legcounts->Draw();
 	c2->SaveAs(nomeOUT);


	TCanvas *c3 = new TCanvas("counts", "counts", 400, 400);
	c3->cd();
	dileptontype->GetXaxis()->SetBinLabel(1, "2#mu");
	dileptontype->GetXaxis()->SetBinLabel(2, "e#mu");
	dileptontype->GetXaxis()->SetBinLabel(3, "2e");
	dileptontype->GetXaxis()->SetBinLabel(4, "");
  	dileptontype->Draw();

	TString nomeOUT2 = imgdir + "dileptontype.png";
   	c3->SaveAs(nomeOUT2);
  
  return 0;
}
