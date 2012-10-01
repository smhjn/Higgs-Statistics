//-----------------------------------------------------------------------------
// File:        common_likelihood.h
// Description: Initiate likelihood class using data and MC
//
//-----------------------------------------------------------------------------
// Created:     11-04-2012  Joe Bochenek

//-----------------------------------------------------------------------------
#ifndef COMMON_H
#define COMMON_H

#include "TPaveText.h"
#include <algorithm>
#include <cmath>
#include <fstream>
#include <iomanip>
#include <iostream>
#include <sstream>
#include <string>
#include <TCanvas.h>
#include <TChain.h>
#include <TFile.h>
#include <TGaxis.h>
#include <TGraph.h>
#include <TGraphErrors.h>
#include <TH2.h>
#include <THStack.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TLorentzVector.h>
#include <TMath.h>
#include <TObject.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TStyle.h>
#include <TTree.h>
#include <vector>

using namespace std;
#include "samples.h"

TString scripttag = "input1D";	

//---------------------------------------------------------------------------
// RETRIEVE HISTOGRAMS
//---------------------------------------------------------------------------
TH2F* histos_sig[maxsig];

TFile file("output/hist_ensemble_noiso.root");
file.cd();

// ZZ background
TH2F *histos_bkg=(TH2F*) gDirectory->Get("bkg2d");
histos_bkg->SetDirectory(0);
histos_bkg->Sumw2();
histos_bkg->GetXaxis()->SetTitle("BNN(x)");
histos_bkg->GetYaxis()->SetTitle("m_{4l}");
double signal_weights[maxsig] = {0};

// Data driven reducible background
TH2F *histos_zjets=(TH2F*) gDirectory->Get("zjets2d");
histos_zjets->SetDirectory(0);
histos_zjets->Sumw2();
histos_zjets->GetXaxis()->SetTitle("BNN(x)");
histos_zjets->GetYaxis()->SetTitle("m_{4l}");

// Data
TH2F *histos_dat=(TH2F*) gDirectory->Get("dat2d");
histos_dat->SetDirectory(0);	
histos_dat->Sumw2();
histos_dat->GetXaxis()->SetTitle("BNN(x)");
histos_dat->GetYaxis()->SetTitle("m_{4l}");

TFile f5("output/hist_ensemble_morph.root");
f5.cd();

//---------------------------------------------------------------------------
// FILL SIGNAL VECTOR
//---------------------------------------------------------------------------

for(int k=0; k < maxsig; k++) 
{
	sigNames[k] = ymin + k;
}


for(Int_t k=0; k < maxsig; k++) 
{
	std::stringstream _name;
//		sigNames[k] = double(2*k + 116);
	_name << "histos_sig_" << sigNames[k]; 
	TString name = _name.str();
	cout << k << "\t" << sigNames[k] << "\t" << name << endl;
	histos_sig[k]= (TH2F*)gDirectory->Get(name);
	histos_sig[k]->Rebin(25);
	histos_sig[k]->Scale(lumi_scale);
	signal_weights[k] = histos_sig[k]->Integral();
	cout << sigNames[k] << ": " << signal_weights[k] << endl;
	cout << "signal[" << k << "]: " << sigNames[k] << "\t" << signal_weights[k] << endl;
	histos_sig[k]->Sumw2();
}





//---------------------------------------------------------------------------
// FILL DATA VECTOR
//---------------------------------------------------------------------------

vector<double> mva_src_d;
vector<double> mva_error_d;			

// BACKGROUND
for(int i=1; i < xbins+1; i++)
{
for(int j=1; j < ybins+1; j++)
{
	double bincontent = histos_dat->GetBinContent(i, j);
	mva_src_d.push_back(bincontent);
	//				cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " exp error: " << exp_error << " exp error: " << tot_error  << endl;
}
}


//---------------------------------------------------------------------------
// FILL BACKGROUND VECTORS
//---------------------------------------------------------------------------

histos_bkg->Scale(lumi_scale);			
// Calculate errors in bins for the backgrounds
vector<double> mva_src_b;
vector<double> mva_error_b;			
// BACKGROUND
// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
for(int i=1; i < xbins+1; i++)
{
for(int j=1; j < ybins+1; j++)
{
	double bincontent = histos_bkg->GetBinContent(i, j);
	double binerror = histos_bkg->GetBinError(i, j);
	double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
	if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
	mva_src_b.push_back(bincontent);
	mva_error_b.push_back(error);
//		cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " eff error: " << bincontent*error << " rel error: " << binerror/bincontent << endl;
	//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
}
}


histos_zjets->Scale(0.0484985);
histos_zjets->Scale(lumi_scale);
histos_zjets->Smooth();
histos_zjets->Smooth();
histos_zjets->Smooth();
histos_zjets->Smooth();
		
// Calculate errors in bins for the backgrounds
vector<double> mva_src_zjets;
vector<double> mva_error_zjets;			
// BACKGROUND
for(int i=1; i < xbins+1; i++)
{
for(int j=1; j < ybins+1; j++)
{
	double bincontent = histos_zjets->GetBinContent(i, j);
	double binerror = histos_zjets->GetBinError(i, j);
	double error = 0;
	if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
	mva_src_zjets.push_back(bincontent);
	mva_error_zjets.push_back(error);
}
}




#endif
