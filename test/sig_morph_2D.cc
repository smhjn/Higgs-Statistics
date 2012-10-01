//-----------------------------------------------------------------------------
// File:        posterior.C
// Description: 
// Bayesian fitting
//-----------------------------------------------------------------------------
// Created:     14-Nov-2011  Joe Bochenek
//-----------------------------------------------------------------------------


// g++ -I $ROOTSYS/include ratioploter.C `root-config --glibs ` `root-config --libs` `root-config --cflags` -Xlinker -lRooFit -zmuldefs -g -o ratioplot

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
#include <TTree.h>
#include <TLatex.h>
#include <TLegend.h>
#include <TGraph.h>
#include <TROOT.h>
#include <TChain.h>
#include <TFile.h>
#include <iostream>
#include "TPaveText.h"
#include <iomanip>
#include "RooNumIntConfig.h"
#include "RooBinning.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooIntegralMorph.h"
#include "RooMomentMorph.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFit.h"
#include "RooTrace.h"

#include "include/samples.h"
#include "Math/Interpolator.h"

using namespace RooFit;
using namespace std;

//-----------------------------------------------------------------------------
// MAIN MACRO
//-----------------------------------------------------------------------------

int main()
{
	
	TFile f2("output/hist_noiso.root");

	TCanvas* c6 = new TCanvas("morph_input","Sigmal Morph",1000,500);
	c6->cd();

	//---------------------------------------------------------------------------
	// Input Signal histograms
	//---------------------------------------------------------------------------	
	
	int Nx, Ny;
	TH1F* histos_sig_x[maxsig_plots];
	TH1F* histos_sig_y[maxsig_plots];
	TH2F* histos_sig[maxsig_plots];
	double weights[maxsig_plots];

	for(int k=0; k < maxsig; k++) 
	{
		sigNames[k] = ymin + k;
	}
	
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_morph_" << sigNames_mc[k]; 
		TString name = _name.str();

		cout << name << endl;
		histos_sig[k]= (TH2F*) gDirectory->Get(name);
		
		Nx = histos_sig[k]->GetXaxis()->GetNbins();
		cout << Nx << endl;
		Ny = histos_sig[k]->GetYaxis()->GetNbins();
		cout << Ny << endl;

		std::stringstream _name1;
		_name1 << "histos_sig_x_" << sigNames_mc[k]; 
		TString name1 = _name1.str();

		std::stringstream _name2;
		_name2 << "histos_sig_y_" << sigNames_mc[k]; 
		TString name2 = _name2.str();

		histos_sig_x[k] = new TH1F(name1, name1, ybins_morph, ymin_morph, ymax_morph);
		histos_sig_y[k] = new TH1F(name2, name2, 100, 0., 1.0);
		
		// Conver 2d to 1d
		for(Int_t i=1; i < Nx+1; i++) 
		{
			histos_sig_x[k]->SetBinContent( i, histos_sig[k]->ProjectionX()->GetBinContent(i) );
		}

		for(Int_t i=1; i < Ny+1; i++) 
		{
			histos_sig_y[k]->SetBinContent( i, histos_sig[k]->ProjectionY()->GetBinContent(i) );
		}
		
		histos_sig_x[k]->SetMaximum(4.0);
		histos_sig_x[k]->SetTitle("Signal Morph X");	
		histos_sig_x[k]->GetXaxis()->SetTitle("m_{4l} (GeV)");
		histos_sig_x[k]->GetYaxis()->SetTitle("Events");
		histos_sig_x[k]->SetLineWidth(2.0);
		histos_sig_x[k]->SetLineColor(2);	
		weights[k] = histos_sig_x[k]->Integral();
		cout << sigNames_mc[k] << "\t" << histos_sig[k]->Integral() << endl;	

		histos_sig_y[k]->SetMaximum(4.0);
		histos_sig_y[k]->SetTitle("Signal Morph Y");	
		histos_sig_y[k]->GetXaxis()->SetTitle("BNN(x)");
		histos_sig_y[k]->GetYaxis()->SetTitle("Events");
		histos_sig_y[k]->SetLineWidth(2.0);
		histos_sig_y[k]->SetLineColor(2);
	}
	

	std::stringstream _name3;
	_name3 << "histos_sig_hist" << 2; 
	TString name3 = _name3.str();


	// output morphed histograms to a different file
	TFile f("output/hist_ensemble_morph.root", "RECREATE");
	f.cd();

	
	// Define variable
	RooRealVar x("x", "x", ymin, ymax);
	RooArgSet nset("nset");
	nset.add(x);

	// Define pdf list
	RooArgList pdfs("pdfs"); 

	// give reference values (so the masses of the input histograms)
	TVectorD mref(maxsig_plots);

	// a vector of pdfs
	RooHistPdf *g[maxsig_plots];
	RooDataHist *sigDataHist[maxsig_plots];

	// convert the input histograms to roofit pdfs
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name2;
		_name2 << "histos_sig_hist" << k; 
		TString name2 = _name2.str();
		sigDataHist[k] = new RooDataHist(name2,name2, RooArgSet(x), histos_sig_x[k]);
		
		std::stringstream _name4;
		_name4 << "g" << k; 
		TString name4 = _name4.str();

		const RooArgSet* dvars = sigDataHist[k]->get() ;
		cout << "vars: " << nset.getSize() << " dvasrs: " << dvars->getSize() << endl;
		g[k] = new RooHistPdf(name4, name4, RooArgSet(x), *sigDataHist[k]);

		pdfs.add(*g[k]);
		mref[k] = sigNames_mc[k];
	}


	RooMomentMorph *massmorph[ybins];





	// Define variable
	RooRealVar y("y", "y", 0., 1.);
	RooArgSet nsety("nsety");
	nsety.add(y);

	// Define pdf list
	RooArgList pdfsy("pdfsy"); 

	// a vector of pdfs
	RooHistPdf *gy[maxsig_plots];
	RooDataHist *sigDataHisty[maxsig_plots];

	// convert the input histograms to roofit pdfs
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name2;
		_name2 << "histos_sig_histy" << k; 
		TString name2 = _name2.str();
		sigDataHisty[k] = new RooDataHist(name2,name2, RooArgSet(y), histos_sig_y[k]);

		std::stringstream _name4;
		_name4 << "g" << k; 
		TString name4 = _name4.str();

		gy[k] = new RooHistPdf(name4, name4, RooArgSet(y), *sigDataHisty[k]);

		pdfsy.add(*gy[k]);
	}
	RooMomentMorph *massmorphy[ybins];



	// print input histograms to image
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		if(k==0) histos_sig_x[k]->Draw(); else histos_sig_x[k]->Draw("same");
		cout << sigNames_mc[k] << "\t" << mref[k] << endl;
		c6->Update();
	}
	c6->SaveAs(".gif");



	// Interpolate the overall yield separately (spline int)
	cout << "Event Yield interpolation" << endl;
	cout << "Input Values: " << endl;
	vector <double> splinevalues;
	vector <double> _signals;
	for(Int_t k=0; k < maxsig_plots; k++) 
	{		
		_signals.push_back(sigNames_mc[k]);
		splinevalues.push_back(histos_sig[k]->Integral());	
		cout << sigNames_mc[k] << "\t" << histos_sig[k]->Integral() << endl;	
	}
	ROOT::Math::Interpolator itp( _signals, splinevalues, ROOT::Math::Interpolation::kCSPLINE);


	// loop through masses
	int count = 0;
	for(Double_t i=ymin; i < ymax; i=i+1) 
	{
			count++;
			RooArgList var("variables"); 
			var.add(x);

			double thismass = double(i);
			RooRealVar alpha("alphaHist","alphaHist",ymin_morph,ymax_morph);
			RooMomentMorph::Setting setting = RooMomentMorph::Linear;
			alpha.setVal(thismass);

			massmorph[count] = new RooMomentMorph("moment morph", "moment morph", alpha, RooArgList(x), pdfs, mref, setting);
			massmorphy[count] = new RooMomentMorph("moment morph y", "moment morph y", alpha, RooArgList(y), pdfsy, mref, setting);

			RooAbsPdf *morph_hist;
			RooAbsPdf *morph_histy;

			morph_hist = massmorph[count]->sumPdf(&RooArgSet(x));
			morph_histy = massmorphy[count]->sumPdf(&RooArgSet(y));

			TH1* hh = (TH1*) morph_hist->createHistogram("hh",x,Binning(ybins)) ;
			TH1* hhy = (TH1*) morph_histy->createHistogram("hh",y,Binning(xbins)) ;

			cout << "Binning N " << N << endl;

			hh->Scale(itp.Eval(thismass)/hh->Integral());

			hh->SetLineColor(1);
			hh->Draw("same");
			c6->Update();

			std::stringstream _name1;
			_name1 << "histos_sig_" << thismass; 
			TString name1 = _name1.str();

			TH2F* output  = new TH2F(name1, name1, ybins, ymin, ymax, xbins, xmin, xmax);	

			int N1 = hh->GetNbinsX();
						

			for(Int_t i=1; i < ybins +1; i++)
			{
			for(Int_t j=1; j < xbins + 1; j++)
			{
				output->SetBinContent(i, j, hh->GetBinContent(i)*hhy->GetBinContent(j) );
			}
			}


			for(Int_t i=1; i < ybins +1; i++)
			{
//				cout << histos_sig_y[1]->GetBinContent(i) << "\t";
				cout << thismass << "\t" << sigNames[i] << "\t" << hh->GetBinContent(i) << endl;
			}


			
			
			output->Write();
			delete output;
	}

	
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
	delete sigDataHist[k];
	delete g[k];
	}


//	RooTrace::dump(cout,kTRUE);

	
	f.Close();

	c6->SaveAs(".gif");
	
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		delete 		histos_sig_x[k];
	}

}
