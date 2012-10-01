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
	
	TFile f2("output/hist_ensemble_stdcuts.root");

	TCanvas* c6 = new TCanvas("morph_input","Sigmal Morph",1000,500);
	c6->cd();

	//---------------------------------------------------------------------------
	// Input Signal histograms
	//---------------------------------------------------------------------------	
	
	int N;
	TH2F* histos_sig[maxsig_plots];
	TH1F* histos_sig_1d[maxsig_plots];
	double weights[maxsig_plots];

	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_morph_" << sigNames_mc[k]; 
		TString name = _name.str();
		cout << name << endl;
		histos_sig[k]= (TH2F*) gDirectory->Get(name);
		cout << k << name << endl;
		N = histos_sig[k]->GetXaxis()->GetNbins();
		cout << N << endl;
		histos_sig_1d[k] = new TH1F(name, name, ybins_morph, ymin_morph, ymax_morph);
		// Conver 2d to 1d
		for(Int_t i=1; i < N; i++) 
		{
		histos_sig_1d[k]->SetBinContent( i, histos_sig[k]->ProjectionX()->GetBinContent(i) );
		}
		histos_sig_1d[k]->SetMaximum(4.0);
		histos_sig_1d[k]->SetTitle("Signal Morph");	
		histos_sig_1d[k]->GetXaxis()->SetTitle("m_{4l} (GeV)");
		histos_sig_1d[k]->GetYaxis()->SetTitle("Events");
		histos_sig_1d[k]->SetLineWidth(2.0);
		histos_sig_1d[k]->SetLineColor(2);	
		weights[k] = histos_sig_1d[k]->Integral();
	}	


	std::stringstream _name3;
	_name3 << "histos_sig_hist" << 2; 
	TString name3 = _name3.str();


	// output morphed histograms to a different file
	TFile f("output/hist_ensemble_morph.root", "RECREATE");
	f.cd();

	// print input histograms to image
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		if(k==0) histos_sig_1d[k]->Draw("HIST"); else histos_sig_1d[k]->Draw("HIST same");
		cout << histos_sig_1d[k]->Integral() << "\t" << histos_sig_1d[k]->GetXaxis()->GetNbins() << endl;
		c6->Update();
		histos_sig_1d[k]->Rebin(10);
	}

	c6->SaveAs(".gif");

	// Define variable
	RooRealVar x("x", "x", ymin_morph, ymax_morph);
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
		sigDataHist[k] = new RooDataHist(name2,name2, RooArgSet(x), histos_sig_1d[k]);
		
		std::stringstream _name4;
		_name4 << "g" << k; 
		TString name4 = _name4.str();

		const RooArgSet* dvars = sigDataHist[k]->get() ;
		cout << "vars: " << nset.getSize() << " dvasrs: " << dvars->getSize() << endl;
		g[k] = new RooHistPdf(name4, name4, RooArgSet(x), *sigDataHist[k]);

		pdfs.add(*g[k]);
		mref[k] = sigNames_mc[k];
	}


	RooMomentMorph *massmorph[nbins];

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
	for(Int_t i=ymin; i < ymax; i=i+2) 
	{
			count++;
		//	RooTrace::dump(cout,kTRUE);
		//	RooMomentMorph(const char* name, const char* title, RooAbsReal& _m, const RooArgList& varList, const RooArgList& pdfList, const RooArgList& mrefList, const RooMomentMorph::Setting& setting = NonLinearPosFractions)
	
			RooArgList var("variables"); 
			var.add(x);

			double thismass = double(i);
			RooRealVar alpha("alphaHist","alphaHist",ymin_morph,ymax_morph);
			RooMomentMorph::Setting setting = RooMomentMorph::NonLinear;
			alpha.setVal(thismass) ;


			massmorph[count] = new RooMomentMorph("moment morph", "moment morph", alpha, RooArgList(x), pdfs, mref, setting);


			cout << "Interpolated Mass: " << thismass << endl;	
//			cout << "Value: " << massmorph->getVal() << endl;		

			RooAbsPdf *morph_hist;
			morph_hist = massmorph[count]->sumPdf(&RooArgSet(x));

			TH1* hh = (TH1*) morph_hist->createHistogram("hh",x,Binning(N)) ;
			cout << "Morph Integral: " << itp.Eval(thismass) << endl;
			hh->Scale(itp.Eval(thismass)/hh->Integral());

			cout << "after" << endl;
			
			
			hh->SetLineColor(1);
			cout << "Integral: " << hh->Integral() << endl;
//			hh->Scale(1/hh->Integral());
			hh->Draw("same");
			c6->Update();
			hh->Rebin(10);



			std::stringstream _name1;
			_name1 << "histos_sig_" << thismass; 
			TString name1 = _name1.str();

			TH2F* output  = new TH2F(name1, name1, xbins, xmin, xmax, ybins, ymin, ymax);	


			int N1 = hh->GetNbinsX();
			cout << "N1" << N1 << endl;	
			
			int i_init = (ymin - ymin_morph)/2;
			
			for(Int_t i=i_init; i < ybins; i++)
			{
				output->SetBinContent( i-i_init, hh->GetBinContent(i), 1 );
			}

//			output->ProjectionY()->Draw("same");
			
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
		delete 		histos_sig_1d[k];
	}

}
