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
#include "Math/Interpolator.h"

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
#include "Slurper.h"
#include "TPaveText.h"
#include <iomanip>
#include "RooNumIntConfig.h"
#include "RooBinning.h"
#include "RooGlobalFunc.h"
#include "RooRealVar.h"
#include "RooDataHist.h"
#include "RooIntegralMorph.h"
#include "RooHistPdf.h"
#include "RooPlot.h"
#include "RooFit.h"

using namespace RooFit;
using namespace std;

#include "samples.h"
#include "PoissonGammaFit.h"
#include "PoissonGammaIntegral.h"


//-----------------------------------------------------------------------------
// MAIN MACRO
//-----------------------------------------------------------------------------

int main()
{
	//TString html_name = "/afs/cern.ch/user/j/jpb/www/studies/HWW_analysis/plots/index.html";  
	//ofstream event_list;
	//event_list.open(html_name, ios_base::trunc);
	
	TFile f2("output/hist_ensemble.root");
	style();


	TCanvas* c = new TCanvas("morph_output","Sigmal Morph",1000,500);
	c->cd();

	TCanvas* c6 = new TCanvas("morph_input","Sigmal Morph",1000,500);
	c6->cd();

	//---------------------------------------------------------------------------
	// Input Signal histograms
	//---------------------------------------------------------------------------	
	
	int N;
	TH2F* histos_sig[3];
	TH1F* histos_sig_1d[3];
	double weights[3];

	for(Int_t k=0; k < 3; k++) 
	{
		std::stringstream _name;
		_name << "histos_sig_" << sigNames_mc[k]; 
		TString name = _name.str();
		cout << name << endl;
		histos_sig[k]= (TH2F*) gDirectory->Get(name);
		cout << k << name << endl;
		N = histos_sig[k]->GetYaxis()->GetNbins();
		cout << N << endl;
		histos_sig_1d[k] = new TH1F(name, name, ybins_morph, ymin_morph, ymax_morph);
		// Conver 2d to 1d
		for(Int_t i=1; i < N; i++) 
		{
			histos_sig_1d[k]->SetBinContent( i, histos_sig[k]->GetBinContent(1,i) );
		}
		histos_sig_1d[k]->SetMaximum(4.0);
		histos_sig_1d[k]->SetTitle("Signal Morph");	
		histos_sig_1d[k]->GetXaxis()->SetTitle("m_{4l} (GeV)");
		histos_sig_1d[k]->GetYaxis()->SetTitle("Events");
		histos_sig_1d[k]->SetLineWidth(2.0);
		histos_sig_1d[k]->SetLineColor(2);	
		weights[k] = histos_sig_1d[k]->Integral();
		cout << sigNames_mc[k] << "\t" << histos_sig[k]->Integral() << endl;	
	}	

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

	cout << "\nA few test Values: " << endl;
	cout << "125\t" <<  itp.Eval(125) << endl;
	cout << "145\t" <<  itp.Eval(145) << endl;


	TCanvas* c2 = new TCanvas("xsec_morph","Sigmal Yield Morph",1000,500);
	c2->cd();
	TGraph *gr = new TGraph();
	TGraph *gr2 = new TGraph();

	for(Int_t k=0; k < maxsig_plots; k++) 
	{		
	gr2->SetPoint( k, sigNames_mc[k], histos_sig[k]->Integral());
	}

	for(Int_t m=0; m < 1000; m++)  {
	double itpmass = sigNames_mc[0] + m*(sigNames_mc[maxsig_plots-1] - sigNames_mc[0])/1000;
	double itpval = itp.Eval(itpmass);
//	cout << sigNames_mc[0] << "\t" << itpmass << "\t" << itpval << endl;
	gr->SetPoint( m, itpmass, itpval);
	}

	gr->SetTitle("Signal Yield Morph");
	gr->GetXaxis()->SetTitle("m_{4l}");
    gr->GetXaxis()->SetTitleFont(42);
	gr->GetYaxis()->SetTitle("yield");
	gr->SetLineColor(2);
	gr->SetMarkerStyle(21);
	gr->SetMarkerSize(0.5);

	gr2->SetMarkerStyle(21);
	gr2->SetMarkerColor(2);
	gr2->SetMarkerSize(1.0);
	gr2->SetTitle("Signal Yield Morph");
	gr2->GetXaxis()->SetTitle("m_{4l}");
    gr2->GetXaxis()->SetTitleFont(42);
	gr2->GetYaxis()->SetTitle("Events");

	gr->Draw("ACP");
	gr2->Draw("ACP same");

	c2->Update();
	c2->Modified();
	c2->SaveAs("output/yield_morph.pdf");
	
	TCanvas* c4 = new TCanvas("morphtest","Morph and Endpoint PDFs",1000,500);
	c4->cd();

//	return 0;

	TFile f("output/hist_ensemble_morph.root", "UPDATE");


	RooIntegralMorph *massmorph = new RooIntegralMorph[maxsig_plots];
	for(Int_t k=0; k < maxsig_plots-1; k++) 
	{
	
		RooBinning fullRange( N );
		
		RooRealVar x("x", "x", ymin, ymax);
		x.setBinning(fullRange, "fullRange");
		
		//RooFormulaVar alpha("alpha","(@1-@2+@3)/(@1-@0)",   
		RooRealVar alpha("alphaHist","alphaHist",0,1);

		std::stringstream _name2;
		_name2 << "histos_sig_hist" << k; 
		TString name2 = _name2.str();

		std::stringstream _name3;
		_name3 << "histos_sig_hist" << (k+1); 
		TString name3 = _name3.str();

		
		RooDataHist *sigDataHist1 = new RooDataHist(name2,"",  RooArgSet(x), histos_sig_1d[k+1]);
		RooDataHist *sigDataHist2 = new RooDataHist(name3,"",  RooArgSet(x), histos_sig_1d[k]);
		
		// Lower end point shape: a histogram shaped p.d.f.
		RooHistPdf *g1 = new RooHistPdf("g1", "", RooArgSet(x), *sigDataHist1);
		// Upper end point shape: also a histogram shaped p.d.f.
		RooHistPdf *g2 = new RooHistPdf("g2", "", RooArgSet(x), *sigDataHist2);


		std::stringstream _name;
		_name << "sigModel" << k; 
		TString name = _name.str();

		alpha.setBins(N,"cache") ;
		massmorph[k] = RooIntegralMorph(name,"Signal Shape before resolution",*g1,*g2,x,alpha);
		massmorph[k].setCacheAlpha(kTRUE) ;
		
		TCanvas* c4 = new TCanvas("morphtest","Morph and Endpoint PDFs",1000,500);
		c4->cd();
		RooPlot* frameRatio= x.frame(5000);

		for(Int_t i=1; i < 11; i++) 
		{
			
			double thisalpha = 0.1*i;
	
			alpha.setVal(thisalpha);
	//		massmorph.plotOn(frameRatio,  LineColor(kBlack));
			//alpha.setVal(0.8);
			//massmorph.plotOn(frameRatio,  LineColor(kOrange));
			c4->cd();
			frameRatio->Draw();
			//TFile f1("morph.root", "RECREATE");
			//set the precision for the integral of the model here: 10^-6 should be fine, default is 10^-8.
			RooNumIntConfig* cfg = RooAbsReal::defaultIntegratorConfig();
			cfg->setEpsAbs(1E-5);
			cfg->setEpsRel(1E-5);
			massmorph[k].setIntegratorConfig(*cfg);
	
			TH1* hh = massmorph[k].createHistogram("hh",x,Binning(N)) ;
			hh->Draw();
//			hh->Write();
			
			double thismass = sigNames_mc[k] + thisalpha*(sigNames_mc[k+1] - sigNames_mc[k]);
			
			std::stringstream _name1;
			_name1 << "histos_sig_" << thismass; 
			TString name1 = _name1.str();


			cout << "Histogram before: " << thismass << "\t" << hh->Integral() << endl;
			cout << "Morph Integral: " << itp.Eval(thismass) << endl;
			hh->Scale(itp.Eval(thismass)/hh->Integral());
			cout << "Histogram after: " << 			hh->Integral() << endl;

			c->cd();
			hh->SetLineWidth(1.0);
			hh->SetLineColor(1);		
			if(!(i%2)) hh->Draw("c same");
			
			TH2F* output  = new TH2F(name1, name1, xbins, xmin, xmax, N, ymin, ymax);	
			int N1 = hh->GetNbinsX();
			cout << N1 << endl;
			for(Int_t i=1; i < N; i++) 
			{
				output->SetBinContent( 1, i, hh->GetBinContent(i) );
			}

			output->Write();
			cout << hh->GetXaxis()->GetNbins() << endl;
			cout << hh->GetYaxis()->GetNbins() << endl;
	
		}
	
		delete sigDataHist1;
		delete sigDataHist2;
		delete g1;
		delete g2;
	
//		histos_sig[k]->Write();
	
		c4->Update();
		c4->Modified();
		c4->SaveAs("output/morphtest.pdf");
	}
	
	c->SaveAs("output/sig_morphs.pdf");

	histos_sig[0]->Write();

	
	f.Close();
	
	
	return 0;
}