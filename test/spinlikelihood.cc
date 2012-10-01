//-----------------------------------------------------------------------------
// File:        spinlikelihood.C
// Description: Likelihood for scalar vs pseudoscalar Higgs with background
//				
//-----------------------------------------------------------------------------
// Created:     17-09-2012  Joe Bochenek

//-----------------------------------------------------------------------------

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
#include <gsl/gsl_sf_gamma.h>

#include <TGraph2D.h>

using namespace std;

//#include "samples.h"
#include "PoissonGammaFit2.h"
#include "PoissonGammaIntegral.h"

int main(int argc, char* argv[]){

	TFile file("scripts/mvaplots/models.root");
	file.cd();

	TH2F *histos_zz=(TH2F*) gDirectory->Get("htempzz");
	histos_zz->SetDirectory(0);
	histos_zz->Sumw2();
	histos_zz->GetXaxis()->SetTitle("BNN(x)");
	histos_zz->GetYaxis()->SetTitle("m_{4l}");

	// Data driven reducible background
	TH2F *histos_zx=(TH2F*) gDirectory->Get("htempzx");
	histos_zx->SetDirectory(0);
	histos_zx->Sumw2();
	histos_zx->GetXaxis()->SetTitle("BNN(x)");
	histos_zx->GetYaxis()->SetTitle("m_{4l}");

	// 	Data
	//	TH2F *histos_dat=(TH2F*) gDirectory->Get("htempdat");
	//	histos_dat->SetDirectory(0);	
	//	histos_dat->Sumw2();
	//	histos_dat->GetXaxis()->SetTitle("BNN(x)");
	//	histos_dat->GetYaxis()->SetTitle("m_{4l}");
	
	// ZZ
	TH2F *histos_sc=(TH2F*) gDirectory->Get("htemps");
	histos_sc->SetDirectory(0);	
	histos_sc->Sumw2();
	histos_sc->GetXaxis()->SetTitle("BNN(x)");
	histos_sc->GetYaxis()->SetTitle("m_{4l}");
	
	// ZX
	TH2F *histos_ps=(TH2F*) gDirectory->Get("htempps");
	histos_ps->SetDirectory(0);	
	histos_ps->Sumw2();
	histos_ps->GetXaxis()->SetTitle("BNN(x)");
	histos_ps->GetYaxis()->SetTitle("m_{4l}");
	
	

	double lumi_scale = 5.;
	double ybins = 20;
	double xbins = 20;
	
	
	//---------------------------------------------------------------------------
	// FILL ZZ VECTORS
	//---------------------------------------------------------------------------

	histos_zz->Scale(lumi_scale*0.1497*19.0042897346/(0.9152*histos_zz->Integral()));
	histos_zz->Smooth();

	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_zz;
	vector<double> mva_error_zz;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_zz->GetBinContent(i, j);
		double binerror = histos_zz->GetBinError(i, j);
		double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_zz.push_back(bincontent);
		mva_error_zz.push_back(error);
	}
	}	


	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
//		cout << i << ", " << j << ":\t" << histos_zz->GetBinContent(i,j) << "\t" << histos_zx->GetBinContent(i,j) << "\t" << histos_ps->GetBinContent(i,j) << "\t" << histos_sc->GetBinContent(i,j) << "\t" << endl;
	}
	}


	//---------------------------------------------------------------------------
	// FILL ZX VECTORS
	//---------------------------------------------------------------------------

	histos_zx->Scale(0.2*lumi_scale);			
	histos_zx->Smooth();
	histos_zx->Smooth();
	histos_zx->Smooth();
	histos_zx->Smooth();

	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_zx;
	vector<double> mva_error_zx;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_zx->GetBinContent(i, j);
		double binerror = histos_zx->GetBinError(i, j);
		double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_zx.push_back(bincontent);
		mva_error_zx.push_back(error);
//		cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " eff error: " << bincontent*error << " rel error: " << binerror/bincontent << endl;
		//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
	}
	}	

	//---------------------------------------------------------------------------
	// FILL PS VECTORS
	//---------------------------------------------------------------------------

	double psscale = (lumi_scale* 1.1509305566+0.5347226809+ 1.6270259788)/ histos_ps->Integral();
	histos_ps->Scale(psscale);			
	histos_ps->Smooth();


	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_ps;
	vector<double> mva_error_ps;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_ps->GetBinContent(i, j);
		double binerror = histos_ps->GetBinError(i, j);
		double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_ps.push_back(bincontent);
		mva_error_ps.push_back(error);
//		cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " eff error: " << bincontent*error << " rel error: " << binerror/bincontent << endl;
		//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
	}
	}	
	


	//---------------------------------------------------------------------------
	// FILL Scalar VECTORS
	//---------------------------------------------------------------------------

	double scscale = (lumi_scale* 1.1509305566+0.5347226809+ 1.6270259788)/ histos_sc->Integral();
	histos_sc->Scale(scscale);			
	histos_sc->Smooth();

	
	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_sc;
	vector<double> mva_error_sc;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_sc->GetBinContent(i, j);
		double binerror = histos_sc->GetBinError(i, j);
		double error = 0;
//		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_sc.push_back(bincontent);
		mva_error_sc.push_back(error);
//		cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " eff error: " << bincontent*error << " rel error: " << binerror/bincontent << endl;
		//				cout << hmlp[j]->GetBinContent(i) << " +/-" << hmlp[j]->GetBinContent(i)/sqrt(hmlp[j]->GetBinError(i))  << endl;
	}
	}	
	
	TRandom myRandom1(12345);

	cout << "PS: " << histos_ps->Integral() << endl;
	cout << "SC: " << histos_sc->Integral() << endl;
	cout << "ZZ: " << histos_zz->Integral() << endl;
	cout << "ZX: " << histos_zx->Integral() << endl;


	double iter = 10000;

	TH1F* hratio1 = new TH1F("lhratio1", "lhratio1", 100, -10, 10);   
	TH1F* hratio2 = new TH1F("lhratio2", "lhratio2", 100, -10, 10);   

	for(Int_t l=0; l<iter; l++) {
	
		// Data counts:
		vector<double> D1;
		vector<double> D2;
		vector< vector<double> > A1;
		vector< vector<double> > A2;
		vector< vector<double> > f1;
		vector< vector<double> > f2;
		vector<double> p;
		
	
		// Make input vectors for Poisson-Gamma model
		for(int i=1; i < ybins+1; i++)
		{
		for(int j=1; j < xbins+1; j++)
		{
	//		cout << "signal[" << thismass << "]: " << sigNames[thismass] << "\t" << signal_weights[thismass] << endl;
			double pdata_bin = myRandom1.Poisson(histos_zz->GetBinContent(i,j)) + myRandom1.Poisson(histos_zx->GetBinContent(i,j)) + myRandom1.Poisson(histos_ps->GetBinContent(i,j));
	//		double pdata_bin = double(histos_dat->GetBinContent(i,j));
			D1.push_back(pdata_bin);			
		}
		}

		// Make input vectors for Poisson-Gamma model
		for(int i=1; i < ybins+1; i++)
		{
		for(int j=1; j < xbins+1; j++)
		{
	//		cout << "signal[" << thismass << "]: " << sigNames[thismass] << "\t" << signal_weights[thismass] << endl;
			double pdata_bin = myRandom1.Poisson(histos_zz->GetBinContent(i,j)) + myRandom1.Poisson(histos_zx->GetBinContent(i,j)) + myRandom1.Poisson(histos_sc->GetBinContent(i,j));
	//		double pdata_bin = double(histos_dat->GetBinContent(i,j));
			D2.push_back(pdata_bin);			
		}
		}


		// prepare pseudo-scalar vectors
		A1.push_back(mva_src_ps);
		f1.push_back(mva_error_ps);

		A1.push_back(mva_src_zz);
		f1.push_back(mva_error_zz);	

		A1.push_back(mva_src_zx);
		f1.push_back(mva_error_zx);	


		// prepare scalar vectors
		A2.push_back(mva_src_sc);
		f2.push_back(mva_error_sc);

		A2.push_back(mva_src_zz);
		f2.push_back(mva_error_zz);	

		A2.push_back(mva_src_zx);
		f2.push_back(mva_error_zx);	


		p.push_back( 1.0 );
		p.push_back( 1.0 );
		p.push_back( 1.0 );

		double ps1 = pg::poissongamma(D1, p, A1, f1, 1.0, true, false);
		double sc1 = pg::poissongamma(D1, p, A2, f2, 1.0, true, false);

		double ps2 = pg::poissongamma(D2, p, A1, f1, 1.0, true, false);
		double sc2 = pg::poissongamma(D2, p, A2, f2, 1.0, true, false);

		double lhratio1 = -2*log(exp(ps1)/exp(sc1));
		double lhratio2 = -2*log(exp(ps2)/exp(sc2));

//		cout << l << "  sc: " << sc1 << "\t ps: " << ps1 << "\t ratio1: " << lhratio1 << "\t ratio2: " << lhratio2 << endl;
		hratio1->Fill(lhratio1);
		hratio2->Fill(lhratio2);

	}

	TString nomeOUT = "output/lhratio.pdf";
	TCanvas *c3 = new TCanvas("lhratio", "lhratio", 400, 400);
	c3->cd();
	hratio1->SetLineColor(4);
	hratio1->SetTitle("Log Likelihood Ratio");
	hratio1->GetXaxis()->SetTitle("log(L_{sc}/L_{s})");
  	hratio1->Draw();

	hratio2->SetLineColor(2);
	hratio2->SetTitle("Log Likelihood Ratio");
	hratio2->GetXaxis()->SetTitle("-2*log(L_{sc}/L_{s})");
  	hratio2->Draw("same");

	TString nomeOUT2 ="output/lhratio.png";
   	c3->SaveAs(nomeOUT2);
   	
}
