//-----------------------------------------------------------------------------
// File:        data_likelihood.C
// Description: Generate likelihood function in the form of a TGraph using
//		the data histogram.
//				discriminant histogram
//-----------------------------------------------------------------------------
// Created:     11-04-2012  Joe Bochenek

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
#include "samples.h"

//#include "samples.h"
#include "PoissonGammaFit2.h"
#include "PoissonGammaIntegral.h"

int main(int argc, char* argv[]){
	//TString html_name = "/afs/cern.ch/user/j/jpb/www/studies/HWW_analysis/plots/index.html";  
	//ofstream event_list;
	//event_list.open(html_name, ios_base::trunc);
	
	int num = 0;

	if (argc > 1){
	if (std::string(argv[1]) == "-n") {
			if (2 < argc) { // Make sure we aren't at the end of argv!
				if (int num = atoi(argv[2])) { // Make sure we aren't at the end of argv!
					cout << "TGraph name " <<  num << endl;
					cout << "argv: " << argv[0] << "\t" << argv[1] << "\t" << num << endl; 
				} else {
				  std::cerr << "-r option requires an integer." << std::endl;	
				  std::cerr << "Usage: " << argv[0] << "-n somenumber (integer)" << std::endl;
				  return 1;			
				}
			} else { // Uh-oh, there was no argument
				  std::cerr << "-n option requires one argument." << std::endl;
        		  std::cerr << "Usage: " << argv[0] << "-n somenumber (used for batch jobs)" << std::endl;
				  return 1;
			}  
	} }
	
	//---------------------------------------------------------------------------
	// Make Vectors/Histograms
	//---------------------------------------------------------------------------

	const int xsec_plots = 50;
	double lumi_scale = 1.0;
	double this_xsec_factor = 1.0;
	int thismass = 2;

	TString scripttag = "input1D";	
	style();
	
	//---------------------------------------------------------------------------
	// RETRIEVE HISTOGRAMS
	//---------------------------------------------------------------------------

	TH2F* histos_sig[maxsig];

	TFile file("output/hist_noiso.root");
	file.cd();

	// ZZ background
	TH2F *histos_bkg=(TH2F*) gDirectory->Get("bkg2d");
	histos_bkg->SetDirectory(0);
	histos_bkg->Sumw2();
	histos_bkg->GetXaxis()->SetTitle("BNN(x)");
	histos_bkg->GetYaxis()->SetTitle("m_{4l}");
	double signal_weights[maxsig] = {0};

	// Data driven reducible background
	TH2F *histos_zjets=(TH2F*) gDirectory->Get("zjets2d_reweight");
	histos_zjets->SetDirectory(0);
	histos_zjets->Sumw2();
	histos_zjets->GetXaxis()->SetTitle("BNN(x)");
	histos_zjets->GetYaxis()->SetTitle("m_{4l}");


	// Data driven reducible background
	TH2F *histos_zjets_sel=(TH2F*) gDirectory->Get("all_redecible_bkg");
	histos_zjets_sel->SetDirectory(0);
	histos_zjets_sel->Sumw2();
	histos_zjets_sel->GetXaxis()->SetTitle("BNN(x)");
	histos_zjets_sel->GetYaxis()->SetTitle("m_{4l}");


	// Data
	TH2F *histos_dat=(TH2F*) gDirectory->Get("dat2d");
	histos_dat->SetDirectory(0);	
	histos_dat->Sumw2();
	histos_dat->GetXaxis()->SetTitle("BNN(x)");
	histos_dat->GetYaxis()->SetTitle("m_{4l}");

	// Print some histograms for debugging
	
	// pseudodata averaged input histogram
	string name = "bnnvsmt_dat";
	std::stringstream _plotname;
	_plotname << "output/all/" << name << "_data.pdf"; 
//	histos_sig[thismass]->Scale( signal_weights[thismass]/histos_sig[thismass]->Integral() );
	TCanvas *c = new TCanvas(name.c_str(),"Average Likelihood vs. Mass",200,10,700,500);
	histos_dat->Draw("colz");
	c->Update();
	c->Modified();
	TString plotname = _plotname.str();		
	cout << plotname << endl;
	c->SaveAs(plotname);
	delete c;



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
//		histos_sig[k]->Rebin(25);
		histos_sig[k]->Scale(lumi_scale);
		signal_weights[k] = histos_sig[k]->Integral();
		cout << sigNames[k] << ": " << signal_weights[k] << endl;
		cout << "signal[" << k << "]: " << sigNames[k] << "\t" << signal_weights[k] << endl;
		histos_sig[k]->Sumw2();
	}
	
	// print a reference histogram

	char filename_s[256];
	char charname_s[] = "bnnvsmt_sig";
	sprintf(filename_s, "output/all/%s_data.pdf", charname_s);
	c = new TCanvas(charname_s,"Average Likelihood vs. Mass",200,10,700,500);
	cout << "signal integral: " << histos_sig[thismass]->Integral() << endl;
	histos_sig[thismass]->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename_s << endl;
	c->SaveAs(filename_s);
	delete c;
	
	



	//---------------------------------------------------------------------------
	// FILL DATA VECTOR
	//---------------------------------------------------------------------------
	
	vector<double> mva_src_d;
	vector<double> mva_error_d;			

	// BACKGROUND
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_dat->GetBinContent(i, j);
		mva_src_d.push_back(bincontent);
		//				cout << "bin " << i << " bin content " << bincontent << " stat error: " << binerror << " exp error: " << exp_error << " exp error: " << tot_error  << endl;
	}
	}


	//---------------------------------------------------------------------------
	// FILL BACKGROUND VECTORS
	//---------------------------------------------------------------------------

	// Calculate errors in bins for the backgrounds
	histos_bkg->Scale(lumi_scale);			
	histos_bkg->Smooth();

	vector<double> mva_src_b;
	vector<double> mva_error_b;			
	// BACKGROUND
	// if(j>0)	cout << "\n" << names[j] << "\t p:" << p[j] << endl;
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
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


	cout << "Z+jets weight: " << histos_zjets_sel->Integral()/histos_zjets->Integral() << endl;
	histos_zjets->Scale(1.233/histos_zjets->Integral());
	histos_zjets->Scale(lumi_scale);

	char filename[256];
	char charname [] = "bnnvsmt_bkg";
	sprintf(filename, "output/all/%s_data.pdf", charname);
	c = new TCanvas(charname,"Average Likelihood vs. Mass",200,10,700,500);
	histos_bkg->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename << endl;
	c->SaveAs(filename);
	delete c;

			
	// Calculate errors in bins for the backgrounds
	vector<double> mva_src_zjets;
	vector<double> mva_error_zjets;			
	// BACKGROUND
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
		double bincontent = histos_zjets->GetBinContent(i, j);
		double binerror = histos_zjets->GetBinError(i, j);
		double error = 0;
		if(bincontent > 0) error = double(bincontent/(binerror*binerror) );	
		mva_src_zjets.push_back(bincontent);
		mva_error_zjets.push_back(error);
	}
	}
	
	cout << "Z+jets Integral: " << histos_zjets->Integral() << endl;
	cout << "ZZ Integral: " << histos_bkg->Integral() << endl;
	cout << "Signal: " << histos_sig[12]->Integral() << endl;
	cout << "Data Integral: " << histos_dat->Integral() << endl;


	char filename_z[256];
	char charname_z [] = "bnnvsmt_zjets";
	sprintf(filename_z, "output/all/%s_bkg.pdf", charname_z);
	c = new TCanvas(charname_z,"Average Likelihood vs. Mass",200,10,700,500);
	histos_zjets->Draw("colz");
	c->Update();
	c->Modified();
	cout << filename_z << endl;
	c->SaveAs(filename_z);
	delete c;

	histos_zjets->Smooth();
	histos_zjets->Smooth();
	histos_zjets->Smooth();
	histos_zjets->Smooth();
	histos_zjets->Smooth();


	//---------------------------------------------------------------------------
	// START LOOP OVER 2D likelihood
	//---------------------------------------------------------------------------


	// Load Reference prior
	const int Npoints = 1500;
	// Get prior as posterior of previous
	TFile rpfile("../HZZ4l/output/systs.root");
	TGraph2D *refprior = (TGraph2D*) gDirectory->Get("likelihood_syst_0");
	Double_t rppoints[Npoints];
	Double_t *rpPoints = refprior->GetZ();
	Double_t *priorx = refprior->GetX();
	Double_t *priory = refprior->GetY();



	double ave_significance = 0.;
	double ave_significance_ref = 0.;
	double ave_significance_prof = 0.;
	double ave_significance_prof2 = 0.;

	int iter = 100;
	for(Int_t k=0; k<iter; k++) {

	TRandom myRandom1(k);
	TRandom myRandom2(2*(k+10));


    TGraph2D *dt = new TGraph2D();
    TGraph2D *dt_post = new TGraph2D();
    TGraph2D *dt2 = new TGraph2D();
    TGraph2D *dt_null = new TGraph2D();

   // draw a frame to define the range
	TLegend leg(0.8,0.5,1.0,0.9);

	// Data counts:
	vector<double> D;
			
	double xsec_factors[xsec_plots]; 

	// Make input vectors for Poisson-Gamma model
	for(int i=1; i < ybins+1; i++)
	{
	for(int j=1; j < xbins+1; j++)
	{
//		cout << "signal[" << thismass << "]: " << sigNames[thismass] << "\t" << signal_weights[thismass] << endl;
//		double pdata_bin = myRandom1.Poisson(histos_bkg->GetBinContent(i,j) + histos_zjets->GetBinContent(i,j) + this_xsec_factor * histos_sig[thismass]->GetBinContent(i,j));
		double pdata_bin = myRandom1.Poisson(histos_zjets->GetBinContent(i, j)) + myRandom1.Poisson(histos_bkg->GetBinContent(i,j)) + myRandom1.Poisson(this_xsec_factor * histos_sig[thismass]->GetBinContent(i,j));
//		double pdata_bin = double(histos_dat->GetBinContent(i,j));
		D.push_back(pdata_bin);			
	}
	}

	double normlikelihood[maxsig][xsec_plots];
	double likelihood[maxsig][xsec_plots]; 
	double significance[maxsig][xsec_plots]; 
	double likelihood_null[maxsig][xsec_plots]; 
	double posterior[maxsig][xsec_plots];


	double l_integral = 0;
	double p_integral_ref = 0;
	double p_integral_flat = 0;
	double rp_integral = 0;
	
	for(Int_t l=0; l<maxsig; l++) {
		for( int j = 0; j < xsec_plots; j++)
		{						
			double xsecfactor =  double(j)*0.1;
			xsec_factors[j] = xsecfactor;

			// Scale factor
			vector< vector<double> > A;
			vector< vector<double> > f;
			vector<double> p;
			vector<double> mva_src_s;
			vector<double> mva_error_s;

			// SIGNAL
			for(int i=1; i < ybins+1; i++)
			{
			for(int j=1; j < xbins+1; j++)
			{
				double bincontent = xsecfactor * histos_sig[l]->GetBinContent(i,j);
				double binerror = xsecfactor * histos_sig[l]->GetBinError(i,j);				
				double error = 0. ;
				if(bincontent < 0) bincontent = 0;
//				if(bincontent > 0) error = double( binerror*binerror/bincontent );		
				if(bincontent > 0) error = 10 / bincontent;

//				if (l==15 && i==10) cout << sigNames[l] << "/" << ymin + i << "\txsecfactor: "  << xsecfactor << "\tbincontent: " << bincontent << "\t eff events: " << bincontent * 25000 << "\t error: " << error << "\t rel error: " << 1/sqrt(25000) << endl;				
				mva_src_s.push_back(bincontent);
				mva_error_s.push_back(error);
			}
			}
			
//			exit(0);
				
			A.push_back(mva_src_s);
			f.push_back(mva_error_s);

			A.push_back(mva_src_b);
			f.push_back(mva_error_b);	

			A.push_back(mva_src_zjets);
			f.push_back(mva_error_zjets);	

			p.push_back( 1.0 );
			p.push_back( 1.0 );
			p.push_back( 1.0 );
						
//			cout << "Print Mass\t" << sigNames[l] << "\t xsec: " << xsecfactor <<endl;


/*
			double prob = 1;
			for(int i=1; i < xbins+1; i++)
			{
			for(int j=1; j < ybins+1; j++)
			{
				double binsum = xsecfactor *  histos_sig[l]->GetBinContent(i,j) + histos_bkg->GetBinContent(i, j);
				double Di = double(histos_dat->GetBinContent(i,j));
//				double Di = histos_bkg->GetBinContent(i,j) + this_xsec_factor * histos_sig[thismass]->GetBinContent(i,j);
				double fact;
				if(Di > 0) { fact = gsl_sf_gamma(Di+1);} else {fact = 1.;}
			    prob *= pow(binsum, Di) * exp(-binsum) / fact;
//				cout << sigNames[l] << "\txsecfactor: "  << xsecfactor << "\t Di: " << Di << "\t binsum: " << binsum << "\tprob " << prob << endl;

			}
			}
			likelihood[l][j] = prob;
*/			

			double prob = pg::poissongamma(D, p, A, f, 1.0, true, false);
			
			likelihood[l][j] = exp(prob);
//			cout << "Prob Here " << likelihood[l][j] << endl;
//			posterior[l][j] = likelihood[l][j] * refprior->Interpolate(sigNames[l], xsecfactor );

			// flat prior
			posterior[l][j] = 0.1 * likelihood[l][j] / 250;

			l_integral += likelihood[l][j] * 0.1;
			p_integral_ref  += likelihood[l][j] * refprior->Interpolate(sigNames[l], xsecfactor ); // rpPoints[(j+l*xsec_plots)];
			p_integral_flat += likelihood[l][j] / 250;
			rp_integral += rpPoints[(j+l*xsec_plots)] * 0.1;

//			if(l<maxsig) l_integral += exp(prob) * (sigNames[l+1] - sigNames[l]) * 0.1;
			f.clear();		
			A.clear();
		}
		
	}



	//---------------------------------------------------------------------------
	// Calculate Significance
	//---------------------------------------------------------------------------
	// Do NULL (mu = 0) likelihood
	vector< vector<double> > A2;
	vector<double> p2;
	vector< vector<double> > f2;

	A2.push_back(mva_src_b);
	f2.push_back(mva_error_b);	

	A2.push_back(mva_src_zjets);
	f2.push_back(mva_error_zjets);	

	p2.push_back( 1.0 );
	p2.push_back( 1.0 );
	
	double prob_null = exp(pg::poissongamma(D, p2, A2, f2, 1.0, true, false));

	cout << "Null Hypothesis: " << prob_null << endl;
	cout << "Ref Prior Int: " << rp_integral << endl;
	cout << "posterior Int: " << p_integral_flat << endl;
	

	for(Int_t m=0; m < maxsig; m++)for(Int_t j=0; j < xsec_plots; j++)  {
		double signif = 0;
		likelihood[m][j] > prob_null ? signif = sqrt(2*log(likelihood[m][j]/prob_null)) : signif = 0;
//		if(likelihood[m][j] > prob_null ) cout << likelihood[m][j] << " / " << l_integral << "\tnull: " << prob_null << "\t signif" << signif   << endl ;
		dt2->SetPoint( (j+m*xsec_plots), sigNames[m], xsec_factors[j], signif );
	}
	

	TGraph *tmass =  new TGraph();
	for(Int_t j=0; j < xsec_plots; j++)  { tmass->SetPoint(tmass->GetN(), xsec_factors[j], likelihood[thismass][j]); }	
	
	for(Int_t m=0; m < maxsig-1; m++)for(Int_t j=1; j < xsec_plots; j++)  dt->SetPoint( ((j-1)+m*(xsec_plots-1)), sigNames[m], xsec_factors[j], likelihood[m][j]);
	double maxsignif = dt->GetZmax();
	
	double signif = sqrt(2*log(p_integral_flat / ( prob_null)));
	double signif_ref = sqrt(2*log(p_integral_ref  / ( rp_integral * prob_null))) ;
	double signif_prof = sqrt(2*log(maxsignif / ( prob_null)));
	double signif_prof2 = sqrt(2*log(tmass->GetHistogram()->GetMaximum()/ ( prob_null)));

	cout << "iter: " << k << endl;

	cout << tmass->GetHistogram()->GetMaximum() << endl;
	cout << "Global Significance (ref  prior): " << sqrt(2*log(p_integral_ref  / ( rp_integral * prob_null))) << endl;
	cout << "Global Significance (flat prior): " << sqrt(2*log(p_integral_flat / ( prob_null))) << endl;

	if(signif>-1) ave_significance += signif;
	if(signif_ref>-1) ave_significance_ref += signif_ref;
	if(signif_prof>-1) ave_significance_prof += signif_prof;
	if(signif_prof2>-1) ave_significance_prof2 += signif_prof2;

	cout << "Average Signif. (so far) = " << ave_significance / (k + 1) << endl;
	cout << "Average ref Signif. (so far) = " << ave_significance_ref / (k + 1) << endl;
	cout << "Average prof Signif. (so far) = " << ave_significance_prof / (k + 1) << endl;
	cout << "Average prof2 Signif. (so far) = " << ave_significance_prof2 / (k + 1) << endl;


	cout << "-----------------------------------------------" << endl;



	}
	
	exit(0);

//	for(Int_t m=0; m < maxsig; m++)for(Int_t j=0; j < xsec_plots; j++)  likelihood[m][j] /= l_integral;
//	for(Int_t m=0; m < maxsig; m++){
//		for(Int_t j=0; j < xsec_plots; j++)  cout << likelihood[m][j] << "\t";
//		cout << endl;
//	}



}