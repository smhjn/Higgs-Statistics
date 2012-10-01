//-----------------------------------------------------------------------------
// File:        data_likelihood.C
// Description: Generate likelihood function in the form of a TGraph using
//		the data histogram.
//				discriminant histogram
//-----------------------------------------------------------------------------
// Created:     11-04-2012  Joe Bochenek

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
#include <TMath.h>
#include <TTree.h>
#include <TLatex.h>
#include <TROOT.h>
#include <TFile.h>
#include <iostream>
#include <iomanip>
#include <TGaxis.h>
#include <gsl/gsl_randist.h>
#include <TGraph2D.h>

using namespace std;
#include "samples.h"

//#include "samples.h"
#include "IntrinsicPrior.h"
#include "PoissonGammaIntegral.h"



double Min( double x, double y )
{
    return ( x > y ) ? y : x;
}


int combinations_i(std::vector< std::vector<int> > &data, std::vector<int> & histogram, int total, int bin, int totalbins, int bincap, int count)
{
		if( bin > (totalbins - 2) ) {
			if(total < bincap + 1) {
				histogram[bin] = total;
				data.push_back(histogram);
				count++;
				if(!(count%1000000)) cout << count << endl;
//				cout << "bin: " << bin << "/" << totalbins << "\t total: " << total << "\tcount: " << count << endl;
			}
			return count;
		} else
		{
			for(int thisbincount = Min(bincap, total); thisbincount >= 0; thisbincount--){
				histogram[bin] = thisbincount;
//				cout << "\t" <<  thisbincount;
				count = combinations_i(data, histogram, (total-thisbincount), (bin+1), totalbins, bincap, count);
			}
			return count;
		}
}


int combinations(std::vector< std::vector<int> > &data, int max, int totalbins, int bincap, int count)
{
	for(int total = 0; total < max + 1; total++)
	{
//		cout << "total: " << total << endl;
		std::vector<int> histogram(totalbins);
		count = combinations_i(data, histogram, total, 0, totalbins, bincap, count);		
//		cout << "===========" << endl;
	}
	return count;
}




int main(){
	//TString html_name = "/afs/cern.ch/user/j/jpb/www/studies/HWW_analysis/plots/index.html";  
	//ofstream event_list;
	//event_list.open(html_name, ios_base::trunc);
	
	
	//---------------------------------------------------------------------------
	// Make Vectors/Histograms
	//---------------------------------------------------------------------------

	const int xsec_plots = 50;
	double lumi_scale = 1.0;

	TString scripttag = "input1D";	
	
	double total_uncer[2] = {
		5.1/64.5,
		10.0/46.0,
	};
	
	style();
	TFile f3("output/hist_ensemble.root");
	f3.cd();
	
	//---------------------------------------------------------------------------
	// BOOK HISTOGRAMS
	//---------------------------------------------------------------------------

	// Make signal histograms
	double signal_weights[maxsig_plots] = {0};
	TH2F* histos_sig[maxsig_plots];
//	TFile f5("output/hist_ensemble_morph.root");
//	f5.cd();
	
	// The number of signal files available from morphing algorithm
	for(Int_t k=0; k < maxsig_plots; k++) 
	{
		std::stringstream _name;
//		sigNames_mc[k] = double(2*k + 116);
		_name << "histos_sig_" << sigNames_mc[k]; 
		TString name = _name.str();
		cout << name << endl;
		histos_sig[k]= (TH2F*)gDirectory->Get(name);
		histos_sig[k]->Scale(lumi_scale);
		histos_sig[k]->RebinY(5);
		cout << "Sig Bins: " <<	histos_sig[k]->GetYaxis()->GetNbins() << endl;
		signal_weights[k] = histos_sig[k]->Integral();
		cout << sigNames_mc[k] << ": " << signal_weights[k] << endl;
	}
	
	f3.cd();
	// Make background histograms
	TH2F *histos_bkg=(TH2F*) gDirectory->Get("bkg2d");
	histos_bkg->RebinY(5);
	cout << "Bkg Bins: " <<	histos_bkg->GetYaxis()->GetNbins() << endl;

	histos_bkg->SetDirectory(0);
	histos_bkg->Sumw2();
	histos_bkg->GetXaxis()->SetTitle("BNN(x)");
	histos_bkg->GetYaxis()->SetTitle("m_{4l}");
	vector<double> mva_src_b;
	vector<double> mva_error_b;			
	// BACKGROUND
	for(int i=1; i < xbins+1; i++)
	{
	for(int j=1; j < 5+1; j++)
	{
		double bincontent = histos_bkg->GetBinContent(i, j);
		double binerror = histos_bkg->GetBinError(i, j);		
		double error = 0;
		if(bincontent > 0) error = double(binerror*binerror/bincontent);	
		cout << "bin " << j << "\tcontent :" << bincontent << "\terror " << binerror <<  "\trelative: " << 100* binerror/bincontent << "\terror: " << error << endl;
		mva_src_b.push_back(bincontent);
		error = 1;
		mva_error_b.push_back(error);
	}
	}
	histos_bkg->Scale(lumi_scale);
	double weighted = histos_bkg->Integral();
	


	//---------------------------------------------------------------------------
	// Remake  A, f vectors for IntrinsicPrior algorithm
	//---------------------------------------------------------------------------

	vector<double> signals;		

	
	// Loop through m_H and \mu to calculate 2D prior and store the grid in a TGraph2D
	double xsec_factors[xsec_plots]; 
	double count = 0;

    TGraph2D *lh_tgraph = new TGraph2D();
    TGraph2D *lh_tgraph_norm = new TGraph2D();

	clock_t start, end;
	start = clock();



	int bins = 5;
	vector< vector<int> > data;
	combinations(data, 30, bins, 6, 0);
	int size = data.size();
	cout << "Size: " << size << "\t COunt: "<< endl;
	double likelihood_ave[maxsig_plots][xsec_plots];
	for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++)  likelihood_ave[l][j] = 0;
	double final_integral = 0;
	double prob_null_integral = 0;

	// Loop through data permutations
	for(Int_t i_d=0; i_d<size; i_d++)
	{
		double likelihood_sig[maxsig_plots][xsec_plots]; 
		double likelihood_null[maxsig_plots][xsec_plots]; 

		double l_integral_sig = 0;
		vector<double> D;
		int this_size = 0;
		for(Int_t k=0; k<bins; k++)
		{
			D.push_back( double(data[i_d][k]) );	
			this_size += data[i_d][k];
//			cout << data[i_d][k] << endl;
		}

		// First loop through signal mass points
		for(Int_t l=0; l<maxsig_plots; l++) {
//			cout << endl;					
			// Next loop through signal strength modifier points
			for( int j = 0; j < xsec_plots; j++)
			{						
				double xsecfactor =  (j+1)*0.1;
				xsec_factors[j] = xsecfactor;

				vector< vector<double> > A;
				vector< vector<double> > f;
				vector<double> p;		
				vector<double> mva_src_s;
				vector<double> mva_error_s;

//				cout << "i_d: " << i_d << "/" << size << "\tSignal: m_H = " << sigNames_mc[l] << " GeV,\t mu = " << xsecfactor;

				// Build signal vectors to pass to intrinsicprior at this mass point
				// SIGNAL
				for(int i=1; i < xbins+1; i++)
				{
				for(int j=1; j < 5+1; j++)
				{
					double bincontent = histos_sig[l]->GetBinContent(i,j);
					double binerror = histos_sig[l]->GetBinError(i,j);				
					double error = 0;
					if(bincontent > 0) error = double(binerror*binerror/bincontent);		
					mva_src_s.push_back(bincontent * xsecfactor);
					error = 1;
					mva_error_s.push_back(error);
				}
				}
				
				A.push_back(mva_src_s);
				f.push_back(mva_error_s);
	
				A.push_back(mva_src_b);
				f.push_back(mva_error_b);	
	
				p.push_back( 1.0 );
				p.push_back( 1.0 );
	
				double prob_sig  = exp(pg::poissongamma(D, p, A, f, true, true));
//				double prob_sig = 0;			
				likelihood_sig[l][j] = prob_sig;
				l==0 ? l_integral_sig += prob_sig * 5.0 * 0.1 : l_integral_sig += prob_sig * 10.0 * 0.1;
				lh_tgraph->SetPoint( (j+l*xsec_plots), sigNames_mc[l], xsecfactor, prob_sig );
				p.clear();
				f.clear();		
				A.clear();
			}
		}

		// Do NULL (mu = 0) likelihood
		vector< vector<double> > A2;
		vector<double> p2;
		vector< vector<double> > f2;

		A2.push_back(mva_src_b);
		f2.push_back(mva_error_b);
		p2.push_back( 1.0 );

//		cout << "\npg: " << endl;
		double prob_null = exp(pg::poissongamma(D, p2, A2, f2, true, true));
//		double prob_null = 1.0;
		
//		cout << "\ngsl: " << endl;
//		double cdf = 1.0;
//		for( int k = 0; k < 5; k++){
//			double thisp = (0.5);///f2[0][k]);
//			double thisA = A2[0][k]; //*f2[0][k];
//			cout <<	"D[" << k << "] = " << D[k] << "\t" << "A[" << k << "] = " << thisA << "\t" << "p[" << k << "] = " << thisp << "\t" <<  gsl_ran_negative_binomial_pdf(D[k], thisp, thisA) << endl;
//			cdf = cdf * gsl_ran_negative_binomial_pdf(D[k], thisp, thisA );
//		}
//		if(!(i_d%10000)) cout << "i_d: " << i_d << "\tsize: " << this_size << ",\tpg = " << prob_null << " gsl: " << cdf << "\t ratio: " << 100*(1 - cdf / prob_null) <<  endl;

		for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++) likelihood_sig[l][j] /= ( l_integral_sig );
		for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++) likelihood_ave[l][j] += prob_null * likelihood_sig[l][j] / 2.2194;

		if(!(i_d%1000)) cout << "i_d: " << i_d << "/" << size << "\t Likelihood = " << prob_null << "\t sig lh int: " << prob_null_integral << endl;
		prob_null_integral += prob_null;
		D.clear();
	}
	
	
	for(Int_t l=0; l < maxsig_plots; l++) {cout <<  "sigNames_mc[" << l << "]:  " << sigNames_mc[l] << "\t";
	for(Int_t j=1; j < xsec_plots; j++) cout << likelihood_ave[l][j] << "\t";
	cout << endl;
	}

	for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++) lh_tgraph_norm->SetPoint( ((j-1)+l*xsec_plots), sigNames_mc[l], xsec_factors[j], likelihood_ave[l][j] );
//	for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++) cout << ((j-1)+l*xsec_plots) << "\t" << sigNames_mc[l] << "\t" << xsec_factors[j] << "\t" << likelihood_ave[l][j] <<  endl;
	for(Int_t l=0; l < maxsig_plots; l++) for(Int_t j=0; j < xsec_plots; j++) l==0 ? final_integral += likelihood_ave[l][j] * 5.0 * 0.1 : final_integral += likelihood_ave[l][j] * 10.0 * 0.1;
	cout << "null dist integral: " << prob_null_integral << "\tfinal integral: " << final_integral << endl;

	count+= 1.0;
	end = clock();
	
	cout << "Time required for execution: "
	<< (double)(end-start)/CLOCKS_PER_SEC
	<< " seconds." << "\n\n";

	std::stringstream name;
	name << "intrinsic_prior_norm";

	std::stringstream plotname;
	plotname << "output/" << name.str() << ".png";

	std::stringstream title;
	title << "Intrinsic Prior, m_{H} vs #mu";

	// Make plot of likelihood
	TCanvas *c = new TCanvas(name.str().c_str(), title.str().c_str(), 200,10,700,500);
	lh_tgraph_norm->SetTitle(title.str().c_str());
    lh_tgraph_norm->GetXaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetXaxis()->SetTitle("m_{H} (GeV)");
    lh_tgraph_norm->GetZaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetZaxis()->SetTitle("#pi_{0}(m_{H},#mu)");
    lh_tgraph_norm->GetYaxis()->SetTitleFont(42);
	lh_tgraph_norm->GetYaxis()->SetTitle("#mu");

	lh_tgraph_norm->Draw("colz");
//	lh_tgraph_norm->SetMinimum(0.0);
//	lh_tgraph_norm->SetMinimum(0.001);

	c->Update();
	c->Modified();
	c->SaveAs( plotname.str().c_str() );
	delete c;

	name.str("");
	plotname.str("");
	title.str("");
}
