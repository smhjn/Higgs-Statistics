{


TFile *_file0 = TFile::Open("/home/jbochenek/work/HZZ4l_2012/output/hist_noiso.root");

// Load Z+X SS Data, selected for isolation
TH2F *zjets = (TH2F*) gDirectory->Get("zjets2d");
TH1F *zjets_y = (TH1F*) zjets2d->ProjectionY();
TH1F *zjets_x = (TH1F*) zjets2d->ProjectionX();

// Load Z+X SS Data, not selected for isolation
TH2F *zjets_sel = (TH2F*) gDirectory->Get("zjets2d_sel");
TH1F *zjets_sel_y = (TH1F*) zjets2d_sel->ProjectionY();
TH1F *zjets_sel_x = (TH1F*) zjets2d_sel->ProjectionX();

// Load Z+X SS Data, not selected for isolation, reweighted for fakes
TH2F *zjets_reweight = (TH2F*) gDirectory->Get("zjets2d_reweight");
TH1F *zjets_reweight_y = (TH1F*) zjets2d_reweight->ProjectionY();
TH1F *zjets_reweight_x = (TH1F*) zjets2d_reweight->ProjectionX();

// Load reducible MC backgrounds
TH2F *redbkg = (TH2F*) gDirectory->Get("all_redecible_bkg");
TH1F *redbkg_y = (TH1F*) redbkg->ProjectionY();
TH1F *redbkg_x = (TH1F*) redbkg->ProjectionX();

// Load combined MC backgrounds
TH2F *allbkg = (TH2F*) gDirectory->Get("bkg2d");
TH1F *allbkg_y = (TH1F*) bkg2d->ProjectionY();
TH1F *allbkg_x = (TH1F*) bkg2d->ProjectionX();

// Load Higgs 125 background
TH2F *sig_125 = (TH2F*) gDirectory->Get("histos_sig_125");
TH1F *sig_125_y = (TH1F*) sig_125->ProjectionY();
TH1F *sig_125_x = (TH1F*) sig_125->ProjectionX();

// Load ZZ background
TH2F *zz = (TH2F*) gDirectory->Get("bkg2d");
TH1F *zz_y = (TH1F*) zz->ProjectionY();
TH1F *zz_x = (TH1F*) zz->ProjectionX();

// Load data
TH2F *dat = (TH2F*) gDirectory->Get("dat2d");
TH1F *dat_y = (TH1F*) dat->ProjectionY();


cout << zjets_reweight_y->Integral() << endl;
cout << redbkg_y->Integral() << endl;



// ---------------------------------------------------
//  bnn_bkgs_sig.pdf - Weighted BNN plot
// ---------------------------------------------------

TLegend legend(0.6, 0.65, 0.95, 0.95);
TCanvas *c1 = new TCanvas("c1","c1",600,500);
gStyle->SetOptStat(0);

scalefactor4 = 1/redbkg_y->Integral();
redbkg_y->SetTitle("BNN Output, Normalized");
redbkg_y->Scale(scalefactor4);
redbkg_y->SetLineColor(4);
redbkg_y->SetLineWidth(2);
redbkg_y->SetMarkerColor(4);
redbkg_y->Rebin(2);
redbkg_y->Smooth(20);
redbkg_y->Draw("hist e");
redbkg_y->SetLineColor(4);

double scalefactor2 = 1/zjets_y->Integral();
zjets_y->Scale(scalefactor2);
zjets_y->SetLineColor(1);
zjets_y->SetLineWidth(2);
zjets_y->SetTitle("BNN Output, Normalized");
zjets_y->GetXaxis()->SetTitle("BNN(x)");
zjets_y->Rebin(2);
zjets_y->Smooth(20);
zjets_y->Draw("hist e same");
zjets_y->SetLineColor(1);

double scalefactor2 = 1/zjets_sel_y->Integral();
zjets_sel_y->Scale(scalefactor2);
zjets_sel_y->SetLineColor(3);
zjets_sel_y->SetLineWidth(2);
zjets_sel_y->SetTitle("BNN Output, Normalized");
zjets_sel_y->GetXaxis()->SetTitle("BNN(x)");
zjets_sel_y->Rebin(2);
zjets_sel_y->Smooth(20);
zjets_sel_y->Draw("hist e same");
zjets_sel_y->SetLineColor(3);



double scalefactor2 = 1/zjets_reweight_y->Integral();
zjets_reweight_y->Scale(scalefactor2);
zjets_reweight_y->SetLineColor(2);
zjets_reweight_y->SetLineWidth(2);
zjets_reweight_y->SetTitle("BNN Output, Normalized");
zjets_reweight_y->GetXaxis()->SetTitle("BNN(x)");
zjets_reweight_y->Rebin(2);
zjets_reweight_y->Smooth(20);
zjets_reweight_y->Draw("hist e same");
zjets_reweight_y->SetLineColor(2);


legend->AddEntry(zjets_y, "Z+X (unselected) from data", "l");
legend->AddEntry(zjets_sel_y, "Z+X (iso selected) from data", "l");
legend->AddEntry(zjets_reweight_y, "Z+X (reweighted) from data", "l");
legend->AddEntry(redbkg_y, "Z+X from MC", "l");
legend->Draw("same");

c1->SaveAs("../output/mc_vs_data_zjets.pdf");


// ---------------------------------------------------
//  bnn_bkgs_sig.pdf - Weighted BNN plot
// ---------------------------------------------------

TLegend legend2(0.6, 0.65, 0.95, 0.95);
TCanvas *c2 = new TCanvas("c1","c1",600,500);
c2->cd()
gStyle->SetOptStat(0);

redbkg_y->Draw("hist e");
redbkg_y->SetLineColor(3);
zjets_reweight_y->Draw("hist e same");

double scalefactor2 = 1/zz_y->Integral();
zz_y->Scale(scalefactor2);
zz_y->SetLineColor(3);
zz_y->SetLineWidth(2);
zz_y->SetTitle("BNN Output, Normalized");
zz_y->GetXaxis()->SetTitle("BNN(x)");
zz_y->Rebin(2);
zz_y->Draw("hist e same");
zz_y->SetLineColor(4);

legend2->AddEntry(zjets_reweight_y, "Z+X (reweighted) from data", "l");
legend2->AddEntry(zz_y, "ZZ MC", "l");
legend2->AddEntry(redbkg_y, "Reducible background MC", "l");

legend2->Draw("same");

c2->SaveAs("../output/zz_zjets.pdf");




// ---------------------------------------------------
//  2D input plots
// ---------------------------------------------------

TCanvas *c2 = new TCanvas("c2","c2",600,500);

zz->SetTitle("ZZ background - BNN(x) vs. m_{4l}");
zz->GetXaxis()->SetTitle("m_{4l}");
zz->GetYaxis()->SetTitle("BNN(x)");
zz->Draw("colz");

c2->SaveAs("../output/ZZ_2d.pdf");


TCanvas *c3 = new TCanvas("c3","c3",600,500);

zz->Smooth();
zz->SetTitle("ZZ background - BNN(x) vs. m_{4l}");
zz->GetXaxis()->SetTitle("m_{4l}");
zz->GetYaxis()->SetTitle("BNN(x)");
zz->Draw("colz");

c3->SaveAs("../output/ZZ_2d_smooth.pdf");


TCanvas *c4 = new TCanvas("c4","c4",600,500);

zjets->SetTitle("Z+X background - BNN(x) vs. m_{4l}");
zjets->GetXaxis()->SetTitle("m_{4l}");
zjets->GetYaxis()->SetTitle("BNN(x)");
zjets->Draw("colz");

c4->SaveAs("../output/Z_jets.pdf");


TCanvas *c5 = new TCanvas("c5","c5",600,500);

histos_sig_125->SetTitle("Higgs (m_{H} = 125) - BNN(x) vs. m_{4l}");
histos_sig_125->GetXaxis()->SetTitle("m_{4l}");
histos_sig_125->GetYaxis()->SetTitle("BNN(x)");
histos_sig_125->Draw("colz");

c5->SaveAs("../output/signal_example.pdf");



cout << "last one" << endl;

// ---------------------------------------------------
//  mass_plot_comb_norm.pdf - Weighted BNN plot
// ---------------------------------------------------

TLegend legend_x(0.13, 0.69, 0.62, 0.88);
TCanvas *cmass = new TCanvas("cmass","cmass",600,500);
gStyle->SetOptStat(0);

//dat_y->Draw("lep");

//scalefactor4 = zjets_sel_y->Integral()/zjets_y->Integral();
zjets_y->Scale(scalefactor4);
zjets_y->SetLineColor(1);
zjets_y->SetLineWidth(2);
zjets_y->SetTitle("BNN Discriminant");
zjets_y->GetXaxis()->SetTitleOffset(0.65);
zjets_y->GetXaxis()->SetTitleSize(0.06);
zjets_y->GetXaxis()->SetTitle("D(x)");

sig_125_y->SetLineColor(3);
sig_125_y->SetLineWidth(2);
sig_125_y->SetMarkerColor(3);
sig_125_y->SetTitle("BNN Discriminant");

sig_125_y->GetYaxis()->SetTitleOffset(1.);
sig_125_y->GetYaxis()->SetTitleSize(0.04);

sig_125_y->GetXaxis()->SetTitle("D(x)");
sig_125_y->GetYaxis()->SetTitle("a.u. (normalized)");

sig_125_y->GetXaxis()->SetTitleSize(0.06);
sig_125_y->GetXaxis()->SetTitleOffset(0.7);
sig_125_y->Draw("hist");

zz_y->SetLineColor(2);
zz_y->SetLineWidth(2);
zz_y->SetMarkerColor(2);
zz_y->Draw("hist same");
zjets_y->Draw("hist same");

dat_y->SetTitle("m_{4l}");
dat_y->GetXaxis()->SetTitle("m_{4l}");
dat_y->SetLineWidth(2);
//dat_y->Draw("lep same");

legend_x->AddEntry(zjets_y, "Z+X (from data)", "l");
legend_x->AddEntry(zz_y, "ZZ", "l");
legend_x->AddEntry(sig_125_y, "Higgs (m_{H} = 125GeV)", "l");
//legend_x->AddEntry(dat_y, "Data (2012 ICHEP)", "l");
//legend_x->AddEntry(dat, "Data (2011, 5.01 fb^{-1}", "lep");
legend_x->Draw("same");

cmass->SaveAs("../output/mass_plot_comb_norm.pdf");




// ---------------------------------------------------
// Weighted BNN plot
// mass_plot_comb.pdf
// ---------------------------------------------------

TCanvas *cx = new TCanvas("c1","c1",600,500);
gStyle->SetOptStat(0);

TLegend legend_x2(0.6, 0.65, 0.95, 0.85);
scalefactor3 = 1/sig_125_x->Integral();
sig_125_x->Scale(scalefactor3);
sig_125_x->SetLineColor(3);
sig_125_x->SetLineWidth(2);
sig_125_x->SetMarkerColor(3);
sig_125_x->SetTitle("Best Higgs Mass");
sig_125_x->SetTitle("Best Higgs Mass");
sig_125_x->GetXaxis()	->SetTitle("m_{4l}");
sig_125_x->GetXaxis()->SetTitleSize(0.06);
sig_125_x->GetXaxis()->SetTitleOffset(0.65);
sig_125_x->Draw("hist");


zjets_x->Draw("hist same");
//zjets_x->Smooth(50);
scalefactor = 1/allbkg_x->Integral();
allbkg_x->Scale(scalefactor);
allbkg_x->SetLineColor(2);
allbkg_x->SetLineWidth(2);
allbkg_x->SetMarkerColor(2);
allbkg_x->Draw("hist same");

legend_x2->AddEntry(zjets_x, "Z+X", "l");
legend_x2->AddEntry(allbkg_x, "ZZ", "l");
legend_x2->AddEntry(sig_125_x, "Higgs (m_{H} = 125GeV)", "l");
legend_x2->Draw("same");

cx->SaveAs("../output/mass_plot_comb.pdf");


// ---------------------------------------------------

TCanvas *cdat = new TCanvas("cdat","cdat",600,500);
gStyle->SetOptStat(0);
dat2d->SetTitle("Data (2011 dataset, 5.1 fb^{-1}) - BNN(x) vs. m_{4l}");
dat2d->GetXaxis()->SetTitle("m_{4l}");
dat2d->GetYaxis()->SetTitle("BNN(x)");
dat2d->SetMarkerSize(10.);
dat2d->Draw("y:x");
cdat->SaveAs("../output/data.pdf");

// ------------------------------------------------------


TCanvas *cdatsig = new TCanvas("cdatsig","cdatsig",600,500);
gStyle->SetOptStat(0);
histos_sig_125->SetTitle("Data (2011 dataset, 5.1 fb^{-1}) - BNN(x) vs. m_{4l}");
histos_sig_125->GetXaxis()->SetTitle("m_{4l}");
histos_sig_125->GetYaxis()->SetTitle("BNN(x)");
histos_sig_125->Draw("colz");
dat2d->SetTitle("Data (2011 dataset, 5.1 fb^{-1}) - BNN(x) vs. m_{4l}");
dat2d->GetXaxis()->SetTitle("m_{4l}");
dat2d->GetYaxis()->SetTitle("BNN(x)");

dat2d->SetMarkerColor(0);
dat2d->SetFillColor(0);
dat2d->Draw("scat same");
cdatsig->SaveAs("../output/datasig.pdf");

// ------------------------------------------------------


gSystem->Sleep(10);

_file0.Close();
_file1.Close();




}



