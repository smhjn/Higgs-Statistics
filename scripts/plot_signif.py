from string import *
from ROOT import *
from time import sleep

gcomb = TGraph()
gcomb.SetName("gcomb")

gme = TGraph()
gme.SetName("gme")


masses_me = [
117,
120,
125,
130,
135,
140,
]


signif_me = [
1.31722,
1.61111,
2.62013,
3.48262,
4.27191,
6.00181,
]

masses_comb = [
115,
120,
125,
130,
135,
140,
150,
160
]
 
signif_comb = [
1.29787,
1.73086,
2.59454,
3.4906,
4.25986,
5.36156,
5.73365,
3.34258 
]

for i in xrange(len(signif_comb)):
	gcomb.SetPoint(gcomb.GetN(), masses_comb[i], signif_comb[i]);

for i in xrange(len(signif_me)):
	gme.SetPoint(gme.GetN(), masses_me[i], signif_me[i]);


ceff = TCanvas( 'ceff', 'eff', 200, 10, 700, 700 )
ceff.cd();
pad1 = TPad("pad1", "pad1", 0.0, 0.25, 1.0, 1.0);
pad2 = TPad("pad2", "pad2", 0.0, 0.0, 1.0, 0.25);
pad1.SetGridy();
pad1.SetBottomMargin(0)
pad1.SetTopMargin(0)

pad1.Draw()
ceff.Update()

pad1.cd();

legend = TLegend(0.12, 0.7, 0.7, 0.9)
gStyle.SetOptStat(0);
gme.SetLineColor(2);
gme.SetLineWidth(2);
#gme.SetTitle("Local Significance");
gme.GetXaxis().SetTitleSize(0.05);
gme.GetYaxis().SetTitleSize(0.06);
gme.GetXaxis().SetLabelSize(0.07);
gme.GetYaxis().SetLabelSize(0.05);

gme.GetXaxis().SetTitleOffset(-0.5);
gme.GetYaxis().SetTitleOffset(0.5);
gme.GetXaxis().SetTitle("m_{4l} [GeV]");
gme.GetYaxis().SetTitle("Expected Significance");
gcomb.SetLineColor(4);
gcomb.SetLineWidth(2);
gcomb.GetXaxis().SetTitle("m_{4l} [GeV]");
gcomb.GetYaxis().SetTitle("Expected Significance");
gme.Draw("ACP");
gcomb.Draw("CP same");


legend.AddEntry(gcomb, "ICHEP datacards w/ MELA", "l");
legend.AddEntry(gme, "2D Binned likelihood w/ BNN", "l");
legend.Draw("same")

pad1.Update()
ceff.Update()

ceff.cd();
pad2.SetGridy();
pad2.SetTopMargin(0)
pad2.SetBottomMargin(4)

pad2.Draw()
pad2.cd()
#hratio = TH1F("ratio", "ratio", 10, 114.5, 142.2)
hratio = TGraph()

for i in xrange(12):
#	hratio.SetBinContent(i+1, gme.Eval(117+2*i)/gcomb.Eval(117+2*i))
	hratio.SetPoint(hratio.GetN(), 117+2*i, gme.Eval(117+2*i)/gcomb.Eval(117+2*i))
#	hratio.SetBinError(i, 0)
hratio.GetYaxis().SetLabelSize(0.09);
hratio.GetXaxis().ImportAttributes(gme.GetXaxis())
hratio.GetXaxis().SetLabelSize(0.13);
hratio.GetXaxis().SetTitleOffset(0.65);
hratio.GetXaxis().SetLimits(114.6, 142.4);
hratio.GetYaxis().SetTitleOffset(0.65);
hratio.SetMarkerStyle(20);
hratio.Draw("acp")
pad2.Update()
ceff.Update()

ceff.Update();
ceff.SaveAs("signif_comp.gif")
ceff.Update();

sleep(10)