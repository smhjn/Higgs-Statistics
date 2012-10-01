from string import *
from ROOT import *
from time import sleep
import os, sys, re
from array import *

gROOT.Reset()
#gROOT.ProcessLine(".L classes/TR130NoDY_6V_blabo.c")

#file = open('/home/jbochenek/work/HZZ4l/HiggsToZZ4Leptons/FSR2012/DATA2011/roottree_leptons_VBF_HToZZTo4L_M-600_8TeV-powheg-pythia6_e_fr.txt')
file = open('/home/jbochenek/data/HZZ4l_2012_data/fakes/Run2012_e_fr.txt')
label = "Run2012_e"

vars = [
'pt',
'eta',
'charge',
'iso',
'sip',
'massZ1',
'nloose',
'mvaid',
'istight',
]

bins = 15
plotmax = 15

density_loose = TH2F("density_loose", "density_loose", bins, -2.5, 2.5, bins, 0, 100)
density_tight = TH2F("density_tight", "density_tight", bins, -2.5, 2.5, bins, 0, 100)
fake_ratio = TH2F("fake_ratio", "fake_ratio", bins, -2.5, 2.5, bins, 0, 100)

bins1d = 9
edges = array('d',[0, 7, 10, 15, 20, 25, 30, 40, 50, 80])

density_barrel_loose = TH1F("density_barrel_loose", "density_barrel_loose", len(edges)-1, edges)
density_endcap_loose = TH1F("density_endcap_loose", "density_endcap_loose", len(edges)-1, edges)
density_barrel_tight = TH1F("density_barrel_tight", "density_barrel_tight", len(edges)-1, edges)
density_endcap_tight = TH1F("density_endcap_tight", "density_endcap_tight", len(edges)-1, edges)

fake_ratio_barrel = TH1F("fake_ratio_barrel", "fake_ratio_barrel", len(edges)-1, edges)
fake_ratio_endcap = TH1F("fake_ratio_endcap", "fake_ratio_endcap", len(edges)-1, edges)

f = TFile('fake_reweight.root', "update")
f.cd()

line = file.readline()
data = line.rsplit()
columns = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns.append(varname); 

numtight = 0
numloose = 0

for line in file:
	data = line.rsplit()
	d = dict(zip(columns, data))

	if(float(d["eta"]) < 1.5):
		density_barrel_loose.Fill(float(d["pt"]), 1.0)
	else:
		density_endcap_loose.Fill(float(d["pt"]), 1.0)

	density_loose.Fill(float(d["eta"]), float(d["pt"]), 1.0)
		
	numloose+=1
		
	notselected = 0
	if((float(d["istight"]) == 0) ):
		notselected+=1
	if((float(d["iso"]) > 0.4) ):
		print "iso problem"
		notselected+=1
	
	if notselected == 0:
		density_tight.Fill(float(d["eta"]), float(d["pt"]), 1.0)
		numtight+=1		
		if(float(d["eta"]) < 1.5):
			density_barrel_tight.Fill(float(d["pt"]), 1.0)
		else:
			density_endcap_tight.Fill(float(d["pt"]), 1.0)			
		
		
#	if(float(d["pt"]) < 5): 
#		notselected+=1
#	if(float(d["pt"]) < 10):
#		if((float(d["eta"]) < 0.8) & (float(d["mvaid"]) < 0.47)):
#			notselected+=1;
#		if((float(d["eta"]) > 0.8) & (float(d["eta"]) < 1.479) & (float(d["mvaid"]) < 0.004)):
#			notselected+=1;
#		if((float(d["eta"]) > 1.479) & (float(d["mvaid"]) < 0.295)):
#			notselected+=1;
#	else:
#		if((float(d["eta"]) < 0.8) & (float(d["mvaid"]) < 0.5)):
#			notselected+=1;
#		if((float(d["eta"]) > 0.8) & (float(d["eta"]) < 1.479) & (float(d["mvaid"]) < 0.12)):
#			notselected+=1;
#		if((float(d["eta"]) > 1.479) & (float(d["mvaid"]) < 0.6)):
#			notselected+=1;	


print str(numtight) + " " + str(numloose)
density_loose.Smooth(20)
density_tight.Smooth(20)

for binx in xrange(bins):
	for biny in xrange(bins):	
		bincontx = density_loose.GetBinContent(binx+1, biny+1)
		if(bincontx): 
			ratio = density_tight.GetBinContent(binx+1, biny+1)/density_loose.GetBinContent(binx+1, biny+1)
			fake_ratio.SetBinContent(binx+1, biny+1, ratio)



fake_ratio_endcap = density_endcap_tight.Clone("h")
density_endcap_loose.Sumw2()
fake_ratio_endcap.Sumw2()
fake_ratio_endcap.Divide(density_endcap_loose)

fake_ratio_barrel = density_barrel_tight.Clone("h")
density_barrel_loose.Sumw2()
fake_ratio_barrel.Sumw2()
fake_ratio_barrel.Divide(density_barrel_loose)


			
fake_ratio.Smooth(1)
			
ctight = TCanvas( 'ctight', 'ctight', 200, 10, 700, 500 )
ctight.cd();
density_tight.SetFillColor(3)
density_tight.SetMarkerColor(3)
density_tight.SetLineColor(3)
#density_tight.SetMaximum(plotmax);

density_tight.Draw("colz");
ctight.Update();
ctight.SaveAs("ctight.gif")
ctight.Update();


cloose = TCanvas( 'cloose', 'cloose', 200, 10, 700, 500 )
cloose.cd();
density_loose.SetFillColor(2)
density_loose.SetMarkerColor(2)
density_loose.SetLineColor(2)
#density_loose.SetMaximum(plotmax);
density_loose.Draw("colz");
cloose.Update();
cloose.SaveAs("cloose.gif")
cloose.Update();

cratio = TCanvas( 'cratio', 'cratio', 200, 10, 700, 500 )
cratio.cd();
fake_ratio.SetName("cratio_" + label)
fake_ratio.SetTitle("Fake ratio " + label)
gStyle.SetOptStat(0);
fake_ratio.GetXaxis().SetTitle("#eta");
fake_ratio.GetYaxis().SetTitle("p_{T}");
fake_ratio.SetFillColor(1)
fake_ratio.SetMarkerColor(1)
fake_ratio.SetLineColor(1)
#fake_ratio.SetMaximum(1.5);
fake_ratio.Draw("colz");
cratio.Update();
cratio.SaveAs("cratio_" + label + ".gif")
cratio.Update();
fake_ratio.Write();

cratio_pt = TCanvas( 'cratio_pt', 'cratio_pt', 200, 10, 700, 500 )
cratio_pt.cd();
fake_ratio.ProjectionY().Draw("HIST");
cratio_pt.Update();
cratio_pt.SaveAs("cratio_pt_" + label + ".gif")
cratio_pt.Update();

c = TCanvas( 'cratio_1d', 'cratio_1d', 200, 10, 700, 500 )
c.cd();
legend = TLegend(0.3, 0.72, 0.68, 0.88)
fake_ratio_barrel.SetMinimum(0.)
fake_ratio_barrel.SetTitle("Fake Rates, Electrons 8TeV")
fake_ratio_barrel.GetYaxis().SetTitle("Fake Ratio")
fake_ratio_barrel.GetXaxis().SetTitle("p_{T}")
fake_ratio_barrel.SetLineColor(4)
fake_ratio_barrel.Draw("LEP");
fake_ratio_endcap.SetLineColor(2)
fake_ratio_endcap.Draw("LEP same");
legend.AddEntry(fake_ratio_barrel, "barrel", "l");
legend.AddEntry(fake_ratio_endcap, "endcap", "l");
legend.Draw("same")
c.Update();
c.SaveAs("cratio_1d_" + label + ".gif")
c.Update();


f.Write()
f.Close()

sleep(10)