from string import *
from ROOT import *
from time import sleep
import os, sys, re
from array import *

gROOT.Reset()

datafile = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/data.dat')
mcfile = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg.dat')

bins = 60
bkgs = 5
thischannel = 3

name = ""
if (thischannel == 2):
	name = "4e"
if (thischannel == 1):
	name = "4#mu"	
if (thischannel == 3):
	name = "2e2#mu"	
	

samesign_data = TH1F("samesign_data", "samesign_data", bins, 100, 600)
oppsign_data = TH1F("oppsign_data", "oppsign_data", bins, 100, 600)

samesign_mc = []
oppsign_mc = []

for k in xrange(bkgs):
	print k
	samesign_mc.append(TH1F("samesign_mc_" + str(k), "samesign_mc_" + str(k), bins, 100, 600))
	oppsign_mc.append(TH1F("oppsign_mc_" + str(k), "oppsign_mc_" + str(k), bins, 100, 600))

#do data

#get header
line = datafile.readline()
data = line.rsplit()
columns = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns.append(varname); 


#loop through file
for line in datafile:
	data = line.rsplit()
	d = dict(zip(columns, data))
	channel = int(d["channel"])
	if(channel == thischannel):
		if(int(d["issamesign"]) > 0):
			samesign_data.Fill(float(d["mass4l"]))
		else:
			oppsign_data.Fill(float(d["mass4l"]))


#get header
line = mcfile.readline()
data = line.rsplit()
columns = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns.append(varname); 

i = 0
#loop through file
for line in mcfile:
	data = line.rsplit()
	d = dict(zip(columns, data))
	channel = int(d["channel"])
	sample = int(d["sample"])-1

	if(channel == thischannel):
		if int(d["issamesign"]):
			samesign_mc[sample].Fill(float(d["mass4l"]), float(d["weight"]))
		else:
			oppsign_mc[sample].Fill(float(d["mass4l"]), float(d["weight"]))
		i+=1
		if( not (i%1000)):
			print i
			print sample


bkgnames = [
"Z+jets",
"Top",
"WZ",
"ZZ",
"WW"
]

bkgcolors = [
8,
16,
9,
38,
46
]

bkgscales = [
1.,
1.,
1.,
1.,
1.
]

c = TCanvas( 'cratio_1d', 'cratio_1d', 200, 10, 700, 500 )
c.cd();
legend = TLegend(0.6, 0.62, 0.98, 0.88)
htotal = THStack("stack","stack");
gStyle.SetEndErrorSize(3);
gStyle.SetErrorX(1.);
gStyle.SetOptStat(0);

samesign_data.SetMinimum(0.)
samesign_data.SetTitle("SS-SF Control Region, " + name)
samesign_data.GetYaxis().SetTitle("Events / 10 GeV")
samesign_data.GetXaxis().SetTitle("m_{4l}")
samesign_data.SetLineColor(1)
samesign_data.SetMarkerStyle(20);
samesign_data.Draw("E1");
legend.AddEntry(samesign_data, "barrel", "lep");


for k in xrange(bkgs):
	print k
	print bkgscales[k]
	print bkgnames[k]
	print "integral: " + str(samesign_mc[k].Integral())
	samesign_mc[k].Scale(bkgscales[k])
	print "integral (rw): " + str(samesign_mc[k].Integral())
	htotal.Add(samesign_mc[k])
	samesign_mc[k].SetFillColor(bkgcolors[k])
	legend.AddEntry(samesign_mc[k], bkgnames[k], "f");
htotal.Draw("HIST same");
samesign_data.Draw("E1 same");

legend.Draw("same")
c.Update();
c.SaveAs("m_4l_samesign_" + str(thischannel) + ".gif")
c.Update()

c = TCanvas( 'cratio_1d', 'cratio_1d', 200, 10, 700, 500 )
c.cd();
legend = TLegend(0.6, 0.62, 0.98, 0.88)
htotal = THStack("stack","stack");

oppsign_data.SetMinimum(0.)
oppsign_data.SetTitle("OS-SF Control Region, " + name)
oppsign_data.GetYaxis().SetTitle("Events / 10 GeV")
oppsign_data.GetXaxis().SetTitle("m_{4l}")
oppsign_data.SetLineColor(1)
oppsign_data.SetMarkerStyle(20);

oppsign_data.Draw("E1");
legend.AddEntry(oppsign_data, "data", "lep");
print "data integral (os): " + str(oppsign_data.Integral())
print "data integral (ss): " + str(samesign_data.Integral())


for k in xrange(bkgs):
	print bkgnames[k] + " integral (os): " + str(oppsign_mc[k].Integral())
	oppsign_mc[k].Scale(bkgscales[k])
	print bkgnames[k] + " integral (ss): " + str(samesign_mc[k].Integral())
	htotal.Add(oppsign_mc[k])
	oppsign_mc[k].SetFillColor(bkgcolors[k])
	legend.AddEntry(oppsign_mc[k], bkgnames[k], "f");
htotal.Draw("HIST same");
oppsign_data.Draw("E1 same");

legend.Draw("same")
c.Update();
c.SaveAs("m_4l_oppsign_" + str(thischannel) + ".gif")
c.Update()
