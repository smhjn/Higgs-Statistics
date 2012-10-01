from string import *
from ROOT import *
from time import sleep
import os, sys, re, bnn_s_vs_ps_7vars_ext, bnn_2e2mu_loosesel4_ext, bnn_4e_loosesel4_ext, bnn_4mu_loosesel4_ext, applycuts_ext
from applycuts_ext import applycuts

gROOT.Reset()


class Complex:
	def __init__(self, realpart, imagpart):
		self.r = realpart
		self.i = imagpart

def cdf(hist):
	n = hist.GetNbinsX()
	c = []
	print n
	sum = 0
	for i in range(n):
		sum = sum + hist.GetBinContent(i)
		c.append(sum)
	return c;

vdouble = vector("double")
file1 = open('../vars_pvps_7vars.txt')
mvavars = []
for line in file1:
	mvavars.append(line.strip())


#sig = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_ps_train.dat')
#bkg = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_s_train.dat')
sig2 = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/sig_bnn_train_ps_uw.dat')
bkg2 = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/sig_bnn_train_s_uw.dat')
sig = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_s_test.dat')
bkg = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_ps_test.dat')
zz = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/bkg_test.dat')
zz2 = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/bkg_mela.dat')
zx = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/zx_test.dat')
zx2 = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/zx_mela.dat')

file1 = open('../mvavars1.txt')
file2 = open('../mvavars2.txt')
file3 = open('../mvavars3.txt')

mvavars1 = []
mvavars2 = []
mvavars3 = []

for line in file1:
	mvavars1.append(line.strip())
for line in file2:
	mvavars2.append(line.strip())
for line in file3:
	print line
	mvavars3.append(line.strip())


line = sig.readline()
data = line.rsplit()
columns = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns.append(varname); 

bins = 200
heff = TH2F('heff_train', 'heff_train', bins, 0, 1.0,bins, 0.0, 1)
heff_m = TH2F('heff_train_m', 'heff_train_m', bins, 0, 1.0,bins, 0.0, 1)

heffmela = TH2F('heff_mela', 'heff_mela', bins, 0, 1.0,bins, 0.0, 1)


hs = TH1F('mva_s', 'mva_s', bins, 0, 1)
hs.SetLineColor(2);
hb = TH1F('mva_b', 'mva_b', bins, 0, 1)
hb.SetLineColor(1);
hs.Draw()
hb.Draw("same")


hs_m = TH1F('mva_s_m', 'mva_s_m', bins, 0, 1)
hs_m.SetLineColor(2);
hb_m = TH1F('mva_b_m', 'mva_b_m', bins, 0, 1)
hb_m.SetLineColor(1);
hs_m.Draw()
hb_m.Draw("same")

i = 0

c1 = TCanvas( 'c1', 'mva output', 200, 10, 700, 500 )
c1.SetGridx()
c1.SetGridy()
c1.Draw()
gStyle.SetOptStat(0);

c5 = TCanvas( 'c5', 'mela output', 200, 10, 700, 500 )
c5.SetGridx()
c5.SetGridy()
c5.Draw()
gStyle.SetOptStat(0);

c2 = TCanvas( 'mass', 'mass', 200, 10, 700, 500 )
c2.cd()
c2.SetGridx()
c2.SetGridy()

c3 = TCanvas( 'z1mass', 'z1mass', 200, 10, 700, 500 )
c3.cd()
c3.SetGridx()
c3.SetGridy()

c4 = TCanvas( 'z2mass', 'z2mass', 200, 10, 700, 500 )
c4.cd()
c4.SetGridx()
c4.SetGridy()

hsmass = TH1F('mva_smass', 'mva_smass', bins, 110, 150)
hsmass.SetLineColor(2);
hbmass = TH1F('mva_bmass', 'mva_bmass', bins, 110, 150)
hbmass.SetLineColor(1);

hsz1 = TH1F('hsz1', 'hsz1', bins, 30, 100)
hsz1.SetLineColor(1);
hbz1 = TH1F('hbz1', 'hbz1', bins, 30, 100)
hbz1.SetLineColor(1);

hsz2 = TH1F('hsz1', 'hsz1', bins, 0, 100)
hsz2.SetLineColor(1);
hbz2 = TH1F('hbz1', 'hbz1', bins, 0, 100)
hbz2.SetLineColor(1);

c1.cd()

line = sig2.readline()	
data = line.rsplit()
columns2 = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns2.append(varname); 

bins = 20
htemps = TH2F('htemps', 'htemps', bins, 0, 1.0,bins, 0.0, 1)
htempps = TH2F('htempps', 'htempps', bins, 0, 1.0,bins, 0.0, 1 )
htempzz = TH2F('htempzz', 'htempzz', bins, 0, 1.0,bins, 0.0, 1)
htempzx = TH2F('htempzx', 'htempzx', bins, 0, 1.0,bins, 0.0, 1)

for line in sig:
	
	data = line.rsplit()
	data_ = map(float, line.rsplit())
	d = dict(zip(columns, data))
	inputs = []
	channel_ = int(d["channel"])
		
	for var in mvavars:
		inputs.append(float(d[var]))

	y = bnn_s_vs_ps_7vars_ext.bnn_s_vs_ps_7vars(inputs, 100, 200);
	hs.Fill(y)
	hs_m.Fill(float(d["pmela"]))
	hsmass.Fill(float(d["mass4l"]))
	hsz1.Fill(float(d["Z1mass"]))
	hsz2.Fill(float(d["Z2mass"]))

	bnnline = sig2.readline()

	data = bnnline.rsplit()
	data_ = map(float, bnnline.rsplit())

	d = dict(zip(columns2, data))
	channel_ = d["channel"]

	y2 = 0
	inputs = []
	if(int(d["channel"]) == 1):
		for var in mvavars1:
			inputs.append(float(d[var]))
		y2 = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
	if(int(d["channel"]) == 2):
		for var in mvavars2:
			inputs.append(float(d[var]))
		y2 = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
	if(int(d["channel"]) == 3):
		for var in mvavars3:
			inputs.append(float(d[var]))
		y2 = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);

	if(applycuts(columns2, data_, int(d["channel"]))):
		btotal_sel += (weight)
	else:
		htempps.Fill(y, y2)
	
	i+=1
	if( not (i%10000)):
		print i
		print y	
		hs.Draw()
		c1.Update()
sig.close()

print "Integral hs: " + str(hs.Integral())
hs.Scale(1/hs.Integral())
hs_m.Scale(1/hs_m.Integral())

hsmass.Scale(1/hsmass.Integral())
hsz1.Scale(1/hsz1.Integral())
hsz2.Scale(1/hsz2.Integral())

#do background
line = bkg.readline()
average = 0
btotal_ = 0
btotal_sel = 0

i = 0

line = bkg2.readline()	
data = line.rsplit()
columns2 = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns2.append(varname); 
	
for line in bkg:

	data = line.rsplit()
	data_ = map(float, line.rsplit())
	d = dict(zip(columns, data))
	channel_ = d["channel"]

	inputs = []
	for var in mvavars:
		inputs.append(float(d[var]))
	y = bnn_s_vs_ps_7vars_ext.bnn_s_vs_ps_7vars(inputs, 100, 200);
	
	hb.Fill(y)
	hb_m.Fill(float(d["pmela"]))
	
	hbmass.Fill(float(d["mass4l"]))
	hbz1.Fill(float(d["Z1mass"]))
	hbz2.Fill(float(d["Z2mass"]))

	bnnline = bkg2.readline()

	data = bnnline.rsplit()
	data_ = map(float, bnnline.rsplit())
	d = dict(zip(columns2, data))
	channel_ = d["channel"]

	y2 = 0
	inputs = []
	if(int(d["channel"]) == 1):
		for var in mvavars1:
			inputs.append(float(d[var]))
		y2 = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
		
	if(int(d["channel"]) == 2):
		for var in mvavars2:
			inputs.append(float(d[var]))
		y2 = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
	if(int(d["channel"]) == 3):
		for var in mvavars3:
			inputs.append(float(d[var]))
		y2 = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);

	if(applycuts(columns2, data_, int(d["channel"]))):
		btotal_sel += (weight)
	else:
		htemps.Fill(y, y2)
	i+=1
	if( not (i%10000)):
		hb.Draw("same")

		print i
		print y	
		c1.Update();
	average += float(d["mass4l"])

print "btotal_: " + str(btotal_)

print "Integral hb: " + str(hb.Integral())

line = zz.readline()	
line = zz2.readline()	

data = line.rsplit()
columns2 = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns2.append(varname); 
	
for line in zz:
	data = line.rsplit()
	data_ = map(float, line.rsplit())
	d = dict(zip(columns, data))
	channel_ = d["channel"]

	inputs = []
	for var in mvavars:
		inputs.append(float(d[var]))
	y = bnn_s_vs_ps_7vars_ext.bnn_s_vs_ps_7vars(inputs, 100, 200);
	
	hb.Fill(y)
	hb_m.Fill(float(d["pmela"]))
	
	hbmass.Fill(float(d["mass4l"]))
	hbz1.Fill(float(d["Z1mass"]))
	hbz2.Fill(float(d["Z2mass"]))

	bnnline = zz2.readline()

	data = bnnline.rsplit()
	data_ = map(float, bnnline.rsplit())
	d = dict(zip(columns2, data))
	channel_ = d["channel"]

	y2 = 0
	inputs = []
	if(int(d["channel"]) == 1):
		for var in mvavars1:
			inputs.append(float(d[var]))
		y2 = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
	if(int(d["channel"]) == 2):
		for var in mvavars2:
			inputs.append(float(d[var]))
		y2 = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
	if(int(d["channel"]) == 3):
		for var in mvavars3:
			inputs.append(float(d[var]))
		y2 = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);

	weight = float(d["weight"])

	if(float(d["mass4l"]) < 122):
		continue
	if(float(d["mass4l"]) > 128):
		continue

	if(applycuts(columns2, data_, int(d["channel"]))):
#			print "notselected"		
		btotal_sel += (weight)
#			continue
	else:
		htempzz.Fill(y, y2, weight)
	i+=1
	if( not (i%10000)):
		hb.Draw("same")

		print i
		print y	
		c1.Update();
	average += float(d["mass4l"])

line = zx.readline()	
line = zx2.readline()	

data = line.rsplit()
columns2 = []
for varname in data:
	print varname
	varname = varname.lstrip(':b')
	columns2.append(varname); 
	
f_fake = TFile('../fake_reweight.root')
f_fake.cd();
cratio_Run2012_e=gDirectory.Get("cratio_Run2012_e");
cratio_Run2012_mu= gDirectory.Get("cratio_Run2012_mu");

	
for line in zx:
	data = line.rsplit()
	data_ = map(float, line.rsplit())
	d = dict(zip(columns, data))
	channel_ = d["channel"]

	inputs = []
	for var in mvavars:
		inputs.append(float(d[var]))
	y = bnn_s_vs_ps_7vars_ext.bnn_s_vs_ps_7vars(inputs, 100, 200);
	
	hb.Fill(y)
	hb_m.Fill(float(d["pmela"]))
	
	hbmass.Fill(float(d["mass4l"]))
	hbz1.Fill(float(d["Z1mass"]))
	hbz2.Fill(float(d["Z2mass"]))

	bnnline = zx2.readline()

	data = bnnline.rsplit()
	data_ = map(float, bnnline.rsplit())
	d = dict(zip(columns2, data))
	channel_ = d["channel"]

	y2 = 0
	reweight = 1
	inputs = []
	if(int(d["channel"]) == 1):
		for var in mvavars1:
			inputs.append(float(d[var]))
		y2 = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
		reweight *= cratio_Run2012_mu.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_mu.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));

	if(int(d["channel"]) == 2):
		for var in mvavars2:
			inputs.append(float(d[var]))
		y2 = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
		reweight *= cratio_Run2012_e.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_e.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));

	if(int(d["channel"]) == 3):
		for var in mvavars3:
			inputs.append(float(d[var]))
		y2 = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);
		reweight *= cratio_Run2012_e.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_e.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));


	d["lept3_pfx"] = "0"
	d["lept4_pfx"] = "0"
	d["lept3_mvaid"] = "1.0"
	d["lept4_mvaid"] = "1.0"

	if(float(d["mass4l"]) < 122):
		continue
	if(float(d["mass4l"]) > 128):
		continue

	for k, v in enumerate(data_):
		data_[k] = float(d[columns2[k]])

	if(applycuts(columns2, data_, int(d["channel"]))):
#			print "notselected"		
		btotal_sel += weight
	else:
		htempzx.Fill(y, y2, reweight)

	i+=1
	if( not (i%10000)):
		hb.Draw("same")

		print i
		print y	
		c1.Update();
	average += float(d["mass4l"])

hb.Scale(1/hb.Integral())
hb_m.Scale(1/hb_m.Integral())

hbmass.Scale(1/hbmass.Integral())
hbz1.Scale(1/hbz1.Integral())
hbz2.Scale(1/hbz2.Integral())

f = TFile('models.root', "recreate")
f.cd()
htemps.Write()
htempps.Write()
htempzz.Write()
htempzx.Write()
f.Close()




c2.cd()
hsmass.SetLineColor(2)
hbmass.SetLineColor(1)
hbmass.Draw()
hsmass.Draw("same")
c2.Update()
c2.SaveAs("hmass_pvsps.gif")

c3.cd()
hsz1.SetLineColor(2)
hbz1.SetLineColor(1)
hbz1.Draw()
hsz1.Draw("same")
c3.Update()
c3.SaveAs("hz1_pvsps.gif")

c4.cd()
hsz2.SetLineColor(2)
hbz2.SetLineColor(1)
hbz2.Draw()
hsz2.Draw("same")
c4.Update()
c4.SaveAs("hz2_pvsps.gif")


# BNN eff

stotal=hs.Integral();
btotal=hb.Integral();
cs = cdf(hs);
cb = cdf(hb);

maxsig = 0;
optmvax = 0;
optmvay = 0;

geff = TGraph()
geff.SetName("bnneff")

for i in range(hs.GetNbinsX()):
	cs1 = stotal - cs[i-1]
	cb1 = btotal - cb[i-1]
	es = cs1 / stotal
	eb = cb1 / btotal
	sign=0
	if(cb1 > 0): 
		sign = cs1/sqrt(cs1 + cb1)
	if(sign > maxsig):
		maxsig = sign;
		optcut = i;	 
		optmvax = eb;
		optmvay = es;
	heff.Fill(eb, es);
	if(i<2):
		continue
	else:
		geff.SetPoint(geff.GetN(), eb, es);

print "\n Maximum Significance " + str(maxsig) + ", is bin: " + str(optcut) + ", which is BNN Cut: " + str(optcut/bins)
print "es: " + str(optmvay) + "\teb: " + str(optmvax)


geff_m = TGraph()
geff_m.SetName("bnneff_m")

stotal=hs_m.Integral();
btotal=hb_m.Integral();
cs = cdf(hs_m);
cb = cdf(hb_m);

maxsig = 0;
optmvax = 0;
optmvay = 0;

for i in range(hs_m.GetNbinsX()):
	cs1 = stotal - cs[i-1]
	cb1 = btotal - cb[i-1]
	es = cs1 / stotal
	eb = cb1 / btotal
	sign=0
	if(cb1 > 0): 
		sign = cs1/sqrt(cs1 + cb1)
	if(sign > maxsig):
		maxsig = sign;
		optcut = i;	 
		optmvax = eb;
		optmvay = es;
	heff_m.Fill(es, eb);
	if(i<2):
		continue
	else:
		geff_m.SetPoint(geff_m.GetN(), es, eb);

print "\n Maximum Significance " + str(maxsig) + ", is bin: " + str(optcut) + ", which is BNN Cut: " + str(optcut/bins)
print "es: " + str(optmvay) + "\teb: " + str(optmvax)


# MELA eff


stotal=hs_m.Integral();
btotal=hb_m.Integral();
cs = cdf(hs_m);
cb = cdf(hb_m);

maxsig = 0;
optmvax = 0;
optmvay = 0;

print "\n Maximum Significance " + str(maxsig) + ", is bin: " + str(optcut) + ", which is BNN Cut: " + str(optcut/bins)
print "es: " + str(optmvay) + "\teb: " + str(optmvax)



gline = TGraph()
for i in range(11):
	gline.SetPoint(gline.GetN(), i*0.1, i*0.1)

ceff = TCanvas( 'ceff', 'eff', 200, 10, 500, 500 )
ceff.cd();
legend = TLegend(0.12, 0.70, 0.45, 0.9)

gStyle.SetOptStat(0);
heff.SetTitle("ROC: BNN & pMELA");
heff.GetXaxis().SetTitle("#epsilon_{0^{+}}");
heff.GetXaxis().SetTitleSize(0.06);
heff.GetYaxis().SetTitleSize(0.06);
heff.GetXaxis().SetTitleOffset(0.65);
heff.GetYaxis().SetTitleOffset(0.65);
heff.GetYaxis().SetTitle("#epsilon_{0^{-}}");

geff.SetLineWidth(2);
geff_m.SetLineWidth(2);

heff.SetLineColor(3)
geff.SetLineColor(3)
heff.Draw()
geff.Draw("same");

heff_m.SetLineColor(2)
heff_m.Draw("same")
geff_m.SetLineColor(2)
geff_m.Draw("same");

gline.Draw("same");
legend.AddEntry(geff, "BNN", "l");
legend.AddEntry(geff_m, "pseudoMELA", "l");
legend.Draw("same")

ceff.Update();
ceff.SaveAs("heff_psvss.gif")
ceff.Update();


c1.cd()
legend = TLegend(0.8, 0.65, 0.95, 0.85)
gStyle.SetOptStat(0);
hb.Rebin(8)
hs.Rebin(8)
hb.SetTitle("BNN Discriminant (SM Higgs vs. Pseudo-Scalar Higgs)");
hb.GetXaxis().SetTitle("D(x)");
hb.GetXaxis().SetTitleSize(0.06);
hb.GetYaxis().SetTitleSize(0.06);
hb.GetXaxis().SetTitleOffset(0.65);
hb.GetYaxis().SetTitleOffset(0.65);
hs.Draw()
hb.Draw("same")
hb.SetLineWidth(2);
hs.SetLineWidth(2);
legend.AddEntry(hs, "0^{+}", "l");
legend.AddEntry(hb, "0^{-}", "l");
legend.Draw("same")
c1.Update()
c1.SaveAs("hdist_bnn_pvsps.gif")

c5.cd()
legend = TLegend(0.8, 0.65, 0.95, 0.85)
gStyle.SetOptStat(0);
hb_m.Rebin(8)
hs_m.Rebin(8)

hb_m.SetTitle("MELA Discriminant (SM Higgs vs. Pseudo-Scalar Higgs)");
hb_m.GetXaxis().SetTitle("D(x)");
hb_m.GetXaxis().SetTitleSize(0.06);
hb_m.GetYaxis().SetTitleSize(0.06);
hb_m.GetXaxis().SetTitleOffset(0.65);
hb_m.GetYaxis().SetTitleOffset(0.65);
hs_m.Draw()
hb_m.Draw("same")
hb_m.SetLineWidth(2);
hs_m.SetLineWidth(2);
legend.AddEntry(hs_m, "0^{+}", "l");
legend.AddEntry(hb_m, "0^{-}", "l");
legend.Draw("same")
c5.Update()
c5.SaveAs("hdist_mela_pvsps.gif")

