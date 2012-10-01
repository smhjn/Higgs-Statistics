from string import *
from ROOT import *
from time import sleep
import os, sys, re, bnn_2e2mu_loosesel4_ext, bnn_4e_loosesel4_ext, bnn_4mu_loosesel4_ext, applycuts_ext
from applycuts_ext import applycuts

gROOT.Reset()
#gROOT.ProcessLine('.L ~/work/MELA/scripts/MELA.C')

#gROOT.ProcessLine(".L classes/TR130NoDY_6V_blabo.c")

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

#weight	weight	ele1_pt	ele1_eta	ele1_phi	ele1_charge	ele1_trackIso	ele1_EcalIso	ele1_HcalIso	ele1_X	ele1_SIP	ele2_pt	ele2_eta	ele2_phi	ele2_charge	ele2_trackIso	ele2_EcalIso	ele2_HcalIso	ele2_X	ele2_SIP	ele3_pt	ele3_eta	ele3_phi	ele3_charge	ele3_trackIso	ele3_EcalIso	ele3_HcalIso	ele3_X	ele3_SIP	ele4_pt	ele4_eta	ele4_phi	ele4_charge	ele4_trackIso	ele4_EcalIso	ele4_HcalIso	ele4_X	ele4_SIP	

#worst_iso_X	second_worst_iso_X	worst_vertex	second_worst_vertex	mZ	mZstar	mbestH	index	channel	sample	end

def eventloop(sig,bkg,channel,evtset):
	line = sig.readline()
	data = line.rsplit()
	columns = []
	for varname in data:
		print varname
		varname = varname.lstrip(':b')
		columns.append(varname); 

	bins = 200
	heff = TH2F('heff_train', 'heff_train', bins, 0, 1.0,bins, 0.0, 1)
	heffmela = TH2F('heff_mela', 'heff_mela', bins, 0, 1.0,bins, 0.0, 1)


	hbnnvsmela_zz = TH2F('hbnnvsmela', 'hbnnvsmela', 20, 0, 1, 20, 0, 1)
	hbnnvsmela_sig = TH2F('hbnnvsmela', 'hbnnvsmela', 20, 0, 1, 20, 0, 1)
	hbnnvsmela_zx = TH2F('hbnnvsmela', 'hbnnvsmela', 20, 0, 1, 20, 0, 1)
	
	hs = TH1F('mva_s', 'mva_s', bins, 0, 1)
	hs.SetLineColor(2);
	hb = TH1F('mva_b', 'mva_b', bins, 0, 1)
	hb.SetLineColor(1);
	hs.Draw()
	hb.Draw("same")
	
	hb_zx = TH1F('hb_zx', 'hb_zx', 20, 0, 1)
	hb_zx.SetLineColor(2);
	hb_zx.Draw()

	hb_zx_rw = TH1F('hb_zx', 'hb_zx', 20, 0, 1)
	hb_zx_rw.SetLineColor(2);
	hb_zx_rw.Draw()

	hb_zx_sel = TH1F('hb_zx', 'hb_zx', 20, 0, 1)
	hb_zx_sel.SetLineColor(2);
	hb_zx_sel.Draw()

	
	hmelas = TH1F('mela_s', 'mela_s', bins, 0, 1)
	hmelas.SetLineColor(2);
	hmelab = TH1F('mela_b', 'mela_b', bins, 0, 1)
	hmelab.SetLineColor(1);
	hmelas.Draw()
	hmelab.Draw("same")
	
	
	i = 0
	
	c1 = TCanvas( 'c1', 'mva output', 200, 10, 700, 500 )
	c1.cd()
	c1.SetGridx()
	c1.SetGridy()
	c1.Draw()
	gStyle.SetOptStat(0);




	#------- sig events ---------

	stotal_ = 0
	stotal_sel = 0
	sigweight = 0
	
	#do signal
	for line in sig:
		data = line.rsplit()
#		data_ = [float(x) for x in line.rsplit()]
		data_ = map(float, line.rsplit())
		d = dict(zip(columns, data))
		inputs = []

		if(float(d["mass4l"]) < 110):
			continue
		if(float(d["mass4l"]) > 150):
			continue
		
		channel_ = int(d["channel"])
		weight = float(d["weight"])
		if(int(d["channel"]) == 1):
			for var in mvavars1:
				inputs.append(float(d[var]))
			y = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
		if(int(d["channel"]) == 2):
			for var in mvavars2:
				inputs.append(float(d[var]))
			y = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
		if(int(d["channel"]) == 3):
			for var in mvavars3:
				inputs.append(float(d[var]))
			y = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);

		weight = float(d["weight"])
		stotal_ += float(d["weight"])
		if(applycuts(columns, data_, channel_)):
			stotal_sel += float(d["weight"])
#			continue
		else:
			hmelas.Fill(float(d["mela"]), weight)
			hs.Fill(y, weight)
			hbnnvsmela_sig.Fill(y, float(d["mela"]), weight)
	
		if(int(d["index"]) == 6):
			sigweight+=weight

		i+=1
		if( not (i%10000)):
			print i
			print y	
			hs.Draw()
			c1.Update()
	sig.close()
	
	print "Higgs 125 weight: " + str(sigweight)

	print "Integral hs: " + str(hs.Integral())
	hs.Scale(1/hs.Integral())
	hmelas.Scale(1/hmelas.Integral())

	#do background
	line = bkg.readline()
	data = line.rsplit()
	columns = []
	for varname in data:
		print varname
		varname = varname.lstrip(':b')
		columns.append(varname); 

	#------- ZZ events ---------

	average = 0
	btotal_ = 0
	btotal_sel = 0


	i = 0
		
	for line in bkg:
		data = line.rsplit()
#		data_ = [float(x) for x in line.rsplit()]
		data_ = map(float, line.rsplit())
		d = dict(zip(columns, data))
		channel_ = d["channel"]

		if(float(d["mass4l"]) < 110):
			continue
		if(float(d["mass4l"]) > 150):
			continue
		
		weight = float(d["weight"])

		y = 0
		inputs = []
		if(int(d["channel"]) == 1):
			for var in mvavars1:
				inputs.append(float(d[var]))
			y = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
		if(int(d["channel"]) == 2):
			for var in mvavars2:
				inputs.append(float(d[var]))
			y = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
		if(int(d["channel"]) == 3):
			for var in mvavars3:
				inputs.append(float(d[var]))
			y = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);
		if(d["mela"] == "-nan"):
			d["mela"]= "0"
		btotal_ += float(weight)
		
		if(applycuts(columns, data_, int(d["channel"]))):
#			print "notselected"		
			btotal_sel += (weight)
#			continue
		else:
			hmelab.Fill(float(d["mela"]), weight)
			hb.Fill(y, weight)
			hbnnvsmela_zz.Fill(y, float(d["mela"]), weight)

		
		i+=1
		if( not (i%10000)):
			print i
			print y	
			c1.Update();
		average += float(d["mass4l"])
#	print "average: " + str(average/i)
	print "btotal_: " + str(btotal_)

	print "Integral hb: " + str(hb.Integral())

	dat = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/zx_mela.dat')

	line = dat.readline()
	data = line.rsplit()
	columns = []
	for varname in data:
		print varname
		varname = varname.lstrip(':b')
		columns.append(varname); 


	#------- Z+X events ---------


	f_fake = TFile('../fake_reweight.root')
	f_fake.cd();
	cratio_Run2012_e=gDirectory.Get("cratio_Run2012_e");
	cratio_Run2012_mu= gDirectory.Get("cratio_Run2012_mu");

	for line in dat:

		data = line.rsplit()
		data_ = map(float, line.rsplit())
		d = dict(zip(columns, data))

		reweight = 1
		inputs = []
		inputs = []
		if(int(d["channel"]) == 1):
			for var in mvavars1:
				inputs.append(float(d[var]))
			y = bnn_4mu_loosesel4_ext.bnn_4mu_loosesel4(inputs, 100, 200);
			reweight *= cratio_Run2012_mu.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_mu.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));

		if(int(d["channel"]) == 2):
			for var in mvavars2:
				inputs.append(float(d[var]))
			y = bnn_4e_loosesel4_ext.bnn_4e_loosesel4(inputs, 100, 200);
			reweight *= cratio_Run2012_e.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_e.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));

		if(int(d["channel"]) == 3):
			for var in mvavars3:
				inputs.append(float(d[var]))
			y = bnn_2e2mu_loosesel4_ext.bnn_2e2mu_loosesel4(inputs, 100, 200);
			reweight *= cratio_Run2012_e.Interpolate(float(d["lept3_eta"]), float(d["lept3_pt"])) * cratio_Run2012_e.Interpolate(float(d["lept4_eta"]), float(d["lept4_pt"]));

		weight = float(d["weight"])

		d["issamesign"] = "0"
		if(d["mela"] == "-nan"):
			d["mela"] = "0"
		btotal_ += float(d["weight"])

		
		d["lept3_pfx"] = "0"
		d["lept4_pfx"] = "0"
		d["lept3_mvaid"] = "1.0"
		d["lept4_mvaid"] = "1.0"


		for k, v in enumerate(data_):
			data_[k] = float(d[columns[k]])
#			print columns[k] + " " + str(v)

		if(applycuts(columns, data_, int(d["channel"]))):
#			print "notselected"		
			btotal_sel += weight
		else:
#			hb_zx_sel.Fill(float(d["mela"]), 1.)
			hmelab.Fill(float(d["mela"]), reweight)
			hbnnvsmela_zx.Fill(y, float(d["mela"]), weight)
		hb.Fill(y, reweight)

		hb_zx_rw.Fill(float(d["mela"]), reweight)
		hb_zx.Fill( float(d["mela"]), 1.)


		i+=1
		if( not (i%10000)):
			print i
			print y	
			hb.Draw("same")
			c1.Update();
		average += float(d["mass4l"])

	hmelas.Fill(0.0, stotal_sel)
	hmelab.Fill(0.0, btotal_sel)

#	hs.Fill(0.0, stotal_sel)
#	hb.Fill(0.0, btotal_sel)

#	print "average: " + str(average/i)
	print "btotal_: " + str(btotal_)
	print "Integral hb: " + str(hb.Integral())

	print "Integral zx: " + str(hb_zx.Integral())
	print "Integral zx (sel.): " + str(hb_zx_sel.Integral())
	print "Integral zx (rw): " + str(hb_zx_rw.Integral())



	#------- Print Z+X plots ---------


	c_zx = TCanvas( 'c_zx', 'mva zx output', 200, 10, 700, 500 )
	c_zx.cd()
	c_zx.SetGridx()
	c_zx.SetGridy()
	c_zx.Draw()
	gStyle.SetOptStat(0);

	hb_zx.Scale(1/hb_zx.Integral())
#	hb_zx_sel.Scale(1/hb_zx_sel.Integral())
	hb_zx_rw.Scale(1/hb_zx_rw.Integral())
	c_zx.cd()
#	c_zx.SetLogy()
	hb_zx.SetTitle("BNN Discriminant: Z+X");
	hb_zx.GetXaxis().SetTitle("D(x)");
	legend = TLegend(0.3, 0.72, 0.68, 0.88)

	hb_zx.SetLineColor(3)
	hb_zx.Draw()
	hb_zx_rw.SetLineColor(2)
	hb_zx_rw.SetLineWidth(2)
	hb_zx.SetLineWidth(2)
	hb_zx_sel.SetLineWidth(2)

	hb_zx_rw.Draw("same")
	hb_zx_sel.SetLineColor(1)
#	hb_zx_sel.Draw("same")
	hb_zx.Draw("hist same")
	legend.AddEntry(hb_zx_rw, "reweight", "l");
	legend.AddEntry(hb_zx, "unweighted", "l");
#	legend.AddEntry(hb_zx_sel, "selected", "l");
	legend.Draw("same")

	c_zx.Update()
	c_zx.SaveAs("zx-reweight_" + evtset + ".gif")



	#------- Print 2D plots ---------

	c_2d = TCanvas( 'c_2d', 'mva zx output', 200, 10, 500, 500 )
	c_2d.cd()
	c_2d.SetGridx()
	c_2d.SetGridy()
	c_2d.Draw()
	gStyle.SetOptStat(0);

	c_2d.cd()
	hbnnvsmela_zz.SetTitle("BNN/MELA");
	hbnnvsmela_zz.GetXaxis().SetTitle("BNN(x)");
	hbnnvsmela_zz.GetYaxis().SetTitle("MELA(x)");

	legend = TLegend(0.22, 0.75, 0.78, 0.88)

	hbnnvsmela_zx.SetMarkerColor(3)
#	hbnnvsmela_zx.Draw()
	hbnnvsmela_zz.SetMarkerColor(2)
	hbnnvsmela_sig.SetMarkerColor(4)
	hbnnvsmela_zz.SetFillColor(2)
	hbnnvsmela_sig.SetFillColor(4)
	hbnnvsmela_zx.SetFillColor(3)
	
	hbnnvsmela_zz.SetLineWidth(2)
	hbnnvsmela_zx.SetLineWidth(2)
	hbnnvsmela_sig.SetLineWidth(2)
	hbnnvsmela_zz.SetMarkerSize(2.);
	hbnnvsmela_sig.SetMarkerSize(2.);
	
	hbnnvsmela_zz.Draw("box")
#	hbnnvsmela_sig.SetMarkerColor(1)
	hbnnvsmela_sig.Draw("same box")
	legend.AddEntry(hbnnvsmela_zz, "ZZ: corr. " + str(round(hbnnvsmela_zz.GetCorrelationFactor(), 3)), "p");
#	legend.AddEntry(hbnnvsmela_zx, "Z+X", "p");
	legend.AddEntry(hbnnvsmela_sig, "sig.: corr. " + str(round(hbnnvsmela_sig.GetCorrelationFactor(), 3)), "p");
	legend.Draw("same")

	c_2d.Update()
	c_2d.SaveAs("bnn_vs_mela_" + evtset + ".gif")


	line = dat.readline()
	data = line.rsplit()
	columns = []
	for varname in data:
		print varname
		varname = varname.lstrip(':b')
		columns.append(varname); 

	print "btotal_: " + str(btotal_)
	print "Integral hb: " + str(hb.Integral())


	bkg.close()
	

	#------- BNN eff ---------
	
	geff = TGraph()
	geff.SetName("bnneff")
	gmelaeff = TGraph()
	gmelaeff.SetName("melaeff")
	
	stotal=hs.Integral();
	btotal=hb.Integral();
	cs = cdf(hs);
	cb = cdf(hb);
	
	maxsig = 0;
	optmvax = 0;
	optmvay = 0;
	
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

	
	#------- MELA eff ---------

	smelatotal=hmelas.Integral()
	bmelatotal=hmelab.Integral()
	
	cmelas = cdf(hmelas);
	cmelab = cdf(hmelab);
	
	maxsig = 0;
	optmvax = 0;
	optmvay = 0;
	
	for i in range(hmelas.GetNbinsX()):
		cs1 = smelatotal - (cmelas[i-1])
		cb1 = bmelatotal - (cmelab[i-1])
		es = cs1 / smelatotal
		eb = cb1 / bmelatotal
		sign=0
		if(cb1 > 0): 
			sign = cs1/sqrt(cs1 + cb1)
		if(sign > maxsig):
			maxsig = sign;
			optcut = i;	 
			optmvax = eb;
			optmvay = es;
		heffmela.Fill(eb, es);
		if(i < 2):
			continue
		else:
			gmelaeff.SetPoint(gmelaeff.GetN(), eb, es);

#		print str(i) + " eb: " + str(eb) + " es: " + str(es)
	print "\n Maximum Significance " + str(maxsig) + ", is bin: " + str(optcut) + ", which is MELA Cut: " + str(optcut/bins)
	print "es: " + str(optmvay) + "\teb: " + str(optmvax)


	print "stotal: " + str(stotal_)
	print "btotal: " + str(btotal_)
	print "stotal_sel: " + str(stotal_sel)
	print "btotal_sel: " + str(btotal_sel)


	ceffmela = TCanvas( 'ceffmela', 'effmela', 200, 10, 700, 500 )
	ceffmela.cd();
	gStyle.SetOptStat(0);
	heffmela.SetMarkerSize(2);

	ceffmela.Update();
	ceffmela.SaveAs("heffmela_" + evtset + ".gif")
	ceffmela.Update();
	
	ceff = TCanvas( 'ceff', 'eff', 200, 10, 700, 500 )
	ceff.cd();
	legend_x = TLegend(0.6, 0.25, 0.85, 0.45)
	gStyle.SetOptStat(0);
	heff.SetLineColor(2);
	geff.SetLineColor(2);
	geff.SetLineWidth(2);

	heff.SetTitle("ROC for BNN & MELA (bkg = ZZ & Z+X)");
	heff.GetXaxis().SetTitle("#epsilon_{bkg.}");
	heff.GetXaxis().SetTitleSize(0.06);
	heff.GetYaxis().SetTitleSize(0.06);
	heff.GetXaxis().SetTitleOffset(0.65);
	heff.GetYaxis().SetTitleOffset(0.65);
	heff.GetYaxis().SetTitle("#epsilon_{sig.}");
	heff.Draw();
	geff.Draw("same");
	heffmela.SetLineColor(3);
	gmelaeff.SetLineColor(3);
	gmelaeff.SetLineWidth(2);

	heffmela.Draw("same");
	gmelaeff.Draw("same")
	legend_x.AddEntry(heffmela, "MELA", "l");
	legend_x.AddEntry(heff, "BNN", "l");
	legend_x.Draw("same")
	ceff.Update();
	ceff.SaveAs("heff_" + evtset + ".gif")
	ceff.Update();


	#------- BNN Dist ---------

#	c1.SetLogy();
	hb.Scale(1/hb.Integral())
	hmelab.Scale(1/hmelab.Integral())

	c1.cd()
	gStyle.SetOptStat(0);
	hs.Rebin(2);
	hb.Rebin(2);
	hs.SetTitle("BNN Discriminant: sig. & bkg (ZZ & Z+X)");
	hs.GetXaxis().SetTitle("D(x)");
	hs.SetLineWidth(2)
	hb.SetLineWidth(2)
	hs.GetXaxis().SetTitleSize(0.06);
	hs.GetYaxis().SetTitleSize(0.06);
	hs.GetXaxis().SetTitleOffset(0.65);
	hs.GetYaxis().SetTitleOffset(0.65);
	hs.Draw()
	hb.Draw("same")
	c1.Update()
	c1.SaveAs("hdist_bnn_" + evtset + ".gif")


	#------- MELA Dist ---------

	c2 = TCanvas( 'c2', 'mva output', 200, 10, 700, 500 )
	c2.cd()
	gStyle.SetOptStat(0);
	c2.SetGridx()
	c2.SetGridy()
	c2.Draw()
#	c2.SetLogy();
#	hmelas.Scale(1/hmelas.Integral())
#	hmelab.Scale(1/hmelab.Integral())
	gStyle.SetOptStat(0);
	hmelas.Rebin(2);
	hmelab.Rebin(2);
	hmelab.SetTitle("MELA Discriminant: sig. & bkg (ZZ & Z+X)");
	hmelab.GetXaxis().SetTitle("D(x)");
	hmelab.SetLineWidth(2)
	hmelas.SetLineWidth(2)
	hmelab.GetXaxis().SetTitleSize(0.06);
	hmelab.GetYaxis().SetTitleSize(0.06);
	hmelab.GetXaxis().SetTitleOffset(0.65);
	hmelab.GetYaxis().SetTitleOffset(0.65);
	hmelab.Draw()
	hmelas.Draw("same")

	c2.Update()
	c2.SaveAs("hdist_mela_" + evtset + ".gif")


	
sig = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_mela.dat')
bkg = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_mela.dat')

eventloop(sig, bkg, "2e2mu", "zx_cuts")