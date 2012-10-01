from string import *
from ROOT import *
from time import sleep

import os, sys, re, bnn_2e2mu_loosesel3_ext, bnn_4e_loosesel3_ext, bnn_4mu_loosesel3_ext, applycuts_ext
from applycuts_ext import applycuts
from ROOT import TMVA
import array



gROOT.Reset()
#gROOT.ProcessLine('.L ~/work/MELA/scripts/MELA.C')

gROOT.LoadMacro(".L TMVAClassification_HZZ_2e2mu_MLP.class.C")




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

reader3 = TMVA.Reader()
var3_ = []
for i, var3 in enumerate(mvavars1):
	var3_.append(array.array('f',[0]))
	reader3.AddVariable(var3,var3_[i])
reader3.BookMVA("MLP","weights/TMVAClassification_HZZ_2e2mu_MLP.weights.xml")


#weight	weight	ele1_pt	ele1_eta	ele1_phi	ele1_charge	ele1_trackIso	ele1_EcalIso	ele1_HcalIso	ele1_X	ele1_SIP	ele2_pt	ele2_eta	ele2_phi	ele2_charge	ele2_trackIso	ele2_EcalIso	ele2_HcalIso	ele2_X	ele2_SIP	ele3_pt	ele3_eta	ele3_phi	ele3_charge	ele3_trackIso	ele3_EcalIso	ele3_HcalIso	ele3_X	ele3_SIP	ele4_pt	ele4_eta	ele4_phi	ele4_charge	ele4_trackIso	ele4_EcalIso	ele4_HcalIso	ele4_X	ele4_SIP	

#worst_iso_X	second_worst_iso_X	worst_vertex	second_worst_vertex	mZ	mZstar	mbestH	index	channel	sample	end

def loop(hweights, sig, hs, c1):


	#------- sig events ---------
	stotal_ = 0
	stotal_sel = 0
	sigweight = 0
	i = 0 
	line = sig.readline()
	data = line.rsplit()
	columns = []
	for varname in data:
		print varname
		varname = varname.lstrip(':b')
		columns.append(varname); 
	
	#do signal
	for line in sig:
		data = line.rsplit()
		data_ = map(float, line.rsplit())
		d = dict(zip(columns, data))
		inputs = []

		if(float(d["mass4l"]) < 110):
			continue
		if(float(d["mass4l"]) > 150):
			continue
	
		
		for i, var3 in enumerate(mvavars1):
			var3_[i][0] = float(d[var3])
#			print var3 + ": " + str(var3_[i][0])


		channel_ = int(d["channel"])
		weight = float(d["weight"])
		if(int(d["channel"]) == 1):
			for var in mvavars1:
				inputs.append(float(d[var]))
			y = bnn_4mu_loosesel3_ext.bnn_4mu_loosesel3(inputs, 100, 200);
		if(int(d["channel"]) == 2):
			for var in mvavars2:
				inputs.append(float(d[var]))
			y = bnn_4e_loosesel3_ext.bnn_4e_loosesel3(inputs, 100, 200);
		if(int(d["channel"]) == 3):
			for var in mvavars3:
				inputs.append(float(d[var]))
			y = bnn_2e2mu_loosesel3_ext.bnn_2e2mu_loosesel3(inputs, 100, 200);

#		y = reader3.EvaluateMVA("MLP")

		weight = float(d["weight"])
		hweights.Fill(float(d["Z2mass"]))
		stotal_ += float(d["weight"])
		if(applycuts(columns, data_, channel_)):
			stotal_sel += float(d["weight"])
#			continue
		else:
			hs.Fill(y, weight)
			if(int(d["index"]) == 6):
				sigweight+=weight

		i+=1
		if( not (i%1000)):
			print i
			print y	
			hs.Draw()
			c1.Update()
	sig.close()
	
	print "Higgs 125 weight: " + str(sigweight)
	print "Integral hs: " + str(hs.Integral())
	
def calceff(geff, heff, hs, hb):
	
	stotal=hs.Integral();
	btotal=hb.Integral();
	cs = cdf(hs);
	cb = cdf(hb);
	
	maxsig = 0;
	optmvax = 0;
	optmvay = 0;
	optcut = 0;
	
	bins = hs.GetNbinsX()
	
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
	
	

def eventloop(sig_test,bkg_test, sig_train, bkg_train,channel,evtset):

	bins = 50
	
	hs_test = TH1F('hs_test', 'hs_test', bins, 0, 1)
	hs_train = TH1F('hs_train', 'hs_train', bins, 0, 1)
	hb_test = TH1F('hb_test', 'hb_test', bins, 0, 1)
	hb_train = TH1F('hb_test', 'hb_test', bins, 0, 1)

	hweights1 = TH1F('weights1', 'weights1', 50, 0, 100)
	hweights2 = TH1F('weights2', 'weights2', 50, 0, 100)
	hweights3 = TH1F('weights2', 'weights2', 50, 0, 100)
	hweights4 = TH1F('weights2', 'weights2', 50, 0, 100)

	hs_test.SetLineColor(2);
	hb_test.SetLineColor(1);

	hs_train.SetLineColor(2);
	hb_train.SetLineColor(1);
	hs_train.Draw("lep same")
	hb_train.Draw("lep same")		
	
	c1 = TCanvas( 'c1', 'mva output', 200, 10, 700, 500 )
	c1.cd()
	c1.SetGridx()
	c1.SetGridy()
#	c1.SetLogy()
	c1.Draw()
	gStyle.SetOptStat(0);


	loop(hweights1, sig_test, hs_test, c1)
	loop(hweights2, sig_train, hs_train, c1)
	loop(hweights3, bkg_test, hb_test, c1)
	loop(hweights4, bkg_train, hb_train, c1)
	
 	print "Integral hb_test: " + str(hb_test.Integral())
 	print "Integral hb_train: " + str(hb_train.Integral())
 	print "Integral hs_test: " + str(hb_test.Integral())
 	print "Integral hs_train: " + str(hb_train.Integral())
	
	
	hs_test.Draw()
	hb_test.Draw("same")		
	hs_train.Draw("lep same")
	hb_train.Draw("lep same")
	c1.Update()
	c1.SaveAs("hdist_bnn_train_vs_test_" + channel + "_" + evtset + ".gif")


	heff_test = TH2F('heff_test', 'heff_test', bins, 0, 1.0,bins, 0.0, 1)
	heff_train = TH2F('heff_train', 'heff_train', bins, 0, 1.0,bins, 0.0, 1)

	geff_test = TGraph()
	geff_test.SetName("bnneff_test")

	geff_train = TGraph()
	geff_train.SetName("bnneff_train")

	calceff(geff_test, heff_test, hs_test, hb_test)
	calceff(geff_train, heff_train, hs_train, hb_train)


	ceff = TCanvas( 'ceff', 'eff', 200, 10, 700, 500 )
	ceff.cd();
	legend = TLegend(0.6, 0.25, 0.85, 0.45)
	gStyle.SetOptStat(0);
	heff_test.SetLineColor(2);
	geff_test.SetLineColor(2);
	geff_test.SetLineWidth(2);

	heff_test.SetTitle("ROC for BNN (overt");
	heff_test.GetXaxis().SetTitle("#epsilon_{bkg.}");
	heff_test.GetXaxis().SetTitleSize(0.06);
	heff_test.GetYaxis().SetTitleSize(0.06);
	heff_test.GetXaxis().SetTitleOffset(0.65);
	heff_test.GetYaxis().SetTitleOffset(0.65);
	heff_test.GetYaxis().SetTitle("#epsilon_{sig.}");
	heff_test.Draw();
	geff_test.Draw("same");
	heff_train.SetLineColor(3);
	geff_train.SetLineColor(3);
	geff_train.SetLineWidth(2);

	heff_train.Draw("same");
	geff_train.Draw("same")
	legend.AddEntry(heff_train, "training events", "l");
	legend.AddEntry(heff_test, "testing events", "l");
	legend.Draw("same")
	ceff.Update();
	ceff.SaveAs("heff_" + channel + "_" + evtset + ".gif")
	ceff.Update();




	sleep(10)	


	
sig_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_4e_test_uw.dat')
bkg_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_4e_test_uw.dat')
sig_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_4e_train_uw.dat')
bkg_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_4e_train_uw.dat')

eventloop(sig_test, bkg_test, sig_train, bkg_train, "4e", "3")


sig_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_4mu_test_uw.dat')
bkg_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_4mu_test_uw.dat')
sig_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_4mu_train_uw.dat')
bkg_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_4mu_train_uw.dat')

eventloop(sig_test, bkg_test, sig_train, bkg_train, "4mu", "3")

sig_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_2e2mu_test_uw.dat')
bkg_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_2e2mu_test_uw.dat')
sig_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/sig_bnn_2e2mu_train_uw.dat')
bkg_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_loosesel/bkg_bnn_2e2mu_train_uw.dat')

eventloop(sig_test, bkg_test, sig_train, bkg_train, "2e2mu", "3")

sig_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_s_test.dat')
bkg_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/new/sig_ps_test.dat')
sig_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/sig_bnn_train_s_uw.dat')
bkg_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_mela/sig_bnn_train_ps_uw.dat')

eventloop(sig_test, bkg_test, sig_train, bkg_train, "p_ps", "yes")
