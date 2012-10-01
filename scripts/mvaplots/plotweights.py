from string import *
from ROOT import *
from time import sleep
import os, sys, re, bnn_2e2mu_fullsel_ext, bnn_4e_fullsel_ext, bnn_4mu_fullsel_ext, applycuts_ext
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
#		data_ = [float(x) for x in line.rsplit()]
		data_ = map(float, line.rsplit())
		d = dict(zip(columns, data))
		inputs = []

		if(float(d["mass4l"]) < 110):
			continue
		if(float(d["mass4l"]) > 150):
			continue

		weight = float(d["weight"])
		hweights.Fill(float(d["Z2mass"]))
		stotal_ += float(d["weight"])
		hs.Fill(weight)

		i+=1
		if( not (i%1000)):
			print i
			hs.Draw()
			c1.Update()
	sig.close()
	
	print "Higgs 125 weight: " + str(sigweight)
	print "Integral hs: " + str(hs.Integral())
	


def eventloop(sig_test,bkg_test, sig_train, bkg_train,channel,evtset):

	bins = 50
	upper = 0.0015
	lower = 0.00001
	
	hs_test = TH1F('hs_test', 'hs_test', bins, lower, upper)
	hs_train = TH1F('hs_train', 'hs_train', bins, lower, upper)
	hb_test = TH1F('hb_test', 'hb_test', bins, lower, upper)
	hb_train = TH1F('hb_test', 'hb_test', bins, lower, upper)

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
#	c1.SetGridy()
#	c1.SetLogy()
#	c1.SetLogx()
	c1.Draw()
	gStyle.SetOptStat(0);


	loop(hweights1, sig_test, hs_test, c1)
	loop(hweights2, sig_train, hs_train, c1)
	loop(hweights3, bkg_test, hb_test, c1)
	loop(hweights4, bkg_train, hb_train, c1)
	
 	print "Integral hb_test: " + str(hb_test.Integral())
 	print "Integral hb_train: " + str(hb_train.Integral())
 	print "Integral hs_test: " + str(hs_test.Integral())
 	print "Integral hs_train: " + str(hs_train.Integral())
	
	
	

	
#	hs_test.Draw()
	hb_test.Draw("")		
#	hs_train.Draw("lep same")
	hb_train.Draw("lep same")
	c1.Update()
	c1.SaveAs("hweights_bkg.gif")


	c1 = TCanvas( 'c1', 'mva output', 200, 10, 700, 500 )
	c1.cd()
	hs_test.Draw()
	hs_train.Draw("lep same")
	c1.SetGridx()
	c1.Draw()
	c1.SaveAs("hweights_sig.gif")



	sleep(10)	


	
sig_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_fullsel/sig_bnn_test2.dat')
bkg_test = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_fullsel/bkg_bnn_test2.dat')
sig_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_fullsel/sig_bnn_train2.dat')
bkg_train = open('/home/jbochenek/data/HZZ4l_2012_data/hzz_fullsel/bkg_bnn_train2.dat')

eventloop(sig_test, bkg_test, sig_train, bkg_train, "2e2mu", "test")