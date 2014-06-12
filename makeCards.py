#!/usr/bin/env python

import optparse
parser = optparse.OptionParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_option("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_option("--signalInjectMass", help="Mass For Injected Signal",type=float,default=125.0)
parser.add_option("--toyData", help="Make Toy Data from PDFs for data_obs",action="store_true",default=False)
parser.add_option("--cutOpt", help="Creates Cards with Different Cut Values",action="store_true",default=False)
parser.add_option("--gaussian", help="Use A Gaussian Signal Template with floating width",type=float,default=-1.0)
parser.add_option("-m","--higgsMass", help="Use This Higgs Mass",type=float,default=-1.0)
parser.add_option("--combinationsOnly", help="Run Only Combinations of Channels",action="store_true",default=False)
parser.add_option("--forCombGroup",help="Removes Debug text to send to combination group",action="store_true",default=False)
args, fakeargs = parser.parse_args()

import math
import ROOT as root
from helpers import *
import datetime
import sys
import os
import os.path
import stat
import copy
import random
import shutil
import multiprocessing
import time
import string
myThread = multiprocessing.Process

from ROOT import gSystem
gSystem.Load('libRooFit')

from effReaderFromFile import *
from signalPars import *
#from signalParsCB import *

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

#from signalfits import getRooFitSignalPars as sigFits
#from signalfitsNoMuScle import getRooFitSignalPars as sigFits
#effReader = EfficiencyReader()

NPROCS = 4
RUNSIMPLELIMITS = False
USEGPANNA = False

SIGNALFIT = [110.,140.]

USETREES=True
HISTNAME="mDiMu"

BAKPARAMUNC=False

if args.cutOpt:
  USEGPANNA = False

from xsec import *

def vetoOutOfBoundsEvents(hist,boundaries=[]):
  xbinLow = None
  xbinHigh = None
  if len(boundaries)==2:
    xbinLow, xbinHigh = getXbinsHighLow(hist,boundaries[0],boundaries[1])
  else:
    print("Error: vetoOutOfBoundsEvents: boundaries must be length 2, exiting.")
    sys.exit(1)
  for i in range(0,xbinLow):
    hist.SetBinContent(i,0.0)
    hist.SetBinError(i,0.0)
  for i in range(xbinHigh+1,hist.GetNbinsX()+2):
    hist.SetBinContent(i,0.0)
    hist.SetBinError(i,0.0)

def convertSigName(name,energyStr):
  if "ggH" in name:
    name="ggH"
  elif "vbfH" in name:
    name="qqH"
  elif "whH" in name:
    name="WH"
  elif "zhH" in name:
    name="ZH"
  name += "_hmm"+energyStr
  return name

def convertUncName(name):
  if "lumi" in name:
    return name
  if "QCDscale" in name:
    return name
  if "CMS_scale" in name:
    return name
  if "CMS_res" in name:
    return name
  if "CMS_eff" in name:
    return name
  if "pdf_gg" == name:
    return name
  if "pdf_qqbar" == name:
    return name
  if "UEPS" == name:
    return name
  name = "CMS_hmm_"+name
  return name

class Param:
  def __init__(self,name,nominal,lowErr,highErr):
    self.name = name
    self.nominal = nominal
    self.lowErr = lowErr
    self.highErr = highErr
  def __str__(self):
    return "{0}: {1:.3g} +{2:.3g} -{3:.3g}".format(self.name,self.nominal,self.lowErr,self.highErr)
  def __repr__(self):
    return str(self)
  def getErrString(self):
    if self.lowErr == self.highErr:
        return "{0:.5g}".format(self.lowErr)
    else:
        return "-{0:.5g}/+{1:.5g}".format(self.lowErr,self.highErr)

###################################################################################

def makePDFBakBernsteinProd(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakBernsteinProd: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    if order == None:
      if "Jets01PassPtG10BB" in name:
        order = 5
      elif "Jets01PassPtG10BO" in name:
        order = 6
      elif "Jets01PassPtG10BE" in name:
        order = 6
      elif "Jets01PassPtG10OO" in name:
        order = 5
      elif "Jets01PassPtG10OE" in name:
        order = 5
      elif "Jets01PassPtG10EE" in name:
        order = 5
      elif "Jets01FailPtG10BB" in name:
        order = 5
      elif "Jets01FailPtG10BO" in name:
        order = 5
      elif "Jets01FailPtG10BE" in name:
        order = 5
      elif "Jets01FailPtG10OO" in name:
        order = 6
      elif "Jets01FailPtG10OE" in name:
        order = 5
      elif "Jets01FailPtG10EE" in name:
        order = 6
      elif "Jet2CutsVBFPass" in name:
        order = 2
      elif "Jet2CutsGFPass" in name:
        order = 5
      elif "Jet2CutsFailVBFGF" in name:
        order = 5
      else:
        print "Error: makePDFBakBernsteinProd: category '"+str(name)+"' not recognized, exiting."
	sys.exit(1)
	

    rooParamList = []
    rooArgList = root.RooArgList()
    for i in range(order+1):
      tmpArg = root.RooRealVar(channelName+"_B"+str(i),"Bernstein Coefficient "+str(i), 0.0, 0., 1.)
      rooArgList.add(tmpArg)
      rooParamList.append(tmpArg)

    debug += "#    Bernstein Order: "+str(order)+"\n"
    debug += "#    pdfArgs: "+dimuonMass.GetName()+" "
    for i in rooParamList:
        debug += i.GetName()+" "
    debug += "\n"

    #print
    #print "Bernstein Order: ",order
    #for i in rooParamList:
    #    i.Print()
    #print
  
    pdfMmumu = root.RooBernstein("bak","Bernstein Order: "+str(order),dimuonMass,rooArgList)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    paramList = [i.GetName() for i in rooParamList]
    for param in paramList:
        param.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    #for i in rooParamList:
    #  debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    #debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, order

def makePDFBakExpMOverSq(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakExpMOverSq: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    #InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    #ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -1., 1.)

    # new values from Anna 13 Jun 2013
    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91.187, 30., 105.)
    ExpMass = root.RooRealVar(channelName+"_ExpMass","ExpMass", 0.0, -2., 2.)
  
    if ('Jet2CutsVBFPass' in name ):
      debug += "###  fixing InvPolMass to Z pdg value\n"
      InvPolMass.setConstant(True)

    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(@0*@2)/(@0-@1)/(@0-@1)",root.RooArgList(dimuonMass,InvPolMass,ExpMass))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [InvPolMass,ExpMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(True)
        if ('Jet2CutsVBFPass' in name ):
          InvPolMass.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
   #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

def makePDFBakExpMOverSqP0(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakExpMOverSqP0: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name


    # new values from Anna 13 Jun 2013
    EXPPOLP0_1 = root.RooRealVar(channelName+"_EXPPOLP0_1","EXPPOLP0_1", 92, 30., 105.)
    EXPPOLP0_2 = root.RooRealVar(channelName+"_EXPPOLP0_2","EXPPOLP0_2", -0.5, -2.5, 2.5)
    EXPPOLP0_3 = root.RooRealVar(channelName+"_EXPPOLP0_3","EXPPOLP0_3", 1., -5., 5.)
  
    #if ('Jet2CutsVBFPass' in name ):
    #  debug += "###  fixing InvPolMass to Z pdg value\n"
    #  InvPolMass.setConstant(True)

    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(-@0*@2*@2)/(@0-@1)*(1./(@0-@1)+@3*@3*@0)",root.RooArgList(dimuonMass,EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(True)
        if ('Jet2CutsVBFPass' in name ):
          InvPolMass.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
   #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

def makePDFBakExpMOverSqP0New(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakExpMOverSqP0New: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name


    # new values from Anna 13 Jun 2013
    EXPPOLP0_1 = root.RooRealVar(channelName+"_EXPPOLP0New_1","EXPPOLP0New_1", 91, 30., 109.)
    EXPPOLP0_2 = root.RooRealVar(channelName+"_EXPPOLP0New_2","EXPPOLP0New_2", 0.2, -2.5, 2.5)
    EXPPOLP0_3 = root.RooRealVar(channelName+"_EXPPOLP0New_3","EXPPOLP0New_3", 1., -20., 20.)
    #EXPPOLP0_4 = root.RooRealVar(channelName+"_EXPPOLP0_4","EXPPOLP0_4", 0.15, -2.5, 2.5)
  

    #v0, p3 is very undefined 
    #pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(-@0*@2*@2)/(@0-@1)*(1./(@0-@1)+@3*@3*@0)",root.RooArgList(dimuonMass,EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3))

    #v1 interesting, may be the best one, good to fix p1 only for VBFtight, for other bad 
    pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(-@0*@2*@2)/(@0-@1)/(@0-@1)+@3*@3*TMath::Exp(-@0*@2*@2)",root.RooArgList(dimuonMass,EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3))
    if ('Jet2CutsVBFPass' in name ):
      EXPPOLP0_1.setConstant(True)

    #v2 <-bad results
    #pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(-@0*@2*@2)/(@0-@1)*(1./(@0-@1)+@3*@3)",root.RooArgList(dimuonMass,EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3))

    #v3 interesting, may be the best one, not good to fix p1 for BO at least
    #pdfMmumu = root.RooGenericPdf("bak","TMath::Exp(-@0*@2*@2)/(@0-@1)/(@0-@1)+@3*@3*TMath::Exp(-@0*@4*@4)",root.RooArgList(dimuonMass,EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3,EXPPOLP0_4))
    #EXPPOLP0_1.setConstant(True)

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3]
    #rooParamList = [EXPPOLP0_1,EXPPOLP0_2,EXPPOLP0_3,EXPPOLP0_4]

    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(True)
        if ('Jet2CutsVBFPass' in name ):
          InvPolMass.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
   #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None


def makePDFBakMOverSq(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakMOverSq: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    InvPolMass = root.RooRealVar(channelName+"_InvPolMass","InvPolMass", 91., 70., 95.)
    pdfMmumu = root.RooGenericPdf("bak","@0/(@0-@1)/(@0-@1)",root.RooArgList(dimuonMass,InvPolMass))

    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [InvPolMass]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
        param.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

def makePDFBakMSSM(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakMSSM: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    bwWidth = root.RooRealVar(channelName+"_bwWidth","bwWidth",5,0,30)
    bwmZ = root.RooRealVar(channelName+"_bwmZ","bwmZ",91,85,95)
    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1e-03,-1e-01,1e-01)
    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)

    phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(dimuonMass,expParam))
    bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(dimuonMass,bwmZ,bwWidth,expParam))
    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(bwExpMmumu,phoExpMmumu),root.RooArgList(mixParam))
    

    # Just For Z-Peak Part
    assert(dimuonMassZ != None)
    assert(rooDatasetZ != None)

    bwMmumuZ  = root.RooGenericPdf("bwMmumuZ","(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(dimuonMassZ,bwmZ,bwWidth))
    bwMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
    bwmZ.setConstant(True)
    bwWidth.setConstant(True)

    ### Debug Time
    #frameZ = dimuonMassZ.frame()
    #frameZ.SetName("bak_PlotZ")
    #rooDatasetZ.plotOn(frameZ)
    #bwMmumuZ.plotOn(frameZ)
    #canvas = root.TCanvas()
    #frameZ.Draw()
    #canvas.SaveAs("debug_"+name+channelName+"_Z.png")

    # Back to everywhere else
   
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)
    fr.Print()

    rooParamList = [bwmZ,bwWidth,expParam,mixParam]
    paramList = [expParam.GetName(),mixParam.GetName()]

    for param in rooParamList:
      param.setConstant(True)

    ## Let them float!
    expParam.setConstant(False)
    mixParam.setConstant(False)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    ### Debug Time
    #frame = dimuonMass.frame()
    #frame.SetName("bak_Plot")
    #rooDataset.plotOn(frame)
    ##pdfMmumu.plotOn(frame,root.RooFit.Range(110,160))
    #pdfMmumu.plotOn(frame,root.RooFit.Range(minMass,maxMass))
    #canvas = root.TCanvas()
    #frame.Draw()
    #canvas.SaveAs("debug_"+name+channelName+".png")

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

def makePDFBakMSSMPrime(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakMSSMPrime: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    bwWidth = root.RooRealVar(channelName+"_bwWidth","bwWidth",5,0,30)
    bwmZ = root.RooRealVar(channelName+"_bwmZ","bwmZ",91,85,95)
    expParam = root.RooRealVar(channelName+"_expParam","expParam",-1e-03,-1e-01,1e-01)
    phoCoef = root.RooRealVar(channelName+"_phoCoef","Photon Term Coef",0.5,0,1)
    bwCoef = root.RooRealVar(channelName+"_bwCoef","Breit-Wigner Term Coef",0.5,0,1)

    phoExpMmumu = root.RooGenericPdf("phoExpMmumu","exp(@0*@1)*pow(@0,-2)",root.RooArgList(dimuonMass,expParam))
    bwExpMmumu  = root.RooGenericPdf("bwExpMmumu","exp(@0*@3)*(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(dimuonMass,bwmZ,bwWidth,expParam))
    phoExpMmumuE = root.RooExtendPdf("phoExpMmumuE",phoExpMmumu,phoCoef)
    bwExpMmumuE = root.RooExtendPdf("bwExpMmumuE",bwExpMmumu,bwCoef)
    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(bwExpMmumuE,phoExpMmumuE))

    # Just For Z-Peak Part
    assert(dimuonMassZ != None)
    assert(rooDatasetZ != None)

    bwMmumuZ  = root.RooGenericPdf("bwMmumuZ","(@2)/(pow(@0-@1,2)+0.25*pow(@2,2))",root.RooArgList(dimuonMassZ,bwmZ,bwWidth))
    bwMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
    bwmZ.setConstant(True)
    bwWidth.setConstant(True)

    ### Debug Time
    #frameZ = dimuonMassZ.frame()
    #frameZ.SetName("bak_PlotZ")
    #rooDatasetZ.plotOn(frameZ)
    #bwMmumuZ.plotOn(frameZ)
    #canvas = root.TCanvas()
    #frameZ.Draw()
    #canvas.SaveAs("debug_"+name+channelName+"_Z.png")

    # Back to everywhere else
   
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)
    fr.Print()

    rooParamList = [bwmZ,bwWidth,expParam,bwCoef,phoCoef]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
      param.setConstant(True)

    bwWidth.setConstant(True)
    bwmZ.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    ### Debug Time
    #frame = dimuonMass.frame()
    #frame.SetName("bak_Plot")
    #rooDataset.plotOn(frame)
    ##pdfMmumu.plotOn(frame,root.RooFit.Range(110,160))
    #pdfMmumu.plotOn(frame,root.RooFit.Range(minMass,maxMass))
    #canvas = root.TCanvas()
    #frame.Draw()
    #canvas.SaveAs("debug_"+name+channelName+".png")

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None


def makePDFBakVoigtPMm2(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakVoigtPMm2: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitWidth.setConstant(True)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",91,85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",1.5,0.0,30.0)
    voitMmumu = root.RooVoigtian(channelName+"bak_voitMmumu","voitMmumu",dimuonMass,voitmZ,voitWidth,voitSig)
    mm2Mmumu = root.RooGenericPdf(channelName+"bak_Mm2","pow(@0,-2)",root.RooArgList(dimuonMass))
    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)
    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,mm2Mmumu),root.RooArgList(mixParam))

    # Just For Z-Peak Part
    assert(dimuonMassZ != None)
    assert(rooDatasetZ != None)
    voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",dimuonMassZ,voitmZ,voitWidth,voitSig)

    voitMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    ### Debug Time
    #frameZ = dimuonMassZ.frame()
    #frameZ.SetName("bak_PlotZ")
    #rooDatasetZ.plotOn(frameZ)
    #voitMmumuZ.plotOn(frameZ)
    #canvas = root.TCanvas()
    #frameZ.Draw()
    #canvas.SaveAs("debug_"+name+channelName+"_Z.png")

    # Back to everywhere else
   
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)
    fr.Print()

    rooParamList = [voitmZ,voitWidth,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
      param.setConstant(True)

    voitWidth.setConstant(True)
    voitmZ.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    ### Debug Time
    #frame = dimuonMass.frame()
    #frame.SetName("bak_Plot")
    #rooDataset.plotOn(frame)
    ##pdfMmumu.plotOn(frame,root.RooFit.Range(110,160))
    #pdfMmumu.plotOn(frame,root.RooFit.Range(minMass,maxMass))
    #canvas = root.TCanvas()
    #frame.Draw()
    #canvas.SaveAs("debug_VoigtPMm2_"+name+channelName+".png")

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None

def makePDFBakVoigtPExpMm2(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    debug = ""
    debug += "### makePDFBakVoigtPExpMm2: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitWidth.setConstant(True)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",91,85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",1.5,0.0,30.0)
    voitMmumu = root.RooVoigtian(channelName+"bak_voitMmumu","voitMmumu",dimuonMass,voitmZ,voitWidth,voitSig)
    expParam = root.RooRealVar(channelName+"_expParam","expParam",0.01,-1,1)
    expMm2Mmumu = root.RooGenericPdf(channelName+"bak_expMm2","exp(-1*@0*@1*@1)*pow(@0,-2)",root.RooArgList(dimuonMass,expParam))
    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)
    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMm2Mmumu),root.RooArgList(mixParam))

    # Just For Z-Peak Part
    assert(dimuonMassZ != None)
    assert(rooDatasetZ != None)
    voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",dimuonMassZ,voitmZ,voitWidth,voitSig)

    voitMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    ### Debug Time
    #frameZ = dimuonMassZ.frame()
    #frameZ.SetName("bak_PlotZ")
    #rooDatasetZ.plotOn(frameZ)
    #voitMmumuZ.plotOn(frameZ)
    #canvas = root.TCanvas()
    #frameZ.Draw()
    #canvas.SaveAs("debug_"+name+channelName+"_Z.png")

    # Back to everywhere else
   
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)
    fr.Print()

    rooParamList = [voitmZ,voitWidth,expParam,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
      param.setConstant(True)

    voitWidth.setConstant(True)
    voitmZ.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    ### Debug Time
    #frame = dimuonMass.frame()
    #frame.SetName("bak_Plot")
    #rooDataset.plotOn(frame)
    ##pdfMmumu.plotOn(frame,root.RooFit.Range(110,160))
    #pdfMmumu.plotOn(frame,root.RooFit.Range(minMass,maxMass))
    #canvas = root.TCanvas()
    #frame.Draw()
    #canvas.SaveAs("debug_VoigtPExpMm2_"+name+channelName+".png")

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None


def makePDFBakOld(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,dimuonMassZ=None,rooDatasetZ=None,order=None):
    #print "GP's debug line *******************************************************\n"
    debug = ""
    debug += "### makePDFBakOld: "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    channelName = name

    voitWidth = root.RooRealVar(channelName+"_voitWidth","voitWidth",2.4952)
    voitWidth.setConstant(True)
    voitmZ = root.RooRealVar(channelName+"_voitmZ","voitmZ",91,85,95)
    voitSig = root.RooRealVar(channelName+"_voitSig","voitSig",1.5,0.0,30.0)
    voitMmumu = root.RooVoigtian("bak_voitMmumu","voitMmumu",dimuonMass,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar(channelName+"_expParam","expParam",0.1,-5,5)
    #expMmumu = root.RooExponential("bak_expMmumu","expMmumu",dimuonMass,expParam*expParam)
    # Generic PDF
    expMmumu = root.RooGenericPdf("bak_expMmumu","expMmumu","TMath::Exp(-1*@0*@1*@1)",root.RooArgList(dimuonMass,expParam))

    mixParam = root.RooRealVar(channelName+"_mixParam","mixParam",0.5,0,1)

    pdfMmumu = root.RooAddPdf("bak","bak",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))

    # Just For Z-Peak Part
    assert(dimuonMassZ != None)
    assert(rooDatasetZ != None)
    voitMmumuZ = root.RooVoigtian("bak_voitMmumuZ","voitMmumuZ",dimuonMassZ,voitmZ,voitWidth,voitSig)

    voitMmumuZ.fitTo(rooDatasetZ,root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

#    ## Debug Time
#    frameZ = dimuonMassZ.frame()
#    frameZ.SetName("bak_PlotZ")
#    dimuonMassZRooDataHist.plotOn(frameZ)
#    voitMmumuZ.plotOn(frameZ)
#    canvas = root.TCanvas()
#    frameZ.Draw()
#    saveAs(canvas,"debug_bakZ")

    # Back to everywhere else

    #expMmumu.fitTo(rooDataset,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expMmumu.fitTo(rooDataset,root.RooFit.Range("exprange"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    #expParam.setConstant(True)
   
    #fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.Range("whole"),root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True))
    fr.SetName("bak"+"_fitResult")
    #chi2 = pdfMmumu.createChi2(rooDataset)

    rooParamList = [voitmZ,voitSig,expParam,mixParam]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]

    for param in rooParamList:
      param.setConstant(True)

    voitWidth.setConstant(True)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      workspaceImportFn(fr)

    #Norm Time
    bakNormTup = None
    if True:
      wholeIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal,low,high"))
      signalIntegral = pdfMmumu.createIntegral(root.RooArgSet(dimuonMass),root.RooFit.Range("signal"))
      signalRangeList = getRooVarRange(dimuonMass,"signal")
      getSidebandString = "dimuonMass < {0} || dimuonMass > {1}".format(*signalRangeList)
      nSideband =  rooDataset.sumEntries(getSidebandString)
      nData =  rooDataset.sumEntries()
      bakNormTup = (nSideband,1.0/(1.0-signalIntegral.getVal()/wholeIntegral.getVal()))
      if nData > 0:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, predicted error: {1:.2%} true error: {2:.2%}".format(getSidebandString,1.0/sqrt(bakNormTup[0]),(bakNormTup[0]*bakNormTup[1] - nData)/nData))
      else:
        print("Gets Bak Norm Assuming Signal region is: {0} GeV, nData=0.0".format(getSidebandString))
    #print("nData: {0}, nPredict: {1}, nSideBand: {2}, alpha: {3}".format(
    #        nData, bakNormTup[0]*bakNormTup[1], bakNormTup[0], bakNormTup[1]))

    #rooDataset2 = rooDataset.reduce(root.RooFit.CutRange("low,signal,high"))
    #rooDataset2.SetName("bak_TemplateNoVeryLow")
    #if workspaceImportFn != None:
    #  workspaceImportFn(rooDataset2)

#    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName("bak_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    canvas.SaveAs("debug_"+name+channelName+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
    debug += "#    Bak Norm Tuple: {0:.2f} {1:.2f}\n".format(*bakNormTup)

    return paramList, bakNormTup, debug, None


def makePDFSigCBPlusGaus(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,channelName,energyStr,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigCBPlusGaus: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    mean = root.RooRealVar(channelName+"_"+name+"_Mean",channelName+"_"+name+"_Mean",125.,100.,150.)
    width = root.RooRealVar(channelName+"_"+name+"_Width",channelName+"_"+name+"_Width",5.0,0.5,20.0)
    width2 = root.RooRealVar(channelName+"_"+name+"_Width2",channelName+"_"+name+"_Width2",5.0,0.1,20.0)
    alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",channelName+"_"+name+"_Alpha",1.0,0.1,10.0)
    n = root.RooRealVar(channelName+"_"+name+"_n",channelName+"_"+name+"_n",1.0,0.1,10.0)
    mix = root.RooRealVar(channelName+"_"+name+"_mix",channelName+"_"+name+"_mix",0.5,0.0,1.0)
    cb = root.RooCBShape(name+"_CB",name+"_CB",dimuonMass,mean,width,alpha,n)
    gaus = root.RooGaussian(name+"_Gaus",name+"_Gaus",dimuonMass,mean,width2)
    pdfMmumu = root.RooAddPdf(convertSigName(name,energyStr),convertSigName(name,energyStr),cb,gaus,mix)
    
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)
    #alpha.setVal(1.4)
    #n.setVal(2.2)

    if forceMean > 0.0:
        mean.setVal(forceMean)

    ## Error time

    rooParamList = [mean,width,width2,alpha,n,mix]
    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      #workspaceImportFn(rooDataset)
      workspaceImportFn(fr)

    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName(name+"_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    now = datetime.datetime.now().isoformat()
#    canvas.SaveAs("debug_"+name+now+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(dimuonMass),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigDG(name,rooDataset,dimuonMass,minMass,maxMass,workspaceImportFn,channelName,energyStr,forceMean=-1.,sigInject=0):

    debug = ""
    debug += "### makePDFSigDG: "+channelName+": "+name+"\n"
    debug += "#    {0:.2f} < {1} < {2:.2f}\n".format(minMass,dimuonMass.GetName(),maxMass)
    debug += "#    {0:.2f} Events in RooDataSet\n".format(rooDataset.sumEntries())

    rooParamList = []

    meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                             channelName+"_"+name+"_MeanG1", 
                             123.,100.,150.)
    meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                             channelName+"_"+name+"_MeanG2", 
                             125.,100.,150.)
    
    widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                             channelName+"_"+name+"_WidthG1", 
                             5.,2.,10.)
    widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                              channelName+"_"+name+"_WidthG2", 
                              1.,0.5,5.)
    
    mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                            channelName+"_"+name+"_mixGG", 
                            0.5,0.,1.)
    gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                             channelName+"_"+name+"_gaus1",
                             dimuonMass,meanG1,widthG1)
    gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                             channelName+"_"+name+"_gaus2",
                             dimuonMass,meanG2,widthG2)
    pdfMmumu = root.RooAddPdf(convertSigName(name,energyStr),
                              name,
                              gaus1,gaus2,mixGG)
    rooParamList += [meanG1,meanG2,widthG1,widthG2,mixGG]
    
    fr = pdfMmumu.fitTo(rooDataset,root.RooFit.SumW2Error(False),PRINTLEVEL,root.RooFit.Save(True),root.RooFit.Range("signalfit"))
    fr.SetName(name+"_fitResult")

    #mean.setVal(125.)
    #width.setVal(2.)
    #alpha.setVal(1.4)
    #n.setVal(2.2)

    if forceMean > 0.0:
        meanG2.setVal(forceMean)

    ## Error time

    paramList = [Param(i.GetName(),i.getVal(),i.getError(),i.getError()) for i in rooParamList]
    print "makePDFSigDG:", paramList
    for i in rooParamList:
       i.setConstant(True)

    if workspaceImportFn != None:
      workspaceImportFn(pdfMmumu)
      #workspaceImportFn(rooDataset)
      workspaceImportFn(fr)

    ## Debug Time
#    frame = dimuonMass.frame()
#    frame.SetName(name+"_Plot")
#    rooDataset.plotOn(frame)
#    pdfMmumu.plotOn(frame)
#    canvas = root.TCanvas()
#    frame.Draw()
#    now = datetime.datetime.now().isoformat()
#    canvas.SaveAs("debug_"+name+now+".png")

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())

    sigInjectDataset = None
    if sigInject > 0:
      sigInjectDataset = pdfMmumu.generate(root.RooArgSet(dimuonMass),sigInject,False)

    return paramList, debug, sigInjectDataset

def makePDFSigNew(channelName,name,dimuonMass,mass,workspaceImportFn,useDG=True):
#def makePDFSigNew(channelName,name,dimuonMass,mass,workspaceImportFn,useDG=False):

    debug = ""
    debug += "### makePDFSigNew: "+channelName+": "+name+"\n"
    debug += "#    useDG: {0}\n".format(useDG)
    matchChannel = re.match(r"(.*)([\d]TeV)",channelName)
    assert(matchChannel)
    category = matchChannel.group(1)
    energy = matchChannel.group(2)
    #massStr = "{0:.0f}".format(mass)
    #if mass % 1 > 0.0:
    massStr = "{0:.1f}".format(mass)
    fitTypeString = "CBG"
    if useDG:
      fitTypeString = "DG"
    #params = sigFits.parameters[energy][massStr][category][fitTypeString]
    prodMode = 'gg'
    if ('Jet2' in category):
      prodMode = 'vbf'

    #print ' ####### prodMode = %s' % prodMode 
    #print ' ####### energy   = %s' % energy
    #print ' ####### category = %s' % category
    #print ' ####### massStr  = %s' % massStr
    sigFits = signalPars('fitresults',
                         prodMode,
                         energy,
                         category)
    
    par_meanG1, par_widthG1, par_meanG2, par_widthG2, par_mixGG = sigFits.getPars()
    if (par_widthG2>par_widthG1):
      par_meanG2, par_widthG2, par_meanG1, par_widthG1, par_mixGG = sigFits.getPars()
      


    #sigFitsCB = signalParsCB('fitresults',
    #                         prodMode,
    #                         energy,
    #                         category)
    #
    #par_mean, par_width, par_alpha, par_n = sigFitsCB.getPars()


    rooParamList = []

    # deal with shape systematic directly here and NOT in
    # nuisanceMap which for non-shape systematics
    paramList = [] # this is the list of shapes parameters to vary
    muonres = {}
    muonres["7TeV"] = 0.03 # this is a percentage
    muonres["8TeV"] = 0.02 # this is a percentage  

    muonscale = {}
    muonscale["7TeV"] = 0.002 # this is a percentage
    muonscale["8TeV"] = 0.002 # this is a percentage  

    if useDG:
      #if len(params.keys())>4:
      # define the Double Gaussian
      meanG1 = root.RooRealVar(channelName+"_"+name+"_MeanG1",
                               channelName+"_"+name+"_MeanG1", 
                               par_meanG1[massStr])
      #params['meanG1'])
      meanG2 = root.RooRealVar(channelName+"_"+name+"_MeanG2",
                               channelName+"_"+name+"_MeanG2", 
                               par_meanG2[massStr])
                               ##par_meanG1[massStr])
      #params['meanG2'])
      
      widthG1 = root.RooRealVar(channelName+"_"+name+"_WidthG1",
                                channelName+"_"+name+"_WidthG1", 
                                par_widthG1[massStr])
      #params['widthG1'])
      widthG2 = root.RooRealVar(channelName+"_"+name+"_WidthG2",
                                channelName+"_"+name+"_WidthG2", 
                                par_widthG2[massStr])
      #params['widthG2'])
        
      mixGG = root.RooRealVar(channelName+"_"+name+"_mixGG",
                              channelName+"_"+name+"_mixGG", 
                              par_mixGG[massStr])
      #params['mixGG'])
      gaus1 = root.RooGaussian(channelName+"_"+name+"_gaus1",
                               channelName+"_"+name+"_gaus1",
                               dimuonMass,meanG1,widthG1)
      gaus2 = root.RooGaussian(channelName+"_"+name+"_gaus2",
                               channelName+"_"+name+"_gaus2",
                               dimuonMass,meanG2,widthG2)
      pdfMmumu = root.RooAddPdf(convertSigName(name,energy),
                                name,
                                gaus1,gaus2,mixGG)

      workspaceImportFn(pdfMmumu)
      rooParamList = [meanG1,meanG2,widthG1,widthG2,mixGG]

      # adding a parameter to this list makes it float in the combination tool
      # for shape systematic uncertainties evaluation
      paramList.append(Param(meanG2.GetName(),meanG2.getVal(),meanG2.getVal()*muonscale[energy],meanG2.getVal()*muonscale[energy]))
      paramList.append(Param(widthG2.GetName(),widthG2.getVal(),widthG2.getVal()*muonres[energy],widthG2.getVal()*muonres[energy]))

    else:
      mean = root.RooRealVar(channelName+"_"+name+"_Mean",
                             channelName+"_"+name+"_Mean",
                             par_mean[massStr])
                             #params['mean'])
      width = root.RooRealVar(channelName+"_"+name+"_Width",
                              channelName+"_"+name+"_Width",
                              par_width[massStr])
                              #params['width1'])
      #width2 = root.RooRealVar(channelName+"_"+name+"_Width2",
      #                         channelName+"_"+name+"_Width2",
      #                         params['width2'])
      alpha = root.RooRealVar(channelName+"_"+name+"_Alpha",
                              channelName+"_"+name+"_Alpha",
                              par_alpha[massStr])
                              #params['Alpha'])
      n = root.RooRealVar(channelName+"_"+name+"_n",
                          channelName+"_"+name+"_n",
                          par_n[massStr])
                          #params['n'])
      #mix = root.RooRealVar(channelName+"_"+name+"_mix",
      #                         channelName+"_"+name+"_mix",
      #                         params['mix'])
      #cb = root.RooCBShape(name+"_CB",name+"_CB",dimuonMass,mean,width,alpha,n)
      #gaus = root.RooGaussian(name+"_Gaus",name+"_Gaus",dimuonMass,mean,width2)
      #pdfMmumu = root.RooAddPdf(convertSigName(name,energy),name,cb,gaus,mix)
      pdfMmumu = root.RooCBShape(convertSigName(name,energy),name,dimuonMass,mean,width,alpha,n)
      #pdfMmumu = root.RooAddPdf(convertSigName(name),name,cb,gaus,mix)
      workspaceImportFn(pdfMmumu)
      rooParamList += [mean,width,alpha,n]

    for i in rooParamList:
      debug += "#    {0:<35}: {1:<8.3f} +/- {2:<8.3f}\n".format(i.GetName(),i.getVal(),i.getError())
      i.setConstant(True)
      
    
    #print "##########################"
    #print paramList
    #print "##########################"
    return paramList, debug


#makePDFSig = makePDFSigCBPlusGaus
makePDFSig = makePDFSigDG
## makePDFBak = makePDFBakOld
#makePDFBak = makePDFBakExpMOverSq
#makePDFBak = makePDFBakBernsteinProd
makePDFBak = makePDFBakMSSM
#makePDFBak = makePDFBakVoigtPMm2
#makePDFBak = makePDFBakVoigtPExpMm2

###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=False,sigInject=0.0,sigInjectMass=125.0,energyStr="8TeV",cutString=""):
    if analysis[0] in "0123456789-.*+/":
       print("Error: Analysis: analysis name '{0}' begins with a forbidden character".format(analysis))
       sys.exit(1)
    getCutHist = getattr(self,"getCutHist")
    doSigInject = getattr(self,"doSigInject")
    self.treename = "outtree"
    self.directory = directory
    self.weightName = "puWeight"
    self.fileList = []
    self.treeList = []
    self.origHistList = []
    self.binSize = 0.5

    higgsPeakMean = args.higgsMass - 0.3
    self.higgsMass = 125.
    if args.higgsMass > 0.0:
      self.higgsMass = args.higgsMass
    higgsMassStr = "{0:.0f}".format(self.higgsMass)
    if self.higgsMass % 1 > 0.0:
      higgsMassStr = "{0:.1f}".format(self.higgsMass)
    self.higgsMassStr = higgsMassStr

    self.energyStr = energyStr
    self.lumi = lumi
    self.sigNames = signalNames
    self.bakNames = backgroundNames
    self.datNames = dataNames
    self.controlRegionVeryLow = controlRegionVeryLow
    self.controlRegionLow = controlRegionLow
    self.controlRegionHigh = controlRegionHigh
    self.analysis = analysis
    self.params = []
    self.debug = "#******************************************************#\n"*3
    self.debug += "#Nominal Higgs Mass: "+str(self.higgsMass) +"\n"
    self.debug += "#Peak Centered at: "+str(higgsPeakMean) +"\n"
    self.debug += "#Using lumi: "+str(lumi)+"\n"
    self.debug += "########################\n"
    self.debug += "#input directory name: "+str(directory)+"\n"
    self.debug += "#Data names:\n"
    for i in dataNames:
      self.debug += "#  "+str(i)+"\n"
    self.debug += "#signal names:\n"
    for i in signalNames:
      self.debug += "#  "+str(i)+"\n"
    self.debug += "#background names:\n"
    for i in backgroundNames:
      self.debug += "#  "+str(i)+"\n"
    self.debug += "########################\n"
    self.debug += "# USEGPANNA: "+str(USEGPANNA)+"\n"
    self.debug += "# SIGNALFIT: "+str(SIGNALFIT)+"\n"
    self.debug += "# USETREES: "+str(USETREES)+"\n"
    self.debug += "########################\n"

    self.workspace = root.RooWorkspace(analysis+energyStr)
    self.workspaceName = analysis+energyStr

    #wImport = getattr(self.workspace,"import")
    def wImport(var):
      if type(var) == root.RooFitResult:
        return
      if type(var) == root.RooDataSet or type(var) == root.RooDataHist:
        getattr(self.workspace,"import")(var,root.RooFit.RenameVariable("dimuonMass","dimuonMass_CMShmm"))
      else:
        getattr(self.workspace,"import")(var,root.RooFit.RenameAllNodes(self.workspaceName+"_CMShmm"),root.RooFit.RenameAllVariables("CMShmm"))

    self.sigInjectWorkspaces = []

    maxMass = controlRegionHigh[1]
    minMass = controlRegionLow[0]
    self.minMass = minMass
    self.maxMass = maxMass
    #dimuonMass = root.RooRealVar("dimuonMass","dimuonMass",minMass,maxMass)
    dimuonMass = root.RooRealVar("dimuonMass","dimuonMass",minMass,maxMass)
    dimuonMass.setRange("exprange",120.,controlRegionHigh[1])
    dimuonMass.setRange("whole",controlRegionLow[0],controlRegionHigh[1])
    dimuonMass.setRange("low",controlRegionLow[0],controlRegionLow[1])
    dimuonMass.setRange("high",controlRegionHigh[0],controlRegionHigh[1])
    dimuonMass.setRange("signal",controlRegionLow[1],controlRegionHigh[0])
    dimuonMass.setRange("signalfit",SIGNALFIT[0],SIGNALFIT[1])
    self.dimuonMass = dimuonMass

    self.rooPUWeight = root.RooRealVar("puWeight","PU Weight",0.0,200.)
    self.observablesForRooDataset = root.RooArgSet(
                            self.dimuonMass,
                            self.rooPUWeight,
                        )

    # Hack to Make makePDFBakOld work
    self.minMassZ = 88.
    self.maxMassZ = 94.
    dimuonMassZ = root.RooRealVar("dimuonMass","dimuonMass",self.minMassZ,self.maxMassZ)
    self.dimuonMassZ = dimuonMassZ
    self.observablesForRooDatasetZ = root.RooArgSet(
                            dimuonMassZ,
                            self.rooPUWeight,
                        )

    # Dealing with cut string
    KDString = '(sigMEPdf*{0}/(bakMEPdf*{1}+sigMEPdf*{0}))'.format(
                   MENormDict[energyStr]["sigMEPdf"],MENormDict[energyStr]["bakMEPdf"]
                )
    self.fullCutString = treeCut(analysis,cutString,eventWeights=False,muonRequirements=True,KDString=KDString)
    self.debug += "#  Full Cut String:\n"
    self.debug += "#  {0}\n".format(self.fullCutString)

    self.sigHistsRaw = []
    if not USEGPANNA:
      for name in signalNames:
        tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
        self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
        self.sigHistsRaw.append(tmpH)

    self.bakHistsRaw = []
    self.bakHistsRawZ = []
    for name in backgroundNames:
      tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
      self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
      self.bakHistsRaw.append(tmpH)
      tmpH = self.getRooDataSample(name,self.observablesForRooDatasetZ,True)
      self.bakHistsRawZ.append(tmpH)

    self.datHists = []
    self.datHistsZ = []
    for name in dataNames:
      tmpH = None
      if "Jet2CutsVBFPass" in self.analysis:
        tmpH = self.getRooDataSet(name,self.observablesForRooDataset)
      else:
        tmpH = self.getRooDataSample(name,self.observablesForRooDataset)
      self.debug += "#  N Events in {0}: {1}\n".format(name,tmpH.sumEntries())
      self.datHists.append(tmpH)
      tmpH = None
      if "Jet2CutsVBFPass" in self.analysis:
        tmpH = self.getRooDataSet(name,self.observablesForRooDatasetZ,True)
      else:
        tmpH = self.getRooDataSample(name,self.observablesForRooDatasetZ,True)
      self.datHistsZ.append(tmpH)

    effMap = {}
    xsecMap = {}
    lowBin = 0
    massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    self.massBounds = massBounds

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    self.bakHists = []
    self.bakHistTotal = None
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      sys.exit(1)
      counts = h.sumEntries()
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      self.xsecBakTotal += xs
      self.xsecBakList.append(xs)
      self.effBakList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])
      self.bakHists.append(h)
      if self.bakHistTotal == None:
        self.bakHistTotal = h.Clone("bak")
      else:
        self.bakHistTotal.Add(h)

    #self.bakHistTotal.Scale(lumi)
    #self.bakHistTotalReal = self.bakHistTotal.Clone("data_obs_"+self.workspaceName+"_CMShmm")

    self.dataCountsTotal = None
    self.datHistTotal = None
    for h,name in zip(self.datHists,dataNames):
      counts = h.sumEntries()
      if self.dataCountsTotal == None:
        self.dataCountsTotal = counts
      else:
        self.dataCountsTotal += counts
      if self.datHistTotal == None:
        self.datHistTotal = h.Clone("data_obs_"+self.workspaceName+"_CMShmm")
      else:
        assert(type(self.datHistTotal)==type(h))
        if type(h) == root.RooDataHist:
          self.datHistTotal.add(h)
        else:
          self.datHistTotal.append(h)
    assert(self.dataCountsTotal == self.datHistTotal.sumEntries())

    self.dataCountsTotalZ = None
    self.datHistTotalZ = None
    for h,name in zip(self.datHistsZ,dataNames):
      counts = h.sumEntries()
      if self.dataCountsTotalZ == None:
        self.dataCountsTotalZ = counts
      else:
        self.dataCountsTotalZ += counts
      if self.datHistTotalZ == None:
        self.datHistTotalZ = h.Clone("data_obs_"+self.workspaceName+"_Z_CMShmm")
      else:
        assert(type(self.datHistTotalZ)==type(h))
        if type(h)==root.RooDataHist:
          self.datHistTotalZ.add(h)
        else:
          self.datHistTotalZ.append(h)
    assert(self.dataCountsTotalZ == self.datHistTotalZ.sumEntries())

    histToUseForBak = self.datHistTotal
    histToUseForBakZ = self.datHistTotalZ
    if histToUseForBak is None:
        print("Background MC rescaling RooDataSet to lumi doesn't work...\n  you must use real data...Exiting")
        sys.exit(1)
        doSigInject(self.bakHistTotal,sigInject,sigInjectMass)
        histToUseForBak = self.bakHistTotal
    else:
        self.dataCountsTotal += doSigInject(self.datHistTotal,sigInject,sigInjectMass)

    tmpBakParams, self.bakNormTup, tmpDebug, tmpOrder = makePDFBak(
                                    analysis+energyStr,histToUseForBak,
                                    dimuonMass, minMass,maxMass,wImport,
                                    dimuonMassZ,histToUseForBakZ
                                       )
    self.params.extend(tmpBakParams)
    self.countsBakTotal = self.bakNormTup[0]*self.bakNormTup[1]
    self.debug += tmpDebug
    
    self.sigParamList = []
    if USEGPANNA:
      sigParams, sigDebug =  makePDFSigNew(analysis+energyStr,"sigPdf",dimuonMass,self.higgsMass,wImport)
      self.sigParamList.append(sigParams)
      self.debug += sigDebug
      self.sigPdf = self.workspace.pdf("sigPdf_hmm"+energyStr+"_"+analysis+energyStr+"_CMShmm")
    else: # Use ggH for non-VBF channels, and VBF shape for VBF
      sigNameToUse = None
      sigHistCounts = []
      for name, hist in zip(signalNames,self.sigHistsRaw):
         sigHistCounts += [hist.sumEntries()]
      maxCounts = max(sigHistCounts)
      maxIndex = sigHistCounts.index(maxCounts)
      self.debug += "# Signal Histogram for Shape: {0} Counts: {1:.1f}\n".format(signalNames[maxIndex],maxCounts)
      sigHistRawToUse = self.sigHistsRaw[maxIndex]
      sigParams, sigDebug, tmpDS = makePDFSig("sigPdf",sigHistRawToUse,dimuonMass,minMass,maxMass,wImport,analysis+energyStr,energyStr)
      self.sigPdf = self.workspace.pdf("sigPdf_hmm"+energyStr+"_"+analysis+energyStr+"_CMShmm")
      self.sigParamList.append(sigParams)
      self.debug += sigDebug
          
    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    if USEGPANNA:
      for name in signalNames:
        nameMatch = re.match(r"(.+)Hmumu[\d.]+_[\d]+TeV",name)
        assert(nameMatch)
        prodMode = nameMatch.group(1)
        #eff, effErr = effReader(self.energyStr,prodMode,analysis,higgsMassStr)
        #print ' ******* EFFICIENCY READER ******* ' 
        #print ' ####### prodMode = %s' % prodMode 
        #print ' ####### energy   = %s' % self.energyStr
        #print ' ####### category = %s' % analysis
        #print ' ####### mass     = {0:.1f}'.format(self.higgsMass) 
        
        effReader = effReaderFromFile('fitresults',
                                      prodMode,
                                      self.energyStr,
                                      analysis)
        efficiencies = effReader.getEff()
        eff = efficiencies[ '{0:.1f}'.format(self.higgsMass) ]
        effErr = 0 #error is not used (syst. dominated)

        #print ' ####### higgsMassStr = %s' % higgsMassStr 
        tmpName = name.replace("125",higgsMassStr)
        #print ' #######  tmpName = %s' % tmpName
        xs = eff*xsec[tmpName]
        self.xsecSigTotal += xs
        self.xsecSigList.append(xs)
        self.effSigList.append(eff)
    else:
      for h,name in zip(self.sigHistsRaw,signalNames):
        counts = h.sumEntries()
        eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
        xs = eff*xsec[name]
        if args.higgsMass > 0.0:
          tmpName = name.replace("125",higgsMassStr)
          xs = eff*xsec[tmpName]
        self.xsecSigTotal += xs
        self.xsecSigList.append(xs)
        self.effSigList.append(eff)

    self.countsSigTotal = self.xsecSigTotal*lumi
    #self.countsBakTotal = self.xsecBakTotal*lumi
    self.countsSigList = [x*lumi for x in self.xsecSigList]
    self.countsBakList = [x*lumi for x in self.xsecBakList]

    if self.dataCountsTotal is None:
      assert(False)
      bakDataTH1 = self.bakHistTotal.Clone("bak_Template"+str(random.randint(0,10000)))
      bakDataTH1 = shrinkTH1(bakDataTH1,minMass,maxMass)
      for i in range(bakDataTH1.GetNbinsX()+2):
        binCenter = bakDataTH1.GetXaxis().GetBinCenter(i)
        contentInt = int("{0:.0f}".format(bakDataTH1.GetBinContent(i)))
        bakDataTH1.SetBinContent(i,0.0)
        bakDataTH1.SetBinError(i,0.0)
        for j in range(contentInt):
          bakDataTH1.Fill(binCenter)
      self.dataCountsTotal = int(bakDataTH1.Integral())
      obsData = root.RooDataHist("data_obs_"+self.workspaceName+"_CMShmm","MC Full-Sim Data",root.RooArgList(dimuonMass),bakDataTH1)
      print "counts: {} obsData: {}".format(self.dataCountsTotal,obsData.sumEntries())
      wImport(obsData)
      self.bakHistTotal.SetName("data_obs_"+self.workspaceName+"_CMShmm_"+analysis)
      #degubf = root.TFile("debug.root","recreate")
      #self.bakHistTotal.Write()
      #degubf.Close()
    elif toyData:
      bakPDF = self.workspace.pdf("bak")
      sigPDFList = [self.workspace.pdf(i) for i in signalNames]
      toyDataset = bakPDF.generate(root.RooArgSet(dimuonMass),int(self.dataCountsTotal))
      doSigInject(toyDataset,sigInject,sigInjectMass)
      #toyDataHist = toyDataset.binnedClone("data_obs_"+self.workspaceName+"_CMShmm","Toy Data")
      self.dataCountsTotal = int(toyDataset.sumEntries())
      wImport(toyDataset)
    else:
      #realDataHist = root.RooDataHist("data_obs_"+self.workspaceName+"_CMShmm","Real Observed Data",root.RooArgList(dimuonMass),self.datHistTotal)
      #wImport(realDataHist)
      #realDataHistNotVeryLow = realDataHist.reduce(root.RooFit.CutRange("low,signal,high"))
      #wImport(realDataHistNotVeryLow)
      self.datHistTotal.SetTitle("Real Observed Data")
      wImport(self.datHistTotal)

    self.debug += "#  Signal Efficiencies: \n"
    for sigName,sigEff in zip(signalNames,self.effSigList):
      self.debug += "#    {0}: {1:.2g}\n".format(sigName,sigEff)

    ### Background Parameterization Uncertainty Part
    if BAKPARAMUNC:
      self.bakParamUncPdfName = "bakParamUnc"
      self.bakParamUncPdfName = "bakParamUnc"+"_hmm"+energyStr
      self.bakParamUncCounts = BakParameterizationUncDict[energyStr][analysis]
      self.bakParamUncPdf = self.sigPdf.Clone(self.bakParamUncPdfName+"_"+self.workspaceName+"_CMShmm")
      getattr(self.workspace,"import")(self.bakParamUncPdf)

      self.bakParamUncNormVar = root.RooRealVar(self.bakParamUncPdfName+"_"+self.workspaceName+"_CMShmm_norm","bakParamUnc_hmm"+energyStr+"_"+self.workspaceName+"_CMShmm_norm",0.,-3.*self.bakParamUncCounts,3.*self.bakParamUncCounts)
      self.bakParamUncNormVar.setConstant(True)
      getattr(self.workspace,"import")(self.bakParamUncNormVar)


  def getBakParamUncName(self,full=False,pdfName=False):
    if pdfName:
      return self.bakParamUncPdfName+"_"+self.workspaceName+"_CMShmm"
    if full:
      # old style example: n_exp_binJets01FailPtG10BO8TeV_proc_bakParamUnc param 0.0 50.0
      #return "n_exp_bin"+self.analysis+self.energyStr+"_proc_"+self.bakParamUncPdfName
      # Mingshui says that just the rooPDF name + "_norm" should do it
      # new style example: shapeBkg_bakParamUnc_hmm8TeV_Jet2CutsVBFPass8TeV__norm
      return "shapeBkg_"+self.bakParamUncPdfName+"_"+self.workspaceName+"__norm"
    else:
      return self.bakParamUncPdfName
  def getBakParamUncCounts(self):
    return self.bakParamUncCounts

  def getCutHist(self,inHist,outHist,cut,cutAbove=False):
    if not inHist.InheritsFrom("TH1"):
      print("Error: Analysis.getCutHist(): inHist is not a hist exiting.")
      sys.exit(1)
    if not inHist.InheritsFrom("TH2"):
      print("Error: Analysis.getCutHist(): inHist is not 2D!! exiting.")
      sys.exit(1)
    if not outHist.InheritsFrom("TH1"):
      print("Error: Analysis.getCutHist(): outHist is not a hist exiting.")
      sys.exit(1)
    if outHist.InheritsFrom("TH2"):
      print("Error: Analysis.getCutHist(): outHist is 2D!! exiting.")
      sys.exit(1)
    cutBin = inHist.GetYaxis().FindBin(cut)
    nBinsY = inHist.GetYaxis().GetNbins()
    nBinsX = inHist.GetXaxis().GetNbins()
    #print("max bin: {0}, cutBin: {1}".format(nBinsY,cutBin))
    for iX in range(0,nBinsX+2):
      mySum = 0.
      if keepBelow:
        for iY in range(0,cutBin):
          mySum += inHist.GetBinContent(iX,iY)
      else:
        for iY in range(cutBin,nBinsY+2):
          mySum += inHist.GetBinContent(iX,iY)
      outHist.SetBinContent(iX,mySum)

  def doSigInject(self,dataHist,sigStrength,sigMass):
    if sigStrength <= 0.0:
        return 0
    dimuonMass = self.dimuonMass
    w = root.RooWorkspace("tmpWorkspace"+str(random.randint(0,10000)))
    self.sigInjectWorkspaces.append(w)
    wImport = getattr(w,"import")
    countsList = []
    for h,name in zip(self.sigHistsRaw,self.sigNames):
      counts = getIntegralAll(h,boundaries=self.massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      if sigMass > 0.0:
        prec = "0"
        if sigMass % 1 > 0.0:
          prec = "1"
        higgsMassString =("{0:."+prec+"f}").format(sigMass)
        tmpName = name.replace("125",higgsMassString)
        xs = eff*xsec[tmpName]
      countsList.append(int(xs*self.lumi*sigStrength))
    for name, hist, counts in zip(self.sigNames,self.sigHistsRaw,countsList):
      if counts == 0:
        continue
      tmpName = name+"TmpSigInject"+str(random.randint(0,10000))
      sigParams, sigDebug, rooData = makePDFSig(tmpName,hist,dimuonMass,sigMass-20,sigMass+20,wImport,tmpName,energyStr,forceMean=sigMass-0.3,sigInject=counts)
      if dataHist.InheritsFrom("RooAbsData"):
        dataHist.append(rooData)
      else:
        tmpHist = dataHist.Clone("sigInjectHistTmp"+str(random.randint(0,10000)))
        tmpHist.Reset()
        rooData.fillHistogram(tmpHist,root.RooArgList(dimuonMass))
        dataHist.Add(tmpHist)
    return sum(countsList)

  def getRooDataSet(self,name,observables,aroundZ=False):
    tmpFLoc = self.directory+name+".root"
    tmpF = root.TFile(tmpFLoc)
    minMass = self.minMass
    maxMass = self.maxMass
    dimuonMass = self.dimuonMass
    if aroundZ:
      minMass = self.minMassZ
      maxMass = self.maxMassZ
      dimuonMass = self.dimuonMassZ

    thetreename  = self.treename
    thetreename += self.analysis

    tmpTree = tmpF.Get(self.treename + self.analysis)
    tmpTree.SetCacheSize(10000000);
    tmpTree.AddBranchToCache("*");

    tmpDataSet = root.RooDataSet(name,name,tmpTree,root.RooArgSet(dimuonMass))
    return tmpDataSet

  def getRooDataSample(self,name,observables,aroundZ=False):
    name = name.replace("whH","wH")
    name = name.replace("zhH","zH")
    tmpFLoc = self.directory+name+".root"
    tmpF = root.TFile(tmpFLoc)
    minMass = self.minMass
    maxMass = self.maxMass
    dimuonMass = self.dimuonMass
    if aroundZ:
      minMass = self.minMassZ
      maxMass = self.maxMassZ
      dimuonMass = self.dimuonMassZ
    if USETREES:
      tmpTree = tmpF.Get(self.treename)
      tmpTree.SetCacheSize(10000000);
      tmpTree.AddBranchToCache("*");
      histName = "hist{0:f}".format(time.time()).replace('.','')
      nBins = int((maxMass-minMass)/self.binSize)
      tmpHist = root.TH1F(histName,"",nBins,minMass,maxMass)
      drawString = "dimuonMass >> {0}".format(histName)
      tmpTree.Draw(drawString,self.fullCutString)
      self.origHistList.append(tmpHist)
      tmpH = root.RooDataHist(name,name,root.RooArgList(dimuonMass),tmpHist)
      return tmpH
    else:
      histName = self.analysis+"/"+HISTNAME
      tmpHist = tmpF.Get(histName)
      tmpHist.GetNbinsX()
      tmpHist = shrinkTH1(tmpHist,minMass,maxMass,True)
      tmpHist.GetNbinsX()
      self.origHistList.append(tmpHist)
      tmpH = root.RooDataHist(name,name,root.RooArgList(dimuonMass),tmpHist)
      return tmpH

  def getSigEff(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.effSigList[i]
    return result
  def getSigXSec(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.xsecSigList[i]
    return result
  def getSigXSecTotal(self):
    return self.xsecSigTotal
  def getBakXSecTotal(self):
    return self.xsecBakTotal
  def getBakXSec(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.xsecBakList[i]
    return result
  def getDataTotal(self):
    return self.dataCountsTotal
  def getParamList(self):
    return self.params
  def getSigCounts(self,name):
    result = -1.0
    if self.sigNames.count(name)>0:
        i = self.sigNames.index(name)
        result = self.countsSigList[i]
    return result
  def getBakCounts(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.countsBakList[i]
    return result
  def getSigCountsTotal(self):
    return self.countsSigTotal
  def getBakCountsTotal(self):
    return self.countsBakTotal
  def getEnergy(self):
    return self.energyStr

###################################################################################

class DataCardMaker:
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,outfilename,lumi,nuisanceMap=None,controlRegionLow=[110.,115],controlRegionHigh=[135,150],controlRegionVeryLow=[80.,110.],sigInject=0.0,sigInjectMass=125.0,toyData=False,energyStr="8TeV",cutString=""):

    ########################
    ## Setup

    if len(analysisNames)==0:
      print("Error: DataCardMaker: len(analysisNames)==0")
      return

    channels = []
    self.is2D = False

    if type(energyStr) == list:
      for analysis in analysisNames:
        cutStringTmp = cutString
        if type(analysis)==list:
          if len(analysis)>1:
            if len(cutStringTmp)>0:
              cutStringTmp += " && "
            cutStringTmp += analysis[1]
          analysis = analysis[0]
        for es,sn,bn,dn,lu in zip(energyStr,signalNames,backgroundNames,dataNames,lumi):
          lu *= 1000.0
          tmp = Analysis(directory,sn,bn,dn,analysis,lu,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,energyStr=es,cutString=cutStringTmp)
          channels.append(tmp)
      self.channels = channels
    else:
      lumi *= 1000.0
      for analysis in analysisNames:
        cutStringTmp = cutString
        if type(analysis)==list:
          if len(analysis)>1:
            if len(cutStringTmp)>0:
              cutStringTmp += " && "
            cutStringTmp += analysis[1]
          analysis = analysis[0]
        tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,lumi,controlRegionVeryLow,controlRegionLow,controlRegionHigh,toyData=toyData,sigInject=sigInject,sigInjectMass=sigInjectMass,energyStr=energyStr,cutString=cutStringTmp)
        channels.append(tmp)
      self.channels = channels

    self.channelNames = [i.workspaceName for i in channels]

    self.nuisance = nuisanceMap
    nuisance = self.nuisance

    self.largestChannelName = 0
    for name in self.channelNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    for channel in channels:
      for name in channel.sigNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
      for name in channel.bakNames:
        if len(name)>self.largestChannelName:
          self.largestChannelName = len(name)
    if self.largestChannelName < 8:
        self.largestChannelName = 8
    self.largestChannelName += 2

    if self.channelNames.count("")>0:
      i = self.channelNames.index("")
      self.channelNames[i] = "Inc"

    self.controlRegionHigh = controlRegionHigh
    self.controlRegionLow = controlRegionLow
    self.controlRegionVeryLow = controlRegionVeryLow
    self.toyData = toyData

    ### ROOT Part
    ##########################################################
    outRootFilename = re.sub(r"\.txt",r".root",outfilename)
    outRootFile = root.TFile(outRootFilename, "RECREATE")
    outRootFile.cd()

    rootDebugString = ""

    for channel in self.channels:
        print(channel.workspace.data("data_obs_"+channel.workspaceName+"_CMShmm").GetTitle())
        rootDebugString += "# "+channel.workspace.data("data_obs_"+channel.workspaceName+"_CMShmm").GetTitle()+"\n"
        channel.workspace.Write()

    outRootFile.Close()

    ### Text Part
    ##########################################################

    print("Writing Card: {0} & {1}".format(outfilename,outRootFilename))
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu shape combine datacard produced by makeCards.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format("*"))
    outfile.write("------------\n")
    ## read in shapes
    for channel,channelName in zip(self.channels,self.channelNames):
      for sigName in channel.sigNames:
        outfile.write("shapes {0} {1} {2} {1}:{3}_{1}_CMShmm\n".format(convertSigName(sigName,channel.getEnergy()),channelName, os.path.basename(outRootFilename),"sigPdf_hmm"+channel.energyStr))
      outfile.write("shapes {0} {1} {2} {1}:{0}_{1}_CMShmm\n".format("bak", channelName,os.path.basename(outRootFilename)))
      if BAKPARAMUNC:
        outfile.write("shapes {0} {1} {2} {1}:{3}\n".format(channel.getBakParamUncName(), channelName, os.path.basename(outRootFilename),channel.getBakParamUncName(pdfName=True)))
      outfile.write("shapes {0} {1} {2} {1}:{0}_{1}_CMShmm\n".format("data_obs", channelName, os.path.basename(outRootFilename)))
    outfile.write("------------\n")
    outfile.write("# Channels, observed N events:\n")
    # Make Channels String
    binFormatString = "bin           "
    observationFormatString = "observation  "
    binFormatList = self.channelNames
    observationFormatList = []
    iParam = 0
    for channel,channelName in zip(self.channels,self.channelNames):
      binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
      binFormatList.append(channelName)
      observationFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
      observedNumber = channel.getDataTotal()
      observationFormatList.append(observedNumber)
      #print("text Observed {0}: {1}".format(channelName,observedNumber))
      iParam += 1
    binFormatString+= "\n"
    observationFormatString+= "\n"
    outfile.write(binFormatString.format(*binFormatList))
    outfile.write(observationFormatString.format(*observationFormatList))
    outfile.write("------------\n")
    outfile.write("# Expected N events:\n")

    binFormatString = "bin           "
    proc1FormatString = "process       "
    proc2FormatString = "process       "
    rateFormatString = "rate          "
    binFormatList = []
    proc1FormatList = []
    proc2FormatList = []
    rateFormatList = []
    iParam = 0
    for channel,channelName in zip(self.channels,self.channelNames):
        ## Signal Time

        iProc = -len(channel.sigNames)+1
        for sigName in channel.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(convertSigName(sigName,channel.getEnergy()))
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigCounts(sigName)
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1

        ## Background Time

        binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        binFormatList.append(channelName)
  
        proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc1FormatList.append("bak")
  
        proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
        proc2FormatList.append(iProc)

        expNum = channel.getBakCountsTotal()
        decimals = ".4f"
        if expNum>1000.0:
          decimals = ".4e"
        rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
        rateFormatList.append(expNum)
    
        iParam += 1
        iProc += 1

        ## Bkg Parameterization Uncertainty Time
        if BAKPARAMUNC:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(channel.getBakParamUncName())
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)

          decimals = ".1f"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(1.0)  # Want in terms of n signal events NOT mu
    
          iParam += 1
          iProc += 1
    binFormatString+= "\n"
    proc1FormatString+= "\n"
    proc2FormatString+= "\n"
    rateFormatString+= "\n"
    outfile.write(binFormatString.format(*binFormatList))
    outfile.write(proc1FormatString.format(*proc1FormatList))
    outfile.write(proc2FormatString.format(*proc2FormatList))
    outfile.write(rateFormatString.format(*rateFormatList))
    outfile.write("------------\n")
    outfile.write("# Uncertainties:\n")

    ## lnN Uncertainties:
    ## Correlated between everything:
    #for nu in nuisance.keysEnergyCorr:
    #  formatString = "{0:<8} {1:^4} "
    #  formatList = [convertUncName(nu),"lnN"]
    #  iParam = 2
    #  for channel,channelName in zip(self.channels,self.channelNames):
    #      channelNameNoEnergy = re.sub(r"[\d]TeV$","",channelName)
    #      for sigName in channel.sigNames:
    #        formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #        value = nuisance(nu,sigName,channelNameNoEnergy,channel.higgsMassStr)
    #        if value == None:
    #          value = "-"
    #        else:
    #          value = abs(value)
    #        formatList.append(value)
    #        iParam += 1
    #      # For background
    #      formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #      value = "-"
    #      formatList.append(value)
    #      iParam += 1
    #      ## Bkg Parameterization Uncertainty Time
    #      if BAKPARAMUNC:
    #        formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #        value = "-"
    #        formatList.append(value)
    #        iParam += 1
    #  formatString += "\n"
    #  #print formatString
    #  #print formatList
    #  outfile.write(formatString.format(*formatList))
    ## Not Correlated between Energies:
    #for nu in nuisance.keysNotEnergyCorr:
    #  energyList = energyStr
    #  if type(energyStr) != list:
    #    energyList = [energyStr]
    #  for energy in energyList:
    #    formatString = "{0:<8} {1:^4} "
    #    formatList = [convertUncName(nu+"_"+energy),"lnN"]
    #    iParam = 2
    #    for channel,channelName in zip(self.channels,self.channelNames):
    #        channelEnergyMatch = re.search(r"[\d]TeV$",channelName)
    #        assert(channelEnergyMatch)
    #        channelEnergy = channelEnergyMatch.group(0)
    #        #print("channelEnergy: %s" % channelEnergy)
    #        channelNameNoEnergy = re.sub(r"[\d]TeV$","",channelName)
    #        for sigName in channel.sigNames:
    #          formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #          value = nuisance(nu,sigName,channelNameNoEnergy,channel.higgsMassStr)
    #          if value == None:
    #            value = "-"
    #          else:
    #            value = abs(value)
    #          if channelEnergy != energy:
    #            value = "-"
    #          formatList.append(value)
    #          iParam += 1
    #        # For background
    #        formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #        value = "-"
    #        formatList.append(value)
    #        iParam += 1
    #        ## Bkg Parameterization Uncertainty Time
    #        if BAKPARAMUNC:
    #          formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #          value = "-"
    #          formatList.append(value)
    #          iParam += 1
    #    formatString += "\n"
    #    #print formatString
    #    #print formatList
    #    outfile.write(formatString.format(*formatList))
    ## Not Correlated between Energies or categories:
    #for nu in nuisance.keysNotCatCorr:
    #  for channel,channelName in zip(self.channels,self.channelNames):
    #    for sigName in channel.sigNames:
    #      formatString = "{0:<8} {1:^4} "
    #      formatList = [convertUncName(nu+"_"+channelName+"_"+sigName),"lnN"]
    #      iParam = 2
    #      for channel2,channelName2 in zip(self.channels,self.channelNames):
    #        channelNameNoEnergy = re.sub(r"[\d]TeV$","",channelName)
    #        value = "-"
    #        for sigName2 in channel.sigNames:
    #          value = "-"
    #          formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #          if channelName == channelName2 and sigName == sigName2:
    #            value = nuisance(nu,sigName,channelNameNoEnergy,channel.higgsMassStr)
    #            if value == None:
    #              value = "-"
    #            else:
    #              value = abs(value)
    #          formatList.append(value)
    #          iParam += 1
    #        # For background
    #        formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #        value = "-"
    #        formatList.append(value)
    #        iParam += 1
    #        ## Bkg Parameterization Uncertainty Time
    #        if BAKPARAMUNC:
    #          formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
    #          value = "-"
    #          formatList.append(value)
    #          iParam += 1
    #      formatString += "\n"
    #      #print formatString
    #      #print formatList
    #      isEmpty = True
    #      for el in formatList[2:]:
    #        if el != "-":
    #          isEmpty = False
    #          break
    #      if isEmpty:
    #        continue
    #      outfile.write(formatString.format(*formatList))

    # Bak norm uncertainty
    for channel1,channel1Name in zip(self.channels,self.channelNames):
      formatString = "{0:<8} {1:^3}"
      formatList = [convertUncName("bkN"+channel1Name),"lnU"]
      iParam = 2
      for channel in self.channels:
          for sigName in channel.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            formatList.append(value)
            iParam += 1
          value = '-'
          tmpString = "{"+str(iParam)+":^"+str(self.largestChannelName)+"}"
          if channel == channel1:
            value = 2.
            tmpString = "{"+str(iParam)+":^"+str(self.largestChannelName)+".2f}"
          formatString += tmpString
          formatList.append(value)
          iParam += 1
          ## Bkg Parameterization Uncertainty Time
          if BAKPARAMUNC:
            value = '-'
            tmpString = "{"+str(iParam)+":^"+str(self.largestChannelName)+"}"
            formatString += tmpString
            formatList.append(value)
            iParam += 1
      formatString += "\n"
      outfile.write(formatString.format(*formatList))

      # Parameter (Shape) Uncertainties For Background
      for channel,channelName in zip(self.channels,self.channelNames):
        for nu in channel.params:
          if type(nu) == str: # Make named parameters flatParam
            nuisanceName = nu + "_CMShmm"
            formatString = "{0:<25} {1:<9}"
            formatList = [nuisanceName,"flatParam"]
            formatString += "\n"
            #print formatString
            #print formatList
            outfile.write(formatString.format(*formatList))

      ## Parameter (Shape) Uncertainties For Signal 
      #for channel,channelName in zip(self.channels,self.channelNames):
      #  for siglist in channel.sigParamList:
      #    for nu in siglist:
      #      #print nu
      #      nuisanceName = nu.name +"_CMShmm"
      #      formatString = "{0:<25} {1:<6} {2:<10.5g} {3:<10}"
      #      formatList = [nuisanceName,"param",nu.nominal,nu.getErrString()]
      #      formatString += "\n"
      #      #print formatString
      #      #print formatList
      #      outfile.write(formatString.format(*formatList))

      ## Parameter Uncertainty for background parameterization uncertainty
      #if BAKPARAMUNC:
      #  for channel,channelName in zip(self.channels,self.channelNames):
      #    nuisanceName = channel.getBakParamUncName(True)
      #    formatString = "{0:<25} {1:<6} {2:<10.1f} {3:<10.5g}"
      #    formatList = [nuisanceName,"param",0.0,channel.getBakParamUncCounts()]
      #    formatString += "\n"
      #    #print formatString
      #    #print formatList
      #    outfile.write(formatString.format(*formatList))

    #Debugging
    if not args.forCombGroup:
      outfile.write("#################################\n")
      for channel,channelName in zip(self.channels,self.channelNames):
          outfile.write("#\n")
          outfile.write("#info: channel {0}: \n".format(channelName))
          outfile.write(channel.debug)
      outfile.write("#\n#\n")
      outfile.write(rootDebugString)
      outfile.write("# Signal Injected: {0:.1f} X SM, mH = {1:.1f}\n".format(sigInject,sigInjectMass))

#      for channel,channelName in zip(self.channels,self.channelNames):
#          outfile.write("#\n")
#          outfile.write("#Efficiency: channel {0}: \n".format(channelName))
#          for sigName in channel.sigNames:
#            outfile.write("#  {0:>20}: {1:4f}\n".format(sigName,channel.getSigEff(sigName)))
    outfile.close()

class ThreadedCardMaker(myThread):
  def __init__(self,*args,**dictArgs):
    myThread.__init__(self)
    self.args = args
    self.dictArgs = dictArgs
    self.started = False
    #print("\nPos Args:\n {0} \nKeyword Args:\n {1}\n\n".format(self.args,self.dictArgs))
  def run(self):
    #try:
      self.started = True
      dataCardMassShape = None
      dataCardMassShape = DataCardMaker(*(self.args),**(self.dictArgs))
    #except Exception as e:
    #  print("Error: Exception: {0}".format(e))
    #  print("  Thread Arguments: {0}, {1}".format(self.args,self.dictArgs))

###################################################################################
###################################################################################
###################################################################################
###################################################################################
###################################################################################

if __name__ == "__main__":
  print "Started makeCards.py"
  print "Higgs mass set to: ", args.higgsMass
  root.gROOT.SetBatch(True)

  directory = getDataStage2Directory()
  #directory = "/data/uftrig01b/jhugon/hmumu/analysisV00-01-10/forGPReRecoMuScleFit/"
  #directory = "/cms/data/store/user/jhugon/hmumu/stage2/"
  #directory = "/afs/cern.ch/work/j/jhugon/public/hmumuNtuplesLevel2/unzipped/"
  outDir = "statsCards/"
  periods = ["7TeV","8TeV"]
  periods = ["8TeV"]
  #periods = ["7TeV"]
  categoriesAll = ["BB","BO","BE","OO","OE","EE"]
  categoriesFF = ["BB","BO","BE","OO","FF"]
  categoriesCC = ["BB","BO","CC","OE","EE"]
  categoriesCCFF = ["BB","BO","CC","FF"]
  categoriesFC = ["BB","BO","FC"]
  categoriesAllCCFF = ["BB","BO","BE","OO","OE","EE","CC","FF"]
  analyses = []

#  analyses += [["Inclusive",""]]
#  analyses += [["IncPresel",""]]
#  analyses += [["IncPreselPtG10",""]]
#  analyses += [["IncPreselPtG10"+i,""] for i in categoriesAll]
#
#  analyses += [["VBFPresel",""]]
#  analyses += [["VBFCutBased",""]]
#  analyses += [["VBFBDTCut",""]]

  combinations = []

#  combinations.append((
#        [["IncPresel"+x] for x in categoriesAll],"IncPreselCat"
#  ))
#  combinations.append((
#        [["IncPreselPtG10"+x] for x in categoriesAll],"IncCutCat"
#  ))
#  combinations.append((
#        [["VBFPresel"],["IncPresel"]],"CombPresel"
#  ))
#  combinations.append((
#        [["VBFCutBased"],["IncPreselPtG10"]],"CombCuts"
#  ))
#  combinations.append((
#        [["VBFCutBased"]]+[["IncPreselPtG10"+x] for x in categoriesAll],"CombCutsCat"
#  ))
#  combinations.append((
#        [["VBFBDTCut"],["IncPreselPtG10"]],"CombBDT"
#  ))
#  combinations.append((
#        [["VBFBDTCut"]]+[["IncPreselPtG10"+x] for x in categoriesAll],"CombBDTCat"
#  ))

  #####################################################
  #####################################################
  ## Baseline++2
  #####################################################
  #####################################################

  jet2PtCuts = " && jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40."
  jet01PtCuts = " && !(jetLead_pt > 40. && jetSub_pt > 30. && ptMiss < 40.)"

#  analyses += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
#  analyses += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
  
  #analyses += [["FirstTry","dimuonPt>110. && dijetMass>60. && dijetMass<110"+jet2PtCuts]]
  #analyses += [["FirstTryDimuonPt120","dimuonPt>120. && dijetMass>60. && dijetMass<110"+jet2PtCuts]]
  #analyses += [["FirstTryDimuonPt110","dimuonPt>110. && dijetMass>60. && dijetMass<110"+jet2PtCuts]]

  analyses += [["Jet2CutsVH","dijetMass>60. && dijetMass<110. && dimuonPt>110. && jetLead_eta<2.7 && jetSub_eta<2.7"+jet2PtCuts]]
  #analyses += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)&& !(dimuonPt>110. && dijetMass>60. && dijetMass<110)"+jet2PtCuts]]
  #analyses += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  #analyses += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #analyses += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #analyses += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]

  # Old Categories
  # You must name a new category "Jet2CutsVH"
  #analyses += [["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts]]
  #analyses += [["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #analyses += [["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts]]
  #analyses += [["Jets01PassPtG10"+x,  "dimuonPt>10." +jet01PtCuts] for x in categoriesAll]
  #analyses += [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]

#
#
#  # Jet 0+1 Pass All Cats
#  combinations.append((
#    [["Jets01PassPtG10"+x,"dimuonPt>10."+jet01PtCuts] for x in categoriesAll]
#    ,"Jets01PassCatAll"
#  ))
# 
#  # Jet 0+1 Fail All Cats
#  combinations.append((
#    [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
#    ,"Jets01FailCatAll"
#  ))
# 
#  # Jet 0+1 Pass All Cats
#  combinations.append((
#    [["Jets01PassPtG10"+x,"dimuonPt>10."+jet01PtCuts] for x in categoriesAll]+
#    [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]
#    ,"Jets01SplitCatAll"
#  ))
# 
#  # Jets >=2 Pass + Fail
#  combinations.append((
#    [  
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#     ["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#     ["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#    ],"Jet2SplitCutsGFSplit"
#  ))
# 
#  # Jets 0,1,>=2 Pass + Fail All
#  combinations.append((
#    [["Jets01PassPtG10"+x,"dimuonPt>10."+jet01PtCuts] for x in categoriesAll]+
#    [["Jets01FailPtG10"+x,"!(dimuonPt>10.)"+jet01PtCuts] for x in categoriesAll]+
#    [
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#     ["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#     ["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#    ],"CombSplitAll"
#  ))


##### New combinations
#
#  combinations.append((
#    [["Jets01PassPtG10"+x,"dimuonPt>10."+jet01PtCuts] for x in categoriesAll]+
#    [
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#     ["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#     ["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#    ],"CombJets01PassJet2Split"
#  ))
#
#
#  combinations.append((
#    [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts],
#     ["Jets01PassPtG10BO",  "dimuonPt>10." +jet01PtCuts],
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#     ["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#     ["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#    ],"CombJets01PassBBBOJet2Split"
#  ))
# 
#  combinations.append((
#    [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts],
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#     ["Jet2CutsGFPass","!(deltaEtaJets>3.5 && dijetMass>650.) && (dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#     ["Jet2CutsFailVBFGF","!(deltaEtaJets>3.5 && dijetMass>650.) && !(dijetMass>250. && dimuonPt>50.)"+jet2PtCuts],
#    ],"CombJets01PassBBJet2Split"
#  ))
#
#  combinations.append((
#    [["Jets01PassPtG10BB",  "dimuonPt>10." +jet01PtCuts],
#     ["Jet2CutsVBFPass","deltaEtaJets>3.5 && dijetMass>650."+jet2PtCuts],
#    ],"CombJets01PassBBJet2VBFPass"
#  ))
 
  ## combinations = []

  # Multi-dimensional Optimization of Cuts
  # First two arguments are just like combinations
  # 3rd arg is a dictionary with the keys 
  # that are variables to be cut on with "L" or "G" appended
  # for cutting Less-than, or Greater-then.
  # Finally the last argument specifies whether to "Split"
  # events rather than cutting them.  If set to True,
  # A category passing the cuts will be made and combined
  # with a category failing the events.  If False, the
  # fail events are discarded.
  combinationsCutOpt = []
  combinationsCutOpt.append((
    [["VBFPreselPlus","ptMiss < 40."]],"BDTCutOptPtMissL40",{
        'bdtVBFG':[21,-0.5,0.5],
        },False
  ))
  combinationsCutOpt.append((
    [["VBFPreselPlus","ptMiss < 40."]],"VBFCutBasedOptPtMissL40",{
        'deltaEtaJetsG':[5,3.0,5.0],
        'dijetMassG':[9,300.,700.],
        },False
  ))
  combinationsCutOpt.append((
    [["IncPreselPlus",""]],"NonVBFCutOpt",{
        'dimuonPtG':[16,0.,75.],
        },True
  ))

  ################################################################

  histPostFix="/mDiMu"
  signalNames=["whHmumu125","zhHmumu125"]
  #signalNames=["ggHmumu125","vbfHmumu125","whHmumu125","zhHmumu125"]
  #signalNames=["ggHmumu125","vbfHmumu125"]
  #signalNames=["whHmumu125"]
  #backgroundNames= ["DYJetsToLL","ttbar"]
  backgroundNames= []
  dataDict = {}
  #dataDict["8TeV"] = [
  #  "SingleMuRun2012Av1",
  #  "SingleMuRun2012Av1Recover",
  #  "SingleMuRun2012Bv1",
  #  "SingleMuRun2012Cv1",
  #  "SingleMuRun2012Cv2",
  #  "SingleMuRun2012D",
  #]
  #dataDict["8TeV"] = [
  #  "SingleMuRun2012A_22Jan2013v1",
  #  "SingleMuRun2012B_22Jan2013v1",
  #  "SingleMuRun2012C_22Jan2013v1",
  #  "SingleMuRun2012D_22Jan2013v1",
  #]
  dataDict["8TeV"] = [
    "SingleMuRun2012Av1-22Jan2013",
    "SingleMuRun2012Bv1-22Jan2013",
    "SingleMuRun2012Cv1-22Jan2013",
    "SingleMuRun2012Dv1-22Jan2013",
  ]
  #dataDict["8TeV"] = []
  dataDict["7TeV"] = [
    "SingleMuRun2011Av1",
    "SingleMuRun2011Bv1"
  ]
  #dataDict["7TeV"] = []
  dataDict["14TeV"] = []
  lumiList = [lumiDict["8TeV"],20,25,30]
  lumiList = [lumiDict["8TeV"]]

  controlRegionVeryLow=[60,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]
  if args.higgsMass > 0.0:
    controlRegionLow=[110,args.higgsMass-5.0]
    controlRegionHigh=[args.higgsMass+5.0,160]

  shape=True
  toyData=args.toyData

  if args.combinationsOnly:
    analyses = []
  threads = []
  if not args.cutOpt:
    print("Simple Analyses to run:")
    for a in analyses:
      print("  {0}".format(a))
    print("Combination Analyses to run:")
    for c in combinations:
      print("  {0}".format(c[1]))
      for a in c[0]:
        print("    {0}".format(a))
  
    print("Creating Threads...")
    for p in periods:
      for i in lumiList:
        if p == "7TeV":
          i = lumiDict[p]
        name = i
        if args.higgsMass > 0.0:
          name = args.higgsMass
        for ana in analyses:
          cutsToDo = ""
          if len(ana) == 2:
            cutsToDo = ana[1]
          anaName = ana[0]
          tmp = ThreadedCardMaker(
            #__init__ args:
            directory,[anaName],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p, cutString=cutsToDo,
            #write args:
            outfilename=outDir+anaName+"_"+p+"_"+str(name)+".txt",lumi=i
            )
          threads.append(tmp)
        for comb in combinations:
         threads.append(
          ThreadedCardMaker(
            #__init__ args:
            directory,
            comb[0],
            appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
            controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
            controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
            energyStr=p,
            #write args:
            outfilename=outDir+comb[1]+"_"+p+"_"+str(name)+".txt",lumi=i
          )
         )
        if p == "7TeV":
          break

    # 2011 + 2012 combination
    #if len(periods)>1:
    #    filenameLumi = str(sum([float(lumiDict[p]) for p in periods]))
    #    if args.higgsMass > 0.0:
    #      filenameLumi = str(args.higgsMass)
    #    filenamePeriod = ""
    #    for p in periods:
    #      filenamePeriod += re.sub("TeV","P",p)
    #    filenamePeriod = filenamePeriod.rstrip("P")
    #    filenamePeriod += "TeV"
    #    for ana in analyses:
    #      cutsToDo = ""
    #      if len(ana) == 2:
    #        cutsToDo = ana[1]
    #      anaName = ana[0]
    #      tmp = ThreadedCardMaker(
    #        #__init__ args:
    #        directory,[anaName],
    #        [appendPeriod(signalNames,p) for p in periods],
    #        [appendPeriod(backgroundNames,p) for p in periods],
    #        [dataDict[p] for p in periods],
    #        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
    #        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
    #        energyStr=periods,cutString=cutsToDo,
    #        #write args:
    #        outfilename=outDir+anaName+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
    #        )
    #      threads.append(tmp)
    #    for comb in combinations:
    #     threads.append(
    #      ThreadedCardMaker(
    #        #__init__ args:
    #        directory,
    #        comb[0],
    #        [appendPeriod(signalNames,p) for p in periods],
    #        [appendPeriod(backgroundNames,p) for p in periods],
    #        [dataDict[p] for p in periods],
    #        controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
    #        controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
    #        energyStr=periods,
    #        #write args:
    #        outfilename=outDir+comb[1]+"_"+filenamePeriod+"_"+filenameLumi+".txt",lumi=[float(lumiDict[p]) for p in periods]
    #      )
    #     )
  
  else: #Do cutOpt Stuff
    print "combinationsCutOpt: {0}".format(combinationsCutOpt)
    for p in periods:
      i = lumiDict[p]
      for comb in combinationsCutOpt:
        print "\ncomb: {0}\n".format(comb)
        cutDict = comb[2]
        print "cutDict: {0}".format(cutDict)
        cutBaseString = "1"
        nameBaseString = comb[1]
        cutDictKeyList = cutDict.keys()
        nCutVars = len(cutDictKeyList)
        cutPoints = []
        cutPointsFail = []
        doSplit=comb[3]
        for key,cutIndex in zip(cutDictKeyList,range(nCutVars)):
          cutIndex = str(cutIndex)
          cut = key[:-1]
          cutType = key[-1]
          cutOp = ''
          cutOpFail = ''
          if cutType == "G":
            cutOp = '>'
            cutOpFail = '<='
          elif cutType == "L":
            cutOp = '<'
            cutOpFail = '>='
          else:
            print("Cut Type Not Recognized...exiting.")
            sys.exit(1)
          cutBaseString += " && "+cut+" "+cutOp+" {"+cutIndex+"}"
          nameBaseString += "_"+key+"{"+cutIndex+"}"
          cutData = cutDict[key]
          print "cutData: {0}".format(cutData)
          cutStep = (cutData[2]-cutData[1])/(cutData[0]-1)
          tmpCuts = [cutData[1]+iC*cutStep for iC in range(cutData[0])]
          cutPoints += [tmpCuts]
        print "cutBaseString: "+cutBaseString
        print "nameBaseString: "+nameBaseString
        print "cutPoints: {0}".format(cutPoints)
        assert(len(cutPoints)==nCutVars)
        
        notDone=True
        indexList = [0 for iC in range(nCutVars)]
        lenList = [len(iC) for iC in cutPoints]
        while notDone:
            cutList = [cL[iC] for iC, cL in zip(indexList,cutPoints)]
            cutString = cutBaseString.format(*cutList)
            cutStringFail = "!( "+cutString+" )"
            nameString = nameBaseString.format(*cutList)
            nameString = nameString.replace('.','p')
            print "cutString: "+cutString
            if doSplit:
              print "cutStringFail: "+cutStringFail
            analysisNames = []
            for ana in comb[0]:
              anaName = ana[0]
              anaCutsOrig = '1'
              if len(ana)>1:
                if ana[1] != "":
                  anaCutsOrig = ana[1]
              anaCuts = anaCutsOrig+" && "+cutString
              analysisNames.append([anaName,anaCuts])
              if doSplit:
                anaCuts = anaCutsOrig+" && "+cutStringFail
                analysisNames.append([anaName+"Fail",anaCuts])
            
            print("analysisNames: {0}".format(analysisNames))
            threads.append(
             ThreadedCardMaker(
               #__init__ args:
               directory,
               analysisNames,
               appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
               controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,
               controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,sigInjectMass=args.signalInjectMass,
               energyStr=p,
               #write args:
               outfilename=outDir+nameString+"_"+p+"_"+str(i)+".txt",lumi=i,
               cutString=""
             )
            )
            for ind,iLen,iInd in reversed(zip(indexList,lenList,range(nCutVars))):
              if ind < iLen-1:
                indexList[iInd] += 1
                break
              else:
                indexList[iInd] = 0
            if sum(indexList) == 0:
              notDone = False

  nThreads = len(threads)
  print("nProcs: {0}".format(NPROCS))
  print("nCards: {0}".format(nThreads))

  threadsNotStarted = copy.copy(threads)
  threadsRunning = []
  threadsDone = []
  while True:
    iThread = 0
    while iThread < len(threadsRunning):
        alive = threadsRunning[iThread].is_alive()
        if not alive:
          tmp = threadsRunning.pop(iThread)
          threadsDone.append(tmp)
        else:
          iThread += 1

    nRunning = len(threadsRunning)
    if nRunning < NPROCS and len(threadsNotStarted) > 0:
        tmp = threadsNotStarted.pop()
        tmp.start()
        threadsRunning.append(tmp)

    nRunning = len(threadsRunning)
    nNotStarted = len(threadsNotStarted)
    if nRunning == 0 and nNotStarted == 0:
          break

    time.sleep(0.1)

  if RUNSIMPLELIMITS:
    shutil.copy("etc/notlxbatch_simple.sh",outDir+"notlxbatch.sh")
    shutil.copy("etc/lxbatch_simple.sh",outDir+"lxbatch.sh")
  else:
    shutil.copy("etc/notlxbatch.sh",outDir+"notlxbatch.sh")
    shutil.copy("etc/lxbatch.sh",outDir+"lxbatch.sh")

  shutil.copy("etc/run.sh",outDir+"run.sh")
  shutil.copy("etc/uftrig01.sh",outDir+"uftrig01.sh")
      
  shutil.copy("etc/nuisanceDiff.py",outDir+"diffNuisances.py")
  shutil.copy("etc/fitNormsToText_mlfit.py",outDir+"fitNormsToText_mlfit.py")
  shutil.copy("etc/myNuisancePrinter.py",outDir+"myNuisancePrinter.py")
  shutil.copy("etc/hpcTemplate.sh",outDir+"hpcTemplate.sh")
  shutil.copy("etc/runHPC.sh",outDir+"runHPC.sh")
  shutil.copy("etc/getStatus2.sh",outDir+"getStatus2.sh")
  shutil.copy("etc/gof.sh",outDir+"gof.sh")
  shutil.copy("etc/gofHPC_Template.sh",outDir+"gofHPC_Template.sh")
  shutil.copy("etc/runHPC_GOF.sh",outDir+"runHPC_GOF.sh")
  shutil.copy("etc/runHPC_Compat.sh",outDir+"runHPC_Compat.sh")
  shutil.copy("etc/compatHPC_Template.sh",outDir+"compatHPC_Template.sh")
  shutil.copy("etc/lxbatch_LEE.sh",outDir+"lxbatch_LEE.sh")
  shutil.copy("etc/runLEE.sh",outDir+"runLEE.sh")
  shutil.copy("etc/combAllText.py",outDir+"combAllText.py")
  shutil.copy("etc/checkWS.py",outDir+"checkWS.py")
  shutil.copy("etc/findExpected.sh",outDir+"findExpected.sh")

  for iexef in glob.glob(outDir+"*.sh")+glob.glob(outDir+"*.py"):
    os.chmod(iexef, os.stat(iexef).st_mode | stat.S_IXUSR | stat.S_IXGRP | stat.S_IXOTH)

  print "makeCards.py done."
