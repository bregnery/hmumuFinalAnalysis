#!/usr/bin/env python

import argparse
parser = argparse.ArgumentParser(description="Makes cards for use in the CMS Combine tool.")
parser.add_argument("--signalInject", help="Inject Signal with Strength into data_obs",type=float,default=0.0)
parser.add_argument("--toyData", help="Make Toy Data from PDFs for data_obs",action="store_true",default=False)
args = parser.parse_args()

import math
import ROOT as root
from helpers import *
import datetime
import sys
import os.path
import copy
import multiprocessing
import time
myThread = multiprocessing.Process

from ROOT import gSystem
gSystem.Load('libRooFit')

root.gErrorIgnoreLevel = root.kWarning
#root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.WARNING)
root.RooMsgService.instance().setGlobalKillBelow(root.RooFit.ERROR)
PRINTLEVEL = root.RooFit.PrintLevel(-1) #For MINUIT

NPROCS = 2

BAKUNC = 0.1

BAKUNCON = True
SIGUNCON = False

SIGGAUS = True

from xsec import *

if scaleHiggsBy != 1.0:
  print("Error: higgs xsec is scaled!!! Return to 1. Exiting.")
  sys.exit(1)
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

def getRooVars(directory,signalNames,histNameBase,analysis):
    hist = None
    is2D = False
    for name in signalNames:
      filename = directory+name+".root"
      histName = histNameBase+analysis
      #print("file name: {0}".format(filename))
      #print("hist name: {0}".format(histName))
      tmpF = root.TFile(filename)
      hist = tmpF.Get(histName)
      break
    if hist.InheritsFrom("TH2"):
      is2D = True

    x = root.RooRealVar('mMuMu','mMuMu',
                    hist.GetXaxis().GetXmin(),
                    hist.GetXaxis().GetXmax()
                    )
    if is2D:
      y = root.RooRealVar('mva','mva',
                    hist.GetYaxis().GetXmin(),
                    hist.GetYaxis().GetXmax()
                    )
    if is2D:
      return [x,y]
    else:
      return [x]

###################################################################################

class MassPDFBakMSSM:
  def __init__(self,name,hist,massLowRange,massHighRange,massVeryLowRange,rooVars=None,smooth=False,hack=True):
    if rooVars == None:
        print("Error: MassPDFBakMSSM requires rooVars list of variables, exiting.")
        sys.exit(1)
    if len(massVeryLowRange)==0:
        print("Error: MassPDFBakMSSM requires verylow mass range, exiting.")
        sys.exit(1)

    self.debug = ""
    self.debug += "### MassPDFBakMSSM: "+name+"\n"

    maxMass = massHighRange[1]
    minMass = massVeryLowRange[0]
    mMuMu = root.RooRealVar("mMuMu","mMuMu",minMass,maxMass)
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])

    bwmZ = root.RooRealVar("bwmZ","bwmZ",85,95)
    bwSig = root.RooRealVar("bwSig","bwSig",0.0,30.0)
    bwLambda = root.RooRealVar("bwLambda","bwLambda",-1e-03,-1e-01,-1e-04)
    bwMmumu = root.RooGenericPdf("bwMmumu","exp(bwLambda*mMuMu)*bwSig/(((mMuMu-bwmZ)*(mMuMu-bwmZ) + bwSig*bwSig*0.25))",root.RooArgList( mMuMu, bwLambda, bwSig, bwmZ))

    phoMmumu = root.RooGenericPdf("phoMmumu","exp(bwLambda*mMuMu)/pow(mMuMu,2)",root.RooArgList(mMuMu,bwLambda))

    mixParam = root.RooRealVar("mixParam","mixParam",0.5,0,1)

    pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(bwMmumu,phoMmumu),root.RooArgList(mixParam))
    
    tmpAxis = hist.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)
    lowBin, highBin = getXbinsHighLow(hist,minMass,maxMass)
    nBinsX = highBin - lowBin + 1
    normalization = getIntegralAll(hist,[massLowRange[0],massHighRange[1]])
    self.debug += "# bak Hist no bounds: {0:.3g}\n".format(hist.Integral())
    self.debug += "# bak Hist bounds:    {0:.3g}\n".format(normalization)

    mMuMuRooDataHist = root.RooDataHist(name+"DataHist",name+"DataHist",root.RooArgList(mMuMu),hist)

    bwMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    bwmZ.setConstant(True)
    bwSig.setConstant(True)

    #phoMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    #bwLambda.setConstant(True)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    plotMmumu = mMuMu.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("low,high"))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("phoMmumu"),root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass),root.RooFit.LineColor(root.kGreen+1))

    mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass)
    nominalHist = pdfMmumu.createHistogram("pdf2dHist",mMuMu,mMuMuBinning)
    nominalHist.Scale(normalization/getIntegralAll(nominalHist,[massLowRange[0],massHighRange[1]]))

    self.name = name
    self.hist = hist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.lowBin = lowBin
    self.highBin = highBin
    self.nBinsX = nBinsX
    self.bwmZ = bwmZ
    self.bwSig = bwSig
    self.bwLambda = bwLambda
    self.bwMmumu = bwMmumu
    self.phoMmumu = phoMmumu
    self.mixParam = mixParam
    self.pdfMmumu = pdfMmumu
    self.mMuMuBinning = mMuMuBinning
    self.nominalHist = nominalHist
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.plotMmumu = plotMmumu
    self.chi2 = chi2

    self.debug += "# nominal Integral: {0:.3g}\n".format(getIntegralAll(nominalHist))
    self.debug += "# BW Sigma: {0:.3g} +/- {1:.3g}\n".format(bwSig.getVal(),bwSig.getError())
    self.debug += "# BW mZ:    {0:.3g} +/- {1:.3g}\n".format(bwmZ.getVal(),bwmZ.getError())
    self.debug += "# Lam:      {0:.3g} +/- {1:.3g}\n".format(bwLambda.getVal(),bwLambda.getError())
    self.debug += "# BW Coef:  {0:.3g} +/- {1:.3g}\n".format(mixParam.getVal(),mixParam.getError())
    self.debug += "# chi2/ndf: {0:.3g}\n".format(chi2.getVal()/(nBinsX-1))

    ## Error time

    self.errNames = []
    self.errHists = {}
    if BAKUNCON:
      for errVar in [bwmZ,bwSig,bwLambda,mixParam]:
        val = errVar.getVal()
        err = errVar.getError()
        varName = errVar.GetName()

        errVar.setVal(val+err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        upHist = pdfMmumu.createHistogram(varName+"Up",mMuMu,mMuMuBinning)

        errVar.setVal(val-err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        downHist = pdfMmumu.createHistogram(varName+"Down",mMuMu,mMuMuBinning)

        errVar.setVal(val)

        upHist.Scale(normalization/getIntegralAll(upHist,[massLowRange[0],massHighRange[1]]))
        downHist.Scale(normalization/getIntegralAll(downHist,[massLowRange[0],massHighRange[1]]))

        setattr(self,varName+"UpHist",upHist)
        setattr(self,varName+"DownHist",downHist)

        self.errNames.append(varName)
        self.errHists[varName+"Up"] = upHist
        self.errHists[varName+"Down"] = downHist

  def writeDebugHistsToCurrTDir(self,compareHist=None):
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    self.plotMmumu.Draw()
    canvas.SetName(self.name+"Canvas")
    canvas.Write()

class MassPDFBak:
  def __init__(self,name,hist,massLowRange,massHighRange,massVeryLowRange,rooVars=None,smooth=False,hack=True):
    if rooVars == None:
        print("Error: MassPDFBak requires rooVars list of variables, exiting.")
        sys.exit(1)
    if len(massVeryLowRange)==0:
        print("Error: MassPDFBak requires verylow mass range, exiting.")
        sys.exit(1)

    self.debug = ""
    self.debug += "### MassPDFBak: "+name+"\n"

    maxMass = massHighRange[1]
    minMass = massVeryLowRange[0]
    mMuMu = root.RooRealVar("mMuMu2","mMuMu",minMass,maxMass)
    mMuMu.setRange("z",88,94)
    mMuMu.setRange("verylow",massVeryLowRange[0],massVeryLowRange[1])
    mMuMu.setRange("low",massLowRange[0],massLowRange[1])
    mMuMu.setRange("high",massHighRange[0],massHighRange[1])
    mMuMu.setRange("signal",massLowRange[1],massHighRange[0])

    voitWidth = root.RooRealVar("voitWidth","voitWidth",2.4952)
    voitmZ = root.RooRealVar("voitmZ","voitmZ",85,95)
    voitSig = root.RooRealVar("voitSig","voitSig",0.0,30.0)
    voitMmumu = root.RooVoigtian("voitMmumu","voitMmumu",mMuMu,voitmZ,voitWidth,voitSig)

    expParam = root.RooRealVar("expParam","expParam",-1,0)
    expMmumu = root.RooExponential("expMmumu","expMmumu",mMuMu,expParam)

    mixParam = root.RooRealVar("mixParam","mixParam",0,1)

    pdfMmumu = root.RooAddPdf("pdfMmumu","pdfMmumu",root.RooArgList(voitMmumu,expMmumu),root.RooArgList(mixParam))
    
    tmpAxis = hist.GetXaxis()
    lowBin = tmpAxis.FindBin(minMass)
    highBin = tmpAxis.FindBin(maxMass)
    lowBin, highBin = getXbinsHighLow(hist,minMass,maxMass)
    nBinsX = highBin - lowBin + 1
    normalization = getIntegralAll(hist,[massLowRange[0],massHighRange[1]])
    self.debug += "# bak Hist no bounds: {0:.3g}\n".format(hist.Integral())
    self.debug += "# bak Hist bounds:    {0:.3g}\n".format(normalization)

    mMuMuRooDataHist = root.RooDataHist(name+"DataHist",name+"DataHist",root.RooArgList(mMuMu),hist)

    voitMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("z"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    voitmZ.setConstant(True)
    voitSig.setConstant(True)

    expMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    expParam.setConstant(True)
    
    pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.Range("low,high"),root.RooFit.SumW2Error(False),PRINTLEVEL)
    chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

    plotMmumu = mMuMu.frame()

    mMuMuRooDataHist.plotOn(plotMmumu)
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Range("low,high"))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass))
    pdfMmumu.plotOn(plotMmumu,root.RooFit.Components("expMmumu"),root.RooFit.LineStyle(2),root.RooFit.Range(minMass,maxMass),root.RooFit.LineColor(root.kGreen+1))

    mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass)
    nominalHist = pdfMmumu.createHistogram("pdf2dHist",mMuMu,mMuMuBinning)
    nominalHist.Scale(normalization/getIntegralAll(nominalHist,[massLowRange[0],massHighRange[1]]))

    self.name = name
    self.hist = hist
    self.mMuMuRooDataHist = mMuMuRooDataHist
    self.lowBin = lowBin
    self.highBin = highBin
    self.nBinsX = nBinsX
    self.voitWidth = voitWidth
    self.voitmZ = voitmZ
    self.voitSig = voitSig
    self.voitMmumu = voitMmumu
    self.expParam = expParam
    self.expMmumu = expMmumu
    self.mixParam = mixParam
    self.pdfMmumu = pdfMmumu
    self.mMuMuBinning = mMuMuBinning
    self.nominalHist = nominalHist
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu
    self.plotMmumu = plotMmumu
    self.chi2 = chi2

    self.debug += "# nominal Integral: {0:.3g}\n".format(getIntegralAll(nominalHist))
    self.debug += "# V Width: {0:.3g} +/- {1:.3g}\n".format(voitWidth.getVal(),voitWidth.getError())
    self.debug += "# V Sigma: {0:.3g} +/- {1:.3g}\n".format(voitSig.getVal(),voitSig.getError())
    self.debug += "# V mZ:    {0:.3g} +/- {1:.3g}\n".format(voitmZ.getVal(),voitmZ.getError())
    self.debug += "# Exp Par: {0:.3g} +/- {1:.3g}\n".format(expParam.getVal(),expParam.getError())
    self.debug += "# V Coef:  {0:.3g} +/- {1:.3g}\n".format(mixParam.getVal(),mixParam.getError())
    self.debug += "# chi2/ndf: {0:.3g}\n".format(chi2.getVal()/(nBinsX-1))

    ## Error time

    self.errNames = []
    self.errHists = {}
    if BAKUNCON:
      for errVar in [voitmZ,voitSig,expParam,mixParam]:
        val = errVar.getVal()
        err = errVar.getError()
        varName = errVar.GetName()

        errVar.setVal(val+err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        upHist = pdfMmumu.createHistogram(varName+"Up",mMuMu,mMuMuBinning)

        errVar.setVal(val-err)
        pdfMmumu.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen-1),root.RooFit.Range(110,150),root.RooFit.LineStyle(3))
        downHist = pdfMmumu.createHistogram(varName+"Down",mMuMu,mMuMuBinning)

        errVar.setVal(val)

        upHist.Scale(normalization/getIntegralAll(upHist,[massLowRange[0],massHighRange[1]]))
        downHist.Scale(normalization/getIntegralAll(downHist,[massLowRange[0],massHighRange[1]]))

        setattr(self,varName+"UpHist",upHist)
        setattr(self,varName+"DownHist",downHist)

        self.errNames.append(varName)
        self.errHists[varName+"Up"] = upHist
        self.errHists[varName+"Down"] = downHist

  def writeDebugHistsToCurrTDir(self,compareHist=None):
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    self.plotMmumu.Draw()
    canvas.SetName(self.name+"Canvas")
    canvas.Write()

class MassPDFSig:
  def __init__(self,name,names,histList,histMapErr,massRange,rooVars=None):
    if rooVars == None:
        print("Error: MassPDFSig requires rooVars list of variables, exiting.")
        sys.exit(1)

    self.debug = ""
    self.debug += "### MassPDFSig: "+name+"\n"

    maxMass = massRange[1]
    minMass = massRange[0]
    mMuMu = root.RooRealVar("mMuMu3","mMuMu",minMass,maxMass)

    self.names = names
    self.name = name
    self.maxMass = maxMass
    self.minMass = minMass
    self.mMuMu = mMuMu

    self.meanList = []
    self.widthList = []
    self.pdfList = []
    self.histList = []
    self.mMuMuRooDataHistList = []
    self.nominalHistList = []
    self.plotMmumuList = []
    self.chi2List = []

    self.meanMapErr = {}
    self.widthMapErr = {}
    self.pdfMapErr = {}
    self.histMapErr = {}

    for i,n in zip(range(len(names)),names):
      hist = histList[i]
      self.debug += "# "+n+":\n"
      width = root.RooRealVar("width_"+n,"width_"+n,0.1,15.)
      mean = root.RooRealVar("mean_"+n,"mean_"+n,100.,150.)

      pdfMmumu = root.RooGaussian("pdfMmumu_"+n,"pdfMmumu_"+n,mMuMu,mean,width)

      tmpAxis = hist.GetXaxis()
      lowBin = tmpAxis.FindBin(minMass)
      highBin = tmpAxis.FindBin(maxMass)
      lowBin, highBin = getXbinsHighLow(hist,minMass,maxMass)
      nBinsX = highBin - lowBin + 1
      normalization = getIntegralAll(hist,[minMass,maxMass])
      self.debug += "#   Int no bounds: {0:.3g}\n".format(hist.Integral())
      self.debug += "#   Int bounds:    {0:.3g}\n".format(normalization)
  
      mMuMuRooDataHist = root.RooDataHist(name+"DataHist_"+n,name+"DataHist_"+n,root.RooArgList(mMuMu),hist)
  
      pdfMmumu.fitTo(mMuMuRooDataHist,root.RooFit.SumW2Error(False),PRINTLEVEL)
      chi2 = pdfMmumu.createChi2(mMuMuRooDataHist)

      plotMmumu = mMuMu.frame()

      mMuMuRooDataHist.plotOn(plotMmumu)
      pdfMmumu.plotOn(plotMmumu)

      mMuMuBinning = root.RooFit.Binning(nBinsX,minMass,maxMass)
      nominalHist = pdfMmumu.createHistogram(name+"_"+n,mMuMu,mMuMuBinning)
      nominalHist.Scale(normalization/getIntegralAll(nominalHist,[minMass,maxMass]))

      self.meanList.append(mean)
      self.widthList.append(mean)
      self.pdfList.append(pdfMmumu)
      self.histList.append(hist)
      self.mMuMuRooDataHistList.append(mMuMuRooDataHist)
      self.nominalHistList.append(nominalHist)
      self.plotMmumuList.append(plotMmumu)
      self.chi2List.append(chi2)

      self.debug += "#   nominal Integral: {0:.3g}\n".format(getIntegralAll(nominalHist))
      self.debug += "#   width: {0:.3g} +/- {1:.3g}\n".format(width.getVal(),width.getError())
      self.debug += "#   mean: {0:.3g} +/- {1:.3g}\n".format(mean.getVal(),mean.getError())
      self.debug += "#   chi2/ndf: {0:.3g}\n".format(chi2.getVal()/(nBinsX-1))

      ## Error time
      for key in histMapErr:
        errHist = histMapErr[key][i]

        widthE = root.RooRealVar("width_"+n+"_"+key,"width_"+n+"_"+key,0.1,15.)
        meanE = root.RooRealVar("mean_"+n+"_"+key,"mean_"+n+"_"+key,100.,150.)
        pdfMmumuE = root.RooGaussian("pdfMmumu_"+n+"_"+key,"pdfMmumu_"+n+"_"+key,mMuMu,meanE,widthE)

        mMuMuRooDataHistE = root.RooDataHist(name+"DataHist_"+n+"_"+key,name+"DataHist_"+n+"_"+key,root.RooArgList(mMuMu),errHist)
        pdfMmumuE.fitTo(mMuMuRooDataHistE,root.RooFit.SumW2Error(False),PRINTLEVEL)
        chi2E = pdfMmumuE.createChi2(mMuMuRooDataHistE)

        histE = pdfMmumuE.createHistogram(name+"_"+n+"_"+key,mMuMu,mMuMuBinning)
        histE.Scale(normalization/getIntegralAll(histE,[minMass,maxMass]))

        pdfMmumuE.plotOn(plotMmumu,root.RooFit.LineColor(root.kGreen+1),root.RooFit.LineStyle(2))

        if self.meanMapErr.has_key(key):
            self.meanMapErr[key].append(meanE)
            self.widthMapErr[key].append(widthE)
            self.pdfMapErr[key].append(pdfMmumuE)
            self.histMapErr[key].append(histE)
        else:
            self.meanMapErr[key] = [meanE]
            self.widthMapErr[key] = [widthE]
            self.pdfMapErr[key] = [pdfMmumuE]
            self.histMapErr[key] = [histE]

        self.debug += "#   Error: {0}\n".format(key)
        self.debug += "#     width: {0:.3g} +/- {1:.3g}\n".format(widthE.getVal(),widthE.getError())
        self.debug += "#     mean: {0:.3g} +/- {1:.3g}\n".format(meanE.getVal(),meanE.getError())
        self.debug += "#     chi2/ndf: {0:.3g}\n".format(chi2E.getVal()/(nBinsX-1))

  def writeDebugHistsToCurrTDir(self,compareHist=None):
    canvas = root.TCanvas("canvas")
    canvas.cd()
    #canvas.SetLogy(1)
    for p,n in zip(self.plotMmumuList,self.names):
      p.Draw()
      canvas.SetName(n+"Canvas")
      canvas.Write()

###################################################################################

class Analysis:
  def __init__(self,directory,signalNames,backgroundNames,dataNames,analysis,x,y,controlRegionLow,controlRegionHigh,histNameBase="mDiMu",bakShape=False,rebin=[],histNameSuffix="",controlRegionVeryLow=[]):
    self.bakShape = bakShape
    self.sigNames = signalNames
    self.bakNames = backgroundNames
    self.datNames = dataNames
    self.controlRegionVeryLow = controlRegionVeryLow
    self.controlRegionLow = controlRegionLow
    self.controlRegionHigh = controlRegionHigh
    self.analysis = analysis

    self.is2D = False
    if y != None:
      self.is2D = True
    self.x = x
    self.y = y
    self.x1d = None

    self.sigFiles = []
    self.sigHistsRaw = []
    for name in signalNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.sigFiles.append(tmpF)
      self.sigHistsRaw.append(tmpH)
      if tmpH.InheritsFrom("TH2"):
        self.is2D = True

    # Signal Shape systematics
    self.sigErrHistsMap = {}
    if SIGUNCON:
      for f in self.sigFiles:
        name = histNameBase+analysis+histNameSuffix
        name = name.split("/")
        tmpDirKey = f.GetKey(name[0]) #Will break in main dir
        tmpDir = tmpDirKey.ReadObj()
        tmpDir.Print()
        for key in tmpDir.GetListOfKeys():
          matchUp = re.match(name[1]+".+Up",key.GetName())
          matchDown = re.match(name[1]+".+Down",key.GetName())
          if matchUp or matchDown:
            self.sigErrHistsMap[re.sub(name[1],"",key.GetName())] = []
        break
      for f in self.sigFiles:
        for sysName in self.sigErrHistsMap:
          tmpHist = f.Get(histNameBase+analysis+histNameSuffix+sysName)
          self.sigErrHistsMap[sysName].append(tmpHist)

    self.bakFiles = []
    self.bakHistsRaw = []
    for name in backgroundNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.bakFiles.append(tmpF)
      self.bakHistsRaw.append(tmpH)

    self.datFiles = []
    self.datHists = []
    for name in dataNames:
      tmpF = root.TFile(directory+name+".root")
      tmpH = tmpF.Get(histNameBase+analysis+histNameSuffix)
      self.datFiles.append(tmpF)
      self.datHists.append(tmpH)

    #Rebin
    rb = rebin
    if type(rb) != list:
      print("Error: Analysis.rebin: argument must be a list!!  Exiting.")
      sys.exit(1)
    if len(rb) == 2 and self.is2D:
        for hist in self.sigHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.bakHistsRaw:
          hist.Rebin2D(*rb)
        for hist in self.datHists:
          hist.Rebin2D(*rb)
        for err in self.sigErrHistsMap:
          for hist in self.sigErrHistsMap[err]:
            hist.Rebin2D(*rb)
    elif len(rb) == 1 and not self.is2D:
        for hist in self.sigHistsRaw:
          hist.Rebin(*rb)
        for hist in self.bakHistsRaw:
          hist.Rebin(*rb)
        for hist in self.datHists:
          hist.Rebin(*rb)
        for err in self.sigErrHistsMap:
          for hist in self.sigErrHistsMap[err]:
            hist.Rebin(*rb)
    elif len(rb) == 0:
      pass
    else:
      print("Error: Analysis.rebin: argument must be len 0, 1, or 2 list!!  Exiting.")
      print("  Must also be same length as dimension of hist, if not 0.")
      sys.exit(1)

    effMap = {}
    xsecMap = {}
    lowBin = 0
    highBin = self.sigHistsRaw[0].GetNbinsX()+1
    massBounds = [controlRegionLow[0],controlRegionHigh[1]]
    self.massBounds = massBounds

    for hist in self.sigHistsRaw:
      vetoOutOfBoundsEvents(hist,boundaries=massBounds)

    self.xsecBakTotal = 0.0
    self.xsecBakList = []
    self.effBakList = []
    self.bakHists = []
    self.bakHistTotal = None
    for h,name in zip(self.bakHistsRaw,backgroundNames):
      counts = getIntegralAll(h,boundaries=massBounds)
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

    self.bakHistTotalReal = self.bakHistTotal.Clone("data_obs")

    self.bakShape = bakShape
    self.bakShapeMkr = None

    self.dataCountsTotal = None
    self.datHistTotal = None
    for h,name in zip(self.datHists,dataNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      if self.dataCountsTotal == None:
        self.dataCountsTotal = counts
      else:
        self.dataCountsTotal += counts
      if self.datHistTotal == None:
        self.datHistTotal = h.Clone("data_obs")
      else:
        self.datHistTotal.Add(h)

    if SIGGAUS:
      self.sigShapeMkr = MassPDFSig(analysis,signalNames,self.sigHistsRaw,
                                self.sigErrHistsMap,
                                [self.controlRegionLow[0],self.controlRegionHigh[1]],
                                rooVars = [self.x,self.y]
                                )
      for i,name in zip(range(len(signalNames)),signalNames):
        tmp = self.sigShapeMkr.nominalHistList[i]
        tmp.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])
        for err in self.sigErrHistsMap:
          tmp = self.sigShapeMkr.histMapErr[err][i]
          tmp.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])

    self.xsecSigTotal = 0.0
    self.xsecSigList = []
    self.effSigList = []
    self.sigHists = []
    for h,name in zip(self.sigHistsRaw,signalNames):
      counts = getIntegralAll(h,boundaries=massBounds)
      eff = counts/nEventsMap[name]*efficiencyMap[getPeriod(name)]
      xs = eff*xsec[name]
      self.xsecSigTotal += xs
      self.xsecSigList.append(xs)
      self.effSigList.append(eff)
      h.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])
      self.sigHists.append(h)

    for err in self.sigErrHistsMap:
      for h,name in zip(self.sigErrHistsMap[err],signalNames):
        h.Scale(xsec[name]/nEventsMap[name]*efficiencyMap[getPeriod(name)])
    self.sigErrNames = set()
    for errLong in self.sigErrHistsMap:
        errLong = re.sub(r"Up$","",errLong)
        errLong = re.sub(r"Down$","",errLong)
        if not errLong in self.sigErrNames:
            self.sigErrNames.add(errLong)

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
  def getSigHist(self,sigName):
    result = -1.0
    if self.sigNames.count(sigName)>0:
        i = self.sigNames.index(sigName)
        if SIGGAUS:
          result = self.sigShapeMkr.nominalHistList[i]
        else:
          result = self.sigHists[i]
    return result
  def getSigSystHist(self,sigName,systName):
    result = -1.0
    l = self.sigErrHistsMap[systName]
    if SIGGAUS:
      l = self.sigShapeMkr.histMapErr[systName]
    if self.sigNames.count(sigName)>0:
        i = self.sigNames.index(sigName)
        result = l[i]
    return result
  def getBakHist(self,bakName):
    result = -1.0
    if self.bakNames.count(bakName)>0:
        i = self.bakNames.index(bakName)
        result = self.bakHists[i]
    return result
  def getBakHistTotal(self,lumi):
    analysis = self.analysis
    bakShape = self.bakShape
    if self.datHistTotal == None:
      self.bakHistTotal = self.bakHistTotalReal.Clone("stuff")
      self.bakHistTotal.Scale(lumi)
    else:
      self.bakHistTotal = self.datHistTotal.Clone("stuff")
    if bakShape and self.is2D:
      bakShapeMkr = MVAvMassPDFBak("pdfHists_"+analysis,
                                self.bakHistTotal,
                                self.controlRegionLow,self.controlRegionHigh,
                                rooVars = [self.x,self.y]
                                )
      self.bakShapeMkr = bakShapeMkr
      self.bakHistTotal = bakShapeMkr.hackHist
      #self.xsecBakTotal = getIntegralAll(self.bakHistTotal,boundaries=massBounds)
    elif bakShape:
      bakShapeMkr = MassPDFBak("pdfHists_"+analysis,
                                self.bakHistTotal,
                                self.controlRegionLow,self.controlRegionHigh,
                                self.controlRegionVeryLow,
                                rooVars = [self.x]
                                )
      self.bakShapeMkr = bakShapeMkr
      self.bakHistTotal = bakShapeMkr.nominalHist
      #self.xsecBakTotal = getIntegralAll(self.bakHistTotal,boundaries=self.massBounds)
    return self.bakHistTotal
  def dump(self):
    print("##########################################")
    print("SigXSecTotal: {0}".format(self.xsecSigTotal))
    print("BakXSecTotal: {0}".format(self.xsecBakTotal))
    print("signames: {0}".format(self.sigNames))
    print("baknames: {0}".format(self.bakNames))
    print("sighists: {0}".format(self.sigHists))
    for i in self.sigHists:
        i.Print()
    print("bakhists: {0}".format(self.bakHists))
    for i in self.bakHists:
        i.Print()
    print("bakHistTotal: {0}".format(self.bakHistTotal))
    self.bakHistTotal.Print()
    print("##########################################")

###################################################################################

class DataCardMaker:
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap=None,histNameBase="",controlRegionLow=[80,115],controlRegionHigh=[135,150],controlRegionVeryLow=[],bakShape=False,rebin=[],histNameSuffix=""):
    channels = []
    self.channelNames = copy.deepcopy(analysisNames)
    self.is2D = False

    x = None
    y = None
    for analysis in analysisNames:
      tmpList = getRooVars(directory,signalNames,histNameBase,analysis+histNameSuffix)
      x = tmpList[0]
      if len(tmpList)==2:
        self.is2D = True
        y = tmpList[1]
    self.x = x
    self.y = y
    self.shape = bakShape

    for analysis in analysisNames:
      tmp = Analysis(directory,signalNames,backgroundNames,dataNames,analysis,x,y,controlRegionLow,controlRegionHigh,controlRegionVeryLow=controlRegionVeryLow,histNameBase=histNameBase,bakShape=bakShape,rebin=rebin,histNameSuffix=histNameSuffix)
      channels.append(tmp)
    self.channels = channels

    if nuisanceMap == None:
      self.nuisance = {}
    else:
      self.nuisance = nuisanceMap

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
    self.sigNames = signalNames
    self.bakNames = backgroundNames

    if self.channelNames.count("")>0:
      i = self.channelNames.index("")
      self.channelNames[i] = "Inc"

  def write(self,outfilename,lumi,sigInject=0.0):
    print("Writing Card: {0}".format(outfilename))
    lumi *= 1000.0
    nuisance = self.nuisance
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format(len(nuisance)))
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
      if channel.dataCountsTotal == None:
        print("Writing Pretend Data Counts")
        counts = channel.getBakXSecTotal()*lumi
        if sigInject != 0.0:
          counts += channel.getSigXSecTotal()*lumi*sigInject
        observationFormatList.append(int(counts))
      else:
        print("Writing Real Data Counts")
        observationFormatList.append(channel.dataCountsTotal)
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
        iProc = -len(channel.sigNames)+1
        for sigName in self.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigXSec(sigName)*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1
        for bakName in self.bakNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(bakName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getBakXSec(bakName)*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
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

    for nu in nuisance:
      thisNu = nuisance[nu]
      formatString = "{0:<8} {1:^4} "
      formatList = [nu,"lnN"]
      iParam = 2
      for channel,channelName in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu.has_key(sigName):
              value = thisNu[sigName]+1.0
            formatList.append(value)
            iParam += 1
          for bakName in self.bakNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu.has_key(bakName):
              value = thisNu[bakName]+1.0
            formatList.append(value)
            iParam += 1
      formatString += "\n"
      #print formatString
      #print formatList
      outfile.write(formatString.format(*formatList))

    #Debugging
    outfile.write("#################################\n")
    for channel,channelName in zip(self.channels,self.channelNames):
        outfile.write("#\n")
        outfile.write("#info: channel {0}: \n".format(channelName))
        for sigName in self.sigNames:
          outfile.write("#  {0} XS, Eff: {1:.3g}, {2:.3%} \n".format(sigName,
                        channel.getSigXSec(sigName), channel.getSigEff(sigName)))
        outfile.write("#   Mass: {0} \n".format(channel.massBounds))
    outfile.close()

###################################################################################

class ShapeDataCardMaker(DataCardMaker):
  def __init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap=None,histNameBase="",rebin=[],useTH1=False,controlRegionLow=[80,115],controlRegionHigh=[135,200],controlRegionVeryLow=[],bakShape=False,histNameSuffix="",toyData=False):
    DataCardMaker.__init__(self,directory,analysisNames,signalNames,backgroundNames,dataNames,nuisanceMap,histNameBase,controlRegionLow,controlRegionHigh,controlRegionVeryLow=controlRegionVeryLow,bakShape=bakShape,rebin=rebin,histNameSuffix=histNameSuffix)

    self.useTH1 = useTH1
    self.controlRegionHigh = controlRegionHigh
    self.controlRegionLow = controlRegionLow
    self.controlRegionVeryLow = controlRegionVeryLow
    self.toyData = toyData

  def makeRFHistWrite(self,channel,hist,thisDir,isData=True,compareHist=None,writeBakShape=False):
    thisDir.cd()
    is2D = hist.InheritsFrom("TH2")
    if self.useTH1 and not is2D:
        hist.Write()
        return
    if not is2D:
      hist = shrinkTH1(hist,self.controlRegionLow[0],self.controlRegionHigh[1])
    x = channel.x
    origHist = hist
    if is2D:
      hist = hist2to1(hist)
      hist.SetName(re.sub("_1d","",hist.GetName()))
      if channel.x1d == None:
        channel.x1d = root.RooRealVar("x1d","x1d",0,hist.GetNbinsX()+2)
      x = channel.x1d
    rfHist = root.RooDataHist(hist.GetName(),hist.GetName(),root.RooArgList(root.RooArgSet(x)),hist)
    rfHistPdf = root.RooHistPdf(hist.GetName(),hist.GetName(),root.RooArgSet(x),rfHist)
    if isData:
      rfHist.Write()
    else:
      rfHistPdf.Write()
    debugDir = thisDir.FindObject('debug')
    if  debugDir==None:
      debugDir = thisDir.mkdir("debug")
    debugDir.cd()
    hist.Write()
    if is2D:
      xBinning = root.RooFit.Binning(origHist.GetNbinsX())
      yBinning = root.RooFit.Binning(origHist.GetNbinsY())
      x = channel.x
      y = channel.y
      rfHistTH2 = rfHist.createHistogram(origHist.GetName()+"rfHist2d",x,xBinning,root.RooFit.YVar(y,yBinning))
      rfHistPdfTH2 = rfHistPdf.createHistogram(origHist.GetName()+"rfHistPdf2d",x,xBinning,root.RooFit.YVar(y,yBinning))
      rfHistTH2.Write()
      rfHistPdfTH2.Write()
      hist.Write()
      if channel.bakShape and compareHist != None:
        channel.bakShapeMkr.writeDebugHistsToCurrTDir(compareHist)
    else:
      plot = x.frame()
      rfHist.plotOn(plot)
      rfHistPdf.plotOn(plot)
      plot.SetName(hist.GetName())
      plot.Write()
      if channel.bakShape and writeBakShape:
        channel.bakShapeMkr.writeDebugHistsToCurrTDir()

  def write(self,outfilename,lumi,sumAllBak=True,sigInject=0.0):
    lumi *= 1000.0
    nuisance = self.nuisance

    ### ROOT Part
    ##########################################################
    outRootFilename = re.sub(r"\.txt",r".root",outfilename)
    outRootFile = root.TFile(outRootFilename, "RECREATE")
    outRootFile.cd()

    rootDebugString = ""

    observedN = {}

    for channel,channelName in zip(self.channels,self.channelNames):
        tmpDir = outRootFile.mkdir(channelName)
        tmpDir.cd()
        sumAllMCHist = None
        sumAllSigMCHist = None
        sumAllBakMCHist = None
        rootDebugString += "# channelName: {0}\n".format(channelName)
        if SIGGAUS:
          rootDebugString += channel.sigShapeMkr.debug
        for sigName in self.sigNames:
          tmpHist = channel.getSigHist(sigName).Clone(sigName)
          tmpHist.Scale(lumi)
          self.makeRFHistWrite(channel,tmpHist,tmpDir)
          #rootDebugString += "#     {0}: {1}\n".format(sigName,getIntegralAll(tmpHist))
          
          if sumAllSigMCHist == None:
            sumAllSigMCHist = tmpHist.Clone("sig")
          else:
            sumAllSigMCHist.Add(tmpHist)
          if sigInject != 0.0:
            tmpHist.Scale(sigInject)
            if sumAllMCHist == None:
              sumAllMCHist = tmpHist.Clone("data_obs")
            else:
              sumAllMCHist.Add(tmpHist)

          # Signal Shape Systematics
          for errHistKey in channel.sigErrHistsMap:
             tmpHist = channel.getSigSystHist(sigName,errHistKey).Clone(sigName+"_"+errHistKey)
             tmpHist.Scale(lumi)
             self.makeRFHistWrite(channel,tmpHist,tmpDir)
  
        if sumAllBak:
          #channel.dump()
          sumAllBakMCHist = channel.getBakHistTotal(lumi).Clone("bak")
          sumAllBakMCHist = shrinkTH1(sumAllBakMCHist,self.controlRegionLow[0],self.controlRegionHigh[1])
          #channel.bakShapeMkr.dump()
          if self.shape:
            rootDebugString += channel.bakShapeMkr.debug
            for nuName in channel.bakShapeMkr.errNames:
                origUp = channel.bakShapeMkr.errHists[nuName+"Up"]
                origDown = channel.bakShapeMkr.errHists[nuName+"Down"]
                tmpUp = origUp.Clone("bak_"+nuName+"Up")
                tmpDown = origDown.Clone("bak_"+nuName+"Down")
                self.makeRFHistWrite(channel,tmpUp,tmpDir)
                self.makeRFHistWrite(channel,tmpDown,tmpDir)

          # for simulated data_obs:
          sumAllBakMCHistReal = channel.bakHistTotalReal.Clone("sumAllBakMCHistReal")
          sumAllBakMCHistReal.Scale(lumi)
          sumAllBakMCHistReal = shrinkTH1(sumAllBakMCHistReal,self.controlRegionLow[0],self.controlRegionHigh[1])

          if self.toyData:
            if sumAllMCHist == None:
              sumAllMCHist = sumAllBakMCHist
            else:
              sumAllMCHist.Add(sumAllBakMCHist)
          else:
            if sumAllMCHist == None:
              sumAllMCHist = sumAllBakMCHistReal
              sumAllMCHist.SetName("data_obs")
            else:
              sumAllMCHist.Add(sumAllBakMCHistReal)
        else:
          for bakName in self.bakNames:
            tmpHist = channel.getBakHist(bakName).Clone(bakName)
            tmpHist.Scale(lumi)
            self.makeRFHistWrite(channel,tmpHist,tmpDir)
            #rootDebugString += "#     {0}: {1}\n".format(bakName,getIntegralAll(tmpHist,boundaries=massBounds))
            
            if sumAllMCHist == None:
                sumAllMCHist = tmpHist.Clone("data_obs")
            else:
                sumAllMCHist.Add(tmpHist)
        massLimits = [self.controlRegionLow[0],self.controlRegionHigh[1]]
        sumAllMCHist.Scale(int(getIntegralAll(sumAllMCHist,boundaries=massLimits))/getIntegralAll(sumAllMCHist,boundaries=massLimits)) # Make Integer
        if channel.datHistTotal == None:
          if self.toyData:
            print("Writing Toy Data Histogram")
            toy = sumAllMCHist.Clone("data_obs")
            toyHistogram(toy)
            observedN[channelName] = getIntegralAll(toy,boundaries=massLimits)
            self.makeRFHistWrite(channel,toy,tmpDir) #Pretend Toy Data
          else:
            print("Writing Pretend Data Histogram")
            observedN[channelName] = getIntegralAll(sumAllMCHist,boundaries=massLimits)
            self.makeRFHistWrite(channel,sumAllMCHist,tmpDir) #Pretend Data
        else:
          print("Writing Real Data Histogram")
          observedN[channelName] = getIntegralAll(channel.datHistTotal,boundaries=massLimits)
          self.makeRFHistWrite(channel,channel.datHistTotal,tmpDir) #Real Data
        self.makeRFHistWrite(channel,sumAllSigMCHist,tmpDir) #Pretend Signal
        self.makeRFHistWrite(channel,sumAllBakMCHist,tmpDir,compareHist=sumAllSigMCHist,writeBakShape=True) #Background Sum
        rootDebugString += "#######################\n"
        rootDebugString += "#     {0}\n".format(channelName)
        rootDebugString += "#     Pretend Obs: {0}\n".format(getIntegralAll(sumAllMCHist,boundaries=massLimits))
        rootDebugString += "#     All Signal:  {0}\n".format(getIntegralAll(sumAllSigMCHist,boundaries=massLimits))
        rootDebugString += "#     All Bak:     {0}\n".format(getIntegralAll(sumAllBakMCHist,boundaries=massLimits))
        #rootDebugString += "#     Pretend Obs: {0}\n".format(getIntegralAll(sumAllMCHist))
        #rootDebugString += "#     All Signal:  {0}\n".format(getIntegralAll(sumAllSigMCHist))
        #rootDebugString += "#     All Bak:     {0}\n".format(getIntegralAll(sumAllBakMCHist))

    outRootFile.Close()

    ### Text Part
    ##########################################################

    print("Writing Card: {0} & {1}".format(outfilename,outRootFilename))
    outfile = open(outfilename,"w")
    outfile.write("# Hmumu shape combine datacard produced by makeTables.py\n")
    now = datetime.datetime.now().replace(microsecond=0).isoformat(' ')
    outfile.write("# {0}\n".format(now))
    outfile.write("############################### \n")
    outfile.write("############################### \n")
    outfile.write("imax {0}\n".format(len(self.channels)))
    #outfile.write("jmax {0}\n".format(len(backgroundNames)))
    outfile.write("jmax {0}\n".format("*"))
    outfile.write("kmax {0}\n".format("*"))
    outfile.write("------------\n")
    outfile.write("shapes * * {0} $CHANNEL/$PROCESS $CHANNEL/$PROCESS_$SYSTEMATIC\n".format( os.path.basename(outRootFilename)))
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
      observedNumber = observedN[channelName]
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
        iProc = -len(channel.sigNames)+1
        for sigName in self.sigNames:
          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
  
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append(sigName)
  
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)
  
          expNum = channel.getSigXSec(sigName)*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
  
          iParam += 1
          iProc += 1

        if sumAllBak:

          binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          binFormatList.append(channelName)
    
          proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc1FormatList.append("bak")
    
          proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
          proc2FormatList.append(iProc)

          expNum = channel.getBakXSecTotal()*lumi
          decimals = ".4f"
          if expNum>1000.0:
            decimals = ".4e"
          rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
          rateFormatList.append(expNum)
      
          iParam += 1
          iProc += 1
        else:
          for bakName in self.bakNames:
            binFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            binFormatList.append(channelName)
    
            proc1FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            proc1FormatList.append(bakName)
    
            proc2FormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            proc2FormatList.append(iProc)
    
            expNum = channel.getBakXSec(bakName)*lumi
            decimals = ".4f"
            if expNum>1000.0:
              decimals = ".4e"
            rateFormatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+decimals+"} "
            rateFormatList.append(expNum)
    
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

    for nu in nuisance:
      thisNu = nuisance[nu]
      formatString = "{0:<8} {1:^4} "
      formatList = [nu,"lnN"]
      iParam = 2
      for channel,channelName in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            if thisNu.has_key(sigName):
              value = thisNu[sigName]+1.0
            formatList.append(value)
            iParam += 1
          if sumAllBak:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              if thisNu.has_key(bakName):
                value = thisNu[bakName]+1.0
              formatList.append(value)
              iParam += 1
          else:
            for bakName in self.bakNames:
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              if thisNu.has_key(bakName):
                value = thisNu[bakName]+1.0
              formatList.append(value)
              iParam += 1
      formatString += "\n"
      #print formatString
      #print formatList
      outfile.write(formatString.format(*formatList))

    # Sig Shape Uncertainties (All Correlated)
    for channel,channelName in zip(self.channels,self.channelNames):
      for nuisanceName in channel.sigErrNames:
        formatString = "{0:<8} {1:^4} "
        formatList = [nuisanceName,"shape"]
        iParam = 2
        for channel2,channelName2 in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "1"
            formatList.append(value)
            iParam += 1
          if sumAllBak:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              formatList.append(value)
              iParam += 1
          else:
            for bakName in self.bakNames:
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              formatList.append(value)
              iParam += 1
        formatString += "\n"
        #print formatString
        #print formatList
        outfile.write(formatString.format(*formatList))
      break

    # Bak Shape Uncertainties (All Correlated)
    for channel,channelName in zip(self.channels,self.channelNames):
      for nuisanceName in channel.bakShapeMkr.errNames:
        formatString = "{0:<8} {1:^4} "
        formatList = [nuisanceName,"shape"]
        iParam = 2
        for channel2,channelName2 in zip(self.channels,self.channelNames):
          for sigName in self.sigNames:
            formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
            value = "-"
            formatList.append(value)
            iParam += 1
          if sumAllBak:
              bakName="bak"
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "1"
              formatList.append(value)
              iParam += 1
          else:
            for bakName in self.bakNames:
              formatString += "{"+str(iParam)+":^"+str(self.largestChannelName)+"} "
              value = "-"
              formatList.append(value)
              iParam += 1
        formatString += "\n"
        #print formatString
        #print formatList
        outfile.write(formatString.format(*formatList))
      break


    #Debugging
    outfile.write("#################################\n")
    for channel,channelName in zip(self.channels,self.channelNames):
        outfile.write("#\n")
        outfile.write("#info: channel {0}: \n".format(channelName))
        outfile.write("#  x var name: {0} \n".format(channel.x.GetName()))
        outfile.write("#  x var range: [{0:.3g},{1:.3g}] \n".format(channel.x.getMin(),channel.x.getMax()))
        if channel.is2D:
          outfile.write("#  y var name: {0} \n".format(channel.y.GetName()))
          outfile.write("#  y var range: [{0:.3g},{1:.3g}] \n".format(channel.y.getMin(),channel.y.getMax()))
    outfile.write(rootDebugString)
    outfile.close()

class ThreadedCardMaker(myThread):
  def __init__(self,*args,**dictArgs):
    myThread.__init__(self)
    self.writeArgs = (dictArgs["outfilename"],dictArgs["lumi"])
    self.writeArgsDict = {}
    if dictArgs.has_key("sumAllBak"):
        self.writeArgsDict["sumAllBak"] = dictArgs["sumAllBak"]
    if dictArgs.has_key("sigInject"):
        self.writeArgsDict["sigInject"] = dictArgs["sigInject"]
    self.shapeDataCardMaker = True
    if dictArgs.has_key("shapeDataCardMaker"):
        self.shapeDataCardMaker = dictArgs["shapeDataCardMaker"]
        dictArgs.pop("shapeDataCardMaker",None)
    self.args = args
    dictArgs.pop("sumAllBak",None)
    dictArgs.pop("sigInject",None)
    dictArgs.pop("outfilename",None)
    dictArgs.pop("lumi",None)
    self.dictArgs = dictArgs
    self.started = False
  def run(self):
    #try:
      self.started = True
      dataCardMassShape = None
      if self.shapeDataCardMaker:
        dataCardMassShape = ShapeDataCardMaker(*(self.args),**(self.dictArgs))
      else:
        dataCardMassShape = DataCardMaker(*(self.args),**(self.dictArgs))
      dataCardMassShape.write(*(self.writeArgs),**(self.writeArgsDict))
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
  root.gROOT.SetBatch(True)

  directory = "input/"
  outDir = "statsCards/"
  periods = ["7TeV","8TeV"]
  analysesInc = ["IncPresel","IncBDTCut"]
  analysesVBF = ["VBFPresel","VBFBDTCut"]
  analyses = analysesInc + analysesVBF
  categoriesInc = ["BB","BO","BE","OO","OE","EE"]
  categoriesVBF = ["BB","NotBB"]
  tmpList = []
  for a in analysesInc:
    for c in categoriesInc:
        tmpList.append(a+c)
  #analyses += tmpList
  tmpList = []
  for a in analysesVBF:
    for c in categoriesVBF:
        tmpList.append(a+c)
  #analyses += tmpList
  combinations = []
  combinationsLong = []
  """
  combinations.append((
        ["IncBDTCut"+x for x in categoriesInc],"IncBDTCutCat"
  ))
  combinations.append((
        ["VBFBDTCut"+x for x in categoriesVBF],"VBFBDTCutCat"
  ))
  combinations.append((
        ["IncPresel"+x for x in categoriesInc],"IncPreselCat"
  ))
  combinations.append((
        ["VBFPresel"+x for x in categoriesVBF],"VBFPreselCat"
  ))
  combinations.append((
        ["IncBDTCut","VBFBDTCut"],"BDTCut"
  ))
  combinations.append((
        ["IncPresel","VBFPresel"],"Presel"
  ))
  combinations.append((
        ["VBFPresel"+x for x in categoriesVBF]+["IncPresel"+x for x in categoriesInc],"PreselCat"
  ))
  combinationsLong.append((
        ["VBFBDTCut"+x for x in categoriesVBF]+["IncBDTCut"+x for x in categoriesInc],"BDTCutCat"
  ))
  """
  histPostFix="/mDiMu"
  #analyses = ["mDiMu"]
  #histPostFix=""
  signalNames=["ggHmumu125","vbfHmumu125","wHmumu125","zHmumu125"]
  backgroundNames= ["DYJetsToLL","ttbar"]
  dataDict = {}
  dataDict["8TeV"] = [
    #"SingleMuRun2012Av1",
    #"SingleMuRun2012Bv1",
    #"SingleMuRun2012Cv1",
    #"SingleMuRun2012Cv2"
  ]
  dataDict["7TeV"] = [
    #"SingleMuRun2011Av1",
    #"SingleMuRun2011Bv1"
  ]
  lumiListLong = [5,10,15,20,25,30,40,50,75,100,200,500,1000,2000,5000]
  lumiListLong = [20,30,50,100,500,1000,5000]
  lumiList = [lumiDict["8TeV"],20,25,30]
  lumiList = [20]
  #lumiListLong = lumiList

  MassRebin = 1 # 4 Bins per GeV originally
  controlRegionVeryLow=[80,110]
  controlRegionLow=[110,120]
  controlRegionHigh=[130,160]

  shape=True
  toyData=args.toyData

  print("Simple Analyses to run:")
  for a in analyses:
    print("  {0}".format(a))
  print("Combination Analyses to run:")
  for c in combinations:
    print("  {0}".format(c[1]))
    for a in c[0]:
      print("    {0}".format(a))

  print("Creating Threads...")
  threads = []
  for p in periods:
    for i in lumiList:
      if p == "7TeV":
        i = lumiDict[p]
      for ana in analyses:
        tmp = ThreadedCardMaker(
          #__init__ args:
          directory,[ana],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin],bakShape=shape,
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+ana+"_"+p+"_"+str(i)+".txt",lumi=i
          )
        threads.append(tmp)
      for comb in combinations:
       threads.append(
        ThreadedCardMaker(
          #__init__ args:
          directory,
          comb[0],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin], bakShape=shape,
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+comb[1]+"_"+p+"_"+str(i)+".txt",lumi=i
        )
       )
#      ## Cut and Count!!!
#      for ana in analyses:
#        tmp = ThreadedCardMaker(
#          #__init__ args:
#          directory,[ana],
#          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
#          rebin=[MassRebin],bakShape=shape,
#          controlRegionLow=[123.0,125.0],controlRegionHigh=[125.0,127.0],histNameSuffix=histPostFix,
#          controlRegionVeryLow=controlRegionVeryLow,nuisanceMap=nuisanceMap,
#          #write args:
#          outfilename=outDir+"CNC_"+ana+"_"+p+"_"+str(i)+".txt",lumi=i,
#  
#          shapeDataCardMaker=False
#          )
#        threads.append(tmp)
      if p == "7TeV":
        break

  for p in periods:
    for i in lumiListLong:
      if p == "7TeV":
        i = lumiDict[p]
      for comb in combinationsLong:
       threads.append(
        ThreadedCardMaker(
          #__init__ args:
          directory,
          comb[0],
          appendPeriod(signalNames,p),appendPeriod(backgroundNames,p),dataDict[p],
          rebin=[MassRebin], bakShape=shape,
          controlRegionLow=controlRegionLow,controlRegionHigh=controlRegionHigh,histNameSuffix=histPostFix,
          controlRegionVeryLow=controlRegionVeryLow,toyData=toyData,nuisanceMap=nuisanceMap,sigInject=args.signalInject,
          #write args:
          outfilename=outDir+comb[1]+"_"+p+"_"+str(i)+".txt",lumi=i
        )
       )
      if p == "7TeV":
        break

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
      

  runFile = open(outDir+"run.sh","w")
  batchString = \
"""#!/bin/bash

chmod +x lxbatch.sh

for i in *.txt; do
    [[ -e "$i" ]] || continue
echo "Running on "$i
bsub lxbatch.sh $i
#bsub -q 1nh lxbatch.sh $i
done
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"lxbatch.sh","w")
  batchString = \
"""#!/bin/bash
echo "Sourcing cmsset_default.sh"
cd /afs/cern.ch/cms/sw
source cmsset_default.sh
export SCRAM_ARCH=slc5_amd64_gcc462
echo "SCRAM_ARCH is $SCRAM_ARCH"
cd $LS_SUBCWD
echo "In Directory: "
pwd
eval `scramv1 runtime -sh`
echo "cmsenv success!"
date

TXTSUFFIX=".txt"
FILENAME=$1
DIRNAME="Dir"$1"Dir"
ROOTFILENAME=${1%$TXTSUFFIX}.root

mkdir $DIRNAME
cp $FILENAME $DIRNAME/
cp $ROOTFILENAME $DIRNAME/
cd $DIRNAME

echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out

cp $FILENAME.out ..


echo "done"
date
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"notlxbatch.sh","w")
  batchString = \
"""#!/bin/bash
echo "running notlxbatch.sh"
date
for i in *.txt; do
    [[ -e "$i" ]] || continue
FILENAME=$i
echo "executing combine -M Asymptotic $FILENAME >& $FILENAME.out"

combine -M Asymptotic $FILENAME >& $FILENAME.out
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M ProfileLikelihood --signifcance $FILENAME -t 100 --expectSignal=1 >& $FILENAME.sig"

combine -M ProfileLikelihood --significance $FILENAME -t 100 --expectSignal=1 >& $FILENAME.expsig
rm -f roostats*
rm -f higgsCombineTest*.root

echo "executing combine -M MaxLikelihoodFit $FILENAME >& $FILENAME.sig"

combine -M MaxLikelihoodFit $FILENAME >& $FILENAME.mu
rm -f roostats*
rm -f higgsCombineTest*.root

done

date
echo "done"
"""
  runFile.write(batchString)
  runFile.close()

  runFile = open(outDir+"getStatus.sh","w")
  batchString = \
"""#!/bin/bash

echo "==========================="
echo "Files Found: `ls *.out | wc -l` of `ls *.txt | wc -l`"
echo "==========================="
for i in *.out; do wc $i; done
echo "==========================="
"""
  runFile.write(batchString)
  runFile.close()
