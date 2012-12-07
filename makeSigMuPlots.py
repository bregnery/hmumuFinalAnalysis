#!/usr/bin/env python

from helpers import *
import ROOT as root
import glob
import re
import os.path

import matplotlib.pyplot as mpl
import numpy

from xsec import *

#######################################

#~48 Charactars Max
titleMap = {
  "AllCat":"H->#mu#mu Catagories Combination",
  "IncCat":"Inclusive Categories Comb.",
  "VBFCat":"VBF Categories Comb.",

  "IncPresel":"Inclusive Preselection",
  "VBFPresel":"VBF Preselection",

  "Pt0to30":"p_{T}^{#mu#mu} #in [0,30]",
  "Pt30to50":"p_{T}^{#mu#mu} #in [30,50]",
  "Pt50to125":"p_{T}^{#mu#mu} #in [50,125]",
  "Pt125to250":"p_{T}^{#mu#mu} #in [125,250]",
  "Pt250":"p_{T}^{#mu#mu}>250",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"H->#mu#mu BDT Combination",
  "IncBDTCut":"Inclusive BDT ",
  "VBFBDTCut":"VBF BDT ",

  "BDTCutCat":"BDT Res. Cat. Combination",
  "IncBDTCutCat":"Inclusive BDT Resolution Categories",
  "VBFBDTCutCat":"VBF BDT Resolution Categories",

  "PreselCat":"Res. Cat. Preselection Combination",
  "IncPreselCat":"Inclusive Resolution Cat. Preselection",
  "VBFPreselCat":"VBF Cat. Resolution Preselection",

  "IncBDTCutBB":"Inclusive BDT BB",
  "IncBDTCutBO":"Inclusive BDT BO",
  "IncBDTCutBE":"Inclusive BDT BE",
  "IncBDTCutOO":"Inclusive BDT OO",
  "IncBDTCutOE":"Inclusive BDT OE",
  "IncBDTCutEE":"Inclusive BDT EE",
  "IncBDTCutNotBB":"Inclusive BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"Inclusive Preselection BB",
  "IncPreselBO":"Inclusive Preselection BO",
  "IncPreselBE":"Inclusive Preselection BE",
  "IncPreselOO":"Inclusive Preselection OO",
  "IncPreselOE":"Inclusive Preselection OE",
  "IncPreselEE":"Inclusive Preselection EE",
  "IncPreselNotBB":"Inclusive Preselection !BB",
  "VBFPreselBB":"VBF Preselection BB",
  "VBFPreselNotBB":"VBF Preselection !BB"
}

comparisonMap = {
  "AllCat":"All Cat. Comb.",
  "IncCat":"Inc. Cat. Comb.",
  "VBFCat":"VBF Cat. Comb.",

  "Presel":"Presel. Comb.",
  "IncPresel":"Inc. Presel.",
  "VBFPresel":"VBF Presel.",

  "Pt0to30":"$p_{T}^{\mu\mu} \in [0,30]$",
  "Pt30to50":"$p_{T}^{\mu\mu} \in [30,50]$",
  "Pt50to125":"$p_{T}^{\mu\mu} \in [50,125]$",
  "Pt125to250":"$p_{T}^{\mu\mu} \in [125,250]$",
  "Pt250":"$p_{T}^{\mu\mu}>250$",

  "VBFLoose":"VBFL",
  "VBFMedium":"VBFM",
  "VBFTight":"VBFT",
  "VBFVeryTight":"VBFVT",

  "BDTCut":"BDT Comb.",
  "IncBDTCut":"Inc. BDT",
  "VBFBDTCut":"VBF BDT",

  "BDTCutCat":"BDT Res. Comb.",
  "IncBDTCutCat":"Inc. BDT Res.",
  "VBFBDTCutCat":"VBF BDT Res.",

  #"BDTCutCat":"Combination",
  #"IncBDTCutCat":"Inclusive",
  #"VBFBDTCutCat":"VBF",

  "PreselCat":"Presel. Res. Comb.",
  "IncPreselCat":"Inc. Res. Presel.",
  "VBFPreselCat":"VBF Res. Presel.",

  "IncBDTCutBB":"Inc. BDT BB",
  "IncBDTCutBO":"Inc. BDT BO",
  "IncBDTCutBE":"Inc. BDT BE",
  "IncBDTCutOO":"Inc. BDT OO",
  "IncBDTCutOE":"Inc. BDT OE",
  "IncBDTCutEE":"Inc. BDT EE",
  "IncBDTCutNotBB":"Inc. BDT !BB",
  "VBFBDTCutBB":"VBF BDT BB",
  "VBFBDTCutNotBB":"VBF BDT !BB",
  "IncPreselBB":"Inc. Presel. BB",
  "IncPreselBO":"Inc. Presel. BO",
  "IncPreselBE":"Inc. Presel. BE",
  "IncPreselOO":"Inc. Presel. OO",
  "IncPreselOE":"Inc. Presel. OE",
  "IncPreselEE":"Inc. Presel. EE",
  "IncPreselNotBB":"Inc. Presel. !BB",
  "VBFPreselBB":"VBF Presel. BB",
  "VBFPreselNotBB":"VBF Presel. !BB"
}

#######################################

# Do match and don't match w/o extensions.  \.expsig and \.sig are added automatically
def getDataSig(fileString,matchString=r"_([-\d.]+)\.txt",dontMatchStrings=[],doSort=True):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  result = []
  fNames =  glob.glob(fileString)
  if doSort:
    fNames.sort(key=sortfun)
  #print fileString
  #print fNames
  for fname in fNames: 
    dontMatch = False
    for dont in dontMatchStrings:
      if re.search(dont,fname):
        dontMatch = True
    if dontMatch:
        continue
    tmpF = open(fname)
    match = re.search(matchString+r"\.sig",fname)
    obs = -10.0
    exp = -10.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
    for line in tmpF:
      obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
    expFname = os.path.splitext(fname)[0]+".expsig"
    try:
      tmpFExp = open(expFname)
      for line in tmpFExp:
        obsMatch = re.search(r"^Significance:[\s]+([.\deE]+)",line)
        if obsMatch:
          exp = obsMatch.group(1)
    except Exception:
      print("Expected Limit Not Found: "+expFname)
    thisPoint = [xNum,obs,exp]
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result

def getDataMu(fileString,matchString=r"_([-\d.]+)\.txt\.mu",dontMatchStrings=[],doSort=True):
  def sortfun(s):
    match = re.search(matchString,s)
    result = 1e12
    if match:
      result = float(match.group(1))
    return result

  result = []
  fNames =  glob.glob(fileString)
  if doSort:
    fNames.sort(key=sortfun)
  #print fileString
  #print fNames
  for fname in fNames: 
    dontMatch = False
    for dont in dontMatchStrings:
      if re.search(dont,fname):
        dontMatch = True
    if dontMatch:
        continue
    tmpF = open(fname)
    match = re.search(matchString,fname)
    obs = -10.0
    low1sig = -10.0
    high1sig = -10.0
    xNum = -10.0
    if match:
      xNum = match.group(1)
    for line in tmpF:
      #obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+([.\deE-]+)/+([.\deE-]+)",line)
      obsMatch = re.search(r"Best fit r:[\s]+([.\deE-]+)[\s]+\-([.\deE-]+)/\+([.\deE-]+)",line)
      if obsMatch:
        obs = obsMatch.group(1)
        low1sig = obsMatch.group(2)
        high1sig = obsMatch.group(3)
        #print("original line: {0}\n  Extracted: {1} -{2} +{3}".format(obsMatch.group(0),obs,low1sig,high1sig))
    thisPoint = [xNum,obs,low1sig,high1sig]
    if thisPoint.count("-10.0")>0:
        continue
    if thisPoint.count(-10.0)>0:
        continue
    #print thisPoint
    result.append(thisPoint)
  return result


class RelativePlot:
  def __init__(self,dataPoints, canvas, legend, caption, ylabel="Error on #sigma/#sigma_{SM}", xlabel="Integrated Luminosity [fb^{-1}]",caption2="",caption3="",ylimits=[],xlimits=[],vertLines=[],showObs=False,energyStr="8TeV"):
    errGraph = root.TGraph()
    #errGraph.SetLineStyle(2)
    errGraph.SetLineColor(root.kRed)
    self.errGraph = errGraph
    obsGraph = root.TGraph()
    obsGraph.SetLineColor(root.kRed)
    #obsGraph.SetLineStyle(2)
    self.obsGraph = obsGraph
    iPoint = 0
    ymax = 1e-20
    ymin = 1e20
    for point in dataPoints:
      value = None
      obs = 0.0
      xNum = float(point[0])
      if len(xlimits)==2:
        if xNum<xlimits[0]:
            continue
        if xNum>xlimits[1]:
            continue
      if len(point)==4: #mu
        #thisPoint = [xNum,obs,low1sig,high1sig]
        high1sig = float(point[3])
        low1sig = float(point[2])
        if high1sig > low1sig:
            value = high1sig
        else:
            value = low1sig
      else: #sig len=3
        #thisPoint = [xNum,obs,exp]
        value = float(point[2])
        obs = float(point[1])
      errGraph.SetPoint(iPoint,xNum,value)
      obsGraph.SetPoint(iPoint,xNum,obs)
      iPoint += 1

    self.vertLines = []
    self.vertLabel = []
    for xPos in vertLines:
      tmp = root.TGraph()
      tmp.SetPoint(0,xPos,ymin)
      tmp.SetPoint(1,xPos,ymax)
      tmp.SetLineColor(root.kRed)
      tmp.SetLineStyle(3)
      self.vertLines.append(tmp)
      self.vertLabel.append(xPos)

    label = root.TLatex()
    #label.SetNDC()
    label.SetTextFont(root.gStyle.GetLabelFont("X"))
    label.SetTextSize(root.gStyle.GetLabelSize("X"))
    label.SetTextColor(root.kRed)
    label.SetTextAlign(22)
    self.label=label

    errGraph.Draw("al")
    errGraph.GetXaxis().SetTitle(xlabel)
    errGraph.GetYaxis().SetTitle(ylabel)
    if len(ylimits)==2:
        errGraph.GetYaxis().SetRangeUser(*ylimits)
    if showObs:
      obsGraph.Draw("l")

    tlatex = root.TLatex()
    tlatex.SetNDC()
    tlatex.SetTextFont(root.gStyle.GetLabelFont())
    #tlatex.SetTextSize(0.05)
    tlatex.SetTextSize(0.04)
    tlatex.SetTextAlign(12)
    tlatex.DrawLatex(gStyle.GetPadLeftMargin(),0.96,"CMS Internal")
    tlatex.SetTextAlign(32)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin(),0.96,caption)

    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.88,caption2)
    tlatex.DrawLatex(1.0-gStyle.GetPadRightMargin()-0.03,0.82,caption3)

    for g in self.vertLines:
      g.Draw("l")

    canvas.RedrawAxis()

class ComparePlot:
  def __init__(self,data,ylabel="Expected 95% CL Limit $\sigma/\sigma_{SM}$",titleMap={},showObs=False):
    fig = mpl.figure()
    self.fig = fig
    ax1 = fig.add_subplot(111)
    self.ax1 = ax1
    ax1bounds = ax1.get_position().bounds
    ax1.set_position([0.25,0.1,0.7,0.85])

    data.sort(key=lambda x: x[0].lower())

    obs = [float(point[1]) for point in data]
    exp = None
    high1sigs = None
    low1sigs = None
    if len(data[0])==4: #mu
      #thisPoint = [xNum,obs,low1sig,high1sig]
      high1sigs = [float(point[3]) for point in data]
      low1sigs = [float(point[2]) for point in data]
    else: #sig len=3
      #thisPoint = [xNum,obs,exp]
      exp = [float(point[2]) for point in data]
      high1sigs = [0.0 for point in data]

    xPos = numpy.arange(len(obs))
    xLabels = [point[0] for point in data]
    xLabels = [re.sub(r".*/","",s) for s in xLabels]
    xLabels = [titleMap[s] if titleMap.has_key(s) else s for s in xLabels]

    ax1.set_yticks(xPos+0.25)
    ax1.set_yticklabels(tuple(xLabels),size="small")
    ax1.set_xlabel(ylabel)
    #ax1.set_xlim([0,20])
    bars = None
    if len(data[0])==4: #mu
      bars = ax1.barh(xPos,obs, 0.5, xerr=[low1sigs,high1sigs],ecolor="k")
    else: #sig
      bars = ax1.barh(xPos,exp, 0.5)
    self.bars = bars
    xPosObs = [x+0.25 for x in xPos]
    if showObs and len(data[0]) != 4: #is sig
      self.obs = ax1.plot(obs,xPosObs,marker="|",color="r",markersize=10,linestyle="None")
    writeInValues = getattr(self,"writeInValues")
    writeInValues(bars,high1sigs)
  def save(self,saveName):
    self.fig.savefig(saveName+".png")
    self.fig.savefig(saveName+".pdf")
  def writeInValues(self,rects,high2sigs):
    axisWidth = self.ax1.get_xlim()[1]-self.ax1.get_xlim()[0]
    size="medium"
    if len(rects)<5:
      size="large"
    elif len(rects)>12:
      size="xx-small"
    elif len(rects)>9:
      size="x-small"
    maxWidth = axisWidth
    for rect,sigUp in zip(rects,high2sigs):
      width = rect.get_width()
      if (width/maxWidth < 0.1): # The bars aren't wide enough to print the ranking inside
        xloc = width + sigUp + maxWidth*0.01 # Shift the text to the right side of the right edge
        clr = 'black' # Black against white background
        align = 'left'
      else:
        xloc = width - maxWidth*0.01 # Shift the text to the left side of the right edge
        clr = 'white' # White on magenta
        align = 'right'
      yloc = rect.get_y()+rect.get_height()/2.0 
      valueStr = "{0:.1f}".format(width)
      self.ax1.text(xloc, yloc, valueStr, horizontalalignment=align,
            verticalalignment='center', color=clr, weight='bold',size=size)

if __name__ == "__main__":

  dirName = "statsInput/"
  
  outDir = "statsOutput/"
  
  root.gErrorIgnoreLevel = root.kWarning
  root.gROOT.SetBatch(True)
  setStyle()
  canvas = root.TCanvas()
  canvas.SetLogx(1)
  canvas.SetLogy(0)
  
  mpl.rcParams["font.family"] = "sans-serif"
  #print mpl.rcParams["backend"]

  ylimits=[]

  lumisToUse={"7TeV":lumiDict["7TeV"],"8TeV":20}
  
  for period in ["7TeV","8TeV","14TeV"]:
    fnToGlob = dirName+"*_"+period+"_*.txt.out"
    allfiles = glob.glob(fnToGlob)

    ## Limit v. Lumi
    energyStr = ""
    plots = set()
    for fn in allfiles:
      match = re.search(r".*/(.+)_(.+)_[.\d]+.txt.out",fn)
      badPlot = re.search(r"Silly",fn)
      badPlot2 = re.search(r"Silly",fn)
      if match and not (badPlot or badPlot2):
        plots.add(match.group(1))
        energyStr = match.group(2)
  
    caption2 = "#sqrt{s}="+energyStr
    legend = root.TLegend(0.58,0.70,0.9,0.9)
    legend.SetFillColor(0)
    legend.SetLineColor(0)
    
    for ytitle,fnPref in zip(["Significance","$\sigma_{Measured}/\sigma_{SM}$"],["sig","mu"]):
      for plotName in plots:
        data = None
        ylabel = None
        if fnPref == "sig":
          data = getDataSig(dirName+plotName+"_"+energyStr+"_*.txt*")
          ylabel="Expected Significance"
        else:
          data = getDataMu(dirName+plotName+"_"+energyStr+"_*.txt*")
          ylabel="Error on #sigma/#sigma_{SM}"
        #print("{0} {1} v. Lumi: {2}".format(period,fnPref, data))
        if len(data)<=1:
          continue
        title = titleMap[plotName]
        title = "Standard Model H#rightarrow#mu#mu"
        incPlot = RelativePlot(data,canvas,legend,title,caption2=caption2,ylabel=ylabel,energyStr=energyStr)
        saveAs(canvas,outDir+fnPref+plotName+"_"+energyStr)
    
    ## Compare all types of limits
    if period == "14TeV":
        continue
    #veto = [r"CNC",r"PM","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","BDT","BB","BO","BE","OO","OE","EE","NotBB"]
    #veto = [r"CNC",r"PM","Presel"]
    #veto = [r"CNC",r"PM","BDT"]
    veto = []

    mustBe = r"(.+)_(.+)_[.\d]+.txt"
    #mustBe = r"(.+Cat)_(.+)_[.\d]+.txt"

    desiredLumiStr=str(lumisToUse[period])
    fnGlobStr = dirName+"*_"+period+"_"+desiredLumiStr+".txt*"
    compareDataSig = getDataSig(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    compareDataMu = getDataMu(fnGlobStr,matchString=mustBe,dontMatchStrings=veto,doSort=False)
    #print compareDataSig
    #print compareDataMu

    for datCase,ytitle,fnPref in zip([compareDataSig,compareDataMu],["Expected Significance","$\sigma/\sigma_{SM}$"],["sig","mu"]):
      if len(datCase)==0:
        print("No Data to Compare for {0} {1}!!".format(period,fnPref))
        continue
      comparePlot = ComparePlot(datCase,titleMap=comparisonMap,showObs=False,ylabel=ytitle)
      comparePlot.fig.text(0.9,0.2,"$\mathcal{L}="+desiredLumiStr+"$ fb$^{-1}$",horizontalalignment="right",size="x-large")
      comparePlot.fig.text(0.9,0.27,"$\sqrt{s}=$"+energyStr,horizontalalignment="right",size="x-large")
      comparePlot.save(outDir+fnPref+"Compare"+"_"+energyStr)