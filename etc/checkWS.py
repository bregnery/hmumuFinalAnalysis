#!/usr/bin/env python
from sys import argv, stdout, stderr, exit
import ROOT
ROOT.gROOT.SetBatch(True)
ROOT.gSystem.Load("libHiggsAnalysisCombinedLimit.so")
root = ROOT

def rooArgSet2List(x):
  itr = x.createIterator()
  result = []
  while True:
    ele = itr.Next()
    if ele:
      result.append(ele)
    else:
      break
  return result

if len(argv) < 2:
  exit(0)

f = root.TFile(argv[1])
w = f.Get("w")
wName = "w"
if not w:
  for k in f.GetListOfKeys():
    if k.GetClassName() == "RooWorkspace":
      w = k.ReadObj()
      wName = k.GetName()
      break

w.Print()
print "Variables:"
for var in rooArgSet2List(w.allVars()):
  var.Print()
print "#"*40
