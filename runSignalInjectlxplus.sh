#!/bin/bash

SIGNALSTRENGTH=1.0

REMOTEDIR=/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_5_2_5/stats/
if ! `ssh lxplus touch $REMOTEDIR/touchfile.txt`; then
  echo "Error: Remote directory $REMOTEDIR doesn't exist or isn't writable"
  exit
fi

rm -f shapes/*
rm -f statsCards/*
rm -f statsInput/*
rm -f statsOutput/*

nice ./makeCards.py --signalInject $SIGNALSTRENGTH --toyData
echo "Removing files in lxplus:$REMOTEDIR"
ssh lxplus "cd /tmp/jhugon/; rm -rf $REMOTEDIR/*;echo \"Contents of dir: \`ls $REMOTEDIR \`\""
echo "Copying input files to lxplus..."
scp statsCards/* lxplus:/afs/cern.ch/user/j/jhugon/work/private/stats/CMSSW_5_2_5/stats/.
echo "Running combine on lxplus..."
ssh lxplus "cd $REMOTEDIR; eval \`scramv1 runtime -sh\`;nice bash notlxbatch.sh;"
echo "Copying output files from lxplus..."
scp lxplus:$REMOTEDIR/*.out statsInput/.
scp lxplus:$REMOTEDIR/*.sig statsInput/.
scp lxplus:$REMOTEDIR/*.expsig statsInput/.
scp lxplus:$REMOTEDIR/*.mu statsInput/.

nice ./makeShapePlots.py
nice ./makeLimitPlots.py
nice ./makeSigMuPlots.py