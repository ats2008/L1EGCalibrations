#!/bin/bash
set -x
source /cvmfs/cms.cern.ch/cmsset_default.sh 
export HOME=@@HOME
export X509_USER_PROXY=@@PROXY
cd @@DIRNAME
eval `scramv1 runtime -sh`
TMPDIR=`mktemp -d`
cd $TMPDIR
cp  @@PWD/applyIsolation.exe .
mv @@RUNSCRIPT @@RUNSCRIPT.busy
./applyIsolation.exe @@CFGFILENAME
if [ $? -eq 0 ]; then 
    mv * @@DIRNAME
    mv @@RUNSCRIPT.busy @@RUNSCRIPT.sucess
    echo SUCESS
else
    mv @@RUNSCRIPT.busy @@RUNSCRIPT
    echo FAIL
fi
