#!/bin/bash

# define the job array (104 jobs)
#BSUB -J testArray[1-104]
#BSUB -q 8nh
#BSUB -o "%J.%I.out"
#BSUB -e "%J.%I.err"

#BSUB -n 1


# setup the working environment
source $HOME/setup/lxplus_setup.sh

# getting to the working dir
cd $TMPDIR
mkdir batch
cd batch


cp /afs/cern.ch/work/j/jloyal/data/Hbb/zaidan_ntuples/nominal_only.py .
cp /afs/cern.ch/work/j/jloyal/data/Hbb/zaidan_ntuples/filelist.txt .

# determine the files to run over
readarray FILES < filelist.txt
./nominal_only.py -i ${FILES[$LSB_JOBINDEX-1]} -o output.MC.OneLeptonMuons.nom.$LSB_JOBINDEX.$LSB_JOBID.root

# copy the output back
MYBATCHDIR=/afs/cern.ch/work/j/jloyal/batch/Egamma/
cp *.root $MYBATCHDIR
