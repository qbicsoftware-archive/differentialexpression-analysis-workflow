#!/bin/bash
module load qbic/anaconda
module load qbic/r/qbic-r-3.2.2


workflowDir=$(cat wfdir)
cp 'GROUPS' $workflowDir"/etc"

#parse using CTDopts and run workflow
Rscript DESeq2.4.R
cp wfdir wfdir2