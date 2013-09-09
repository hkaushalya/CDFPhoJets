#!/bin/tcsh
#########################################
#need run this in the Stntuple base dir.
# change the nEVts, and dataset below.
##########################################

source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 6.1.6.m
setenv PDFTOOLS_DIR ./samantha/pdftools/
setenv DATASET $1
setenv STNTUPLE_CATALOG http://www-cdf.fnal.gov/~cdfopr/Stntuple/cafdfc

echo "DATASET = $1"
echo "Nevts   = $2"
echo "fileset = $3"
echo "cafind  = $4"
#.x samantha/pdftools/pdf_syst.C+(1,1,0)
root -b <<+
.x samantha/pdftools/pdf_syst.C+($2,$3,0,$4)
+
#
#rm -rf *.so
