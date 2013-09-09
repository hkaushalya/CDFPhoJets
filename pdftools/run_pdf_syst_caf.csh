#!/bin/tcsh
#########################################
#need run this in the Stntuple base dir.
# change the nEVts, and dataset below.
##########################################

source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 6.1.4
setenv PDFTOOLS_DIR ./samantha/pdftools/
setenv STNTUPLE_CATALOG http://www-cdf.fnal.gov/~cdfopr/Stntuple/cafdfc

echo "Nevts   = $1"
echo "pdfset  = $2"
echo "cafind  = $3"
echo "cafsections  = $4"
root -b <<+
.x samantha/pdftools/pdf_syst_caf.C+($1,$2,$3,$4)
+
#
#rm -rf *.so
