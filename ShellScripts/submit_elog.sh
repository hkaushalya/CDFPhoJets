#!/bin/sh
#############
#this is script to submit stuff to elog from terminal.
#this take only few seconds as compared to 30-40 min 
# upload ~20 images to elog.
# remove option '-e' to create a new entry
# NOTE: BE VERY CAREFUL with option '-e'.
#you can replace an old entry with this option.
# this mthd ignore the setting NO EDIT to entries after 48rs.
#Type the password in xxxxx
############

elog -h hep05.baylor.edu -d hep -l samantha -v \
-a Type="Information" \
-a Category="Gamma+Jets" \
-a Author="Samantha Hewamanage" \
-a Subject="P1-26 Reweighted sideband: Errors for Method A" \
-x \
-m text.txt \
-f $PWD/plot1_Et_pho.pdf \
-f $PWD/plot1_Et_j1.pdf \
-f $PWD/plot2_Met.pdf \
-f $PWD/plot1_DetEta_j1.pdf \
-n 1 \
-u samantha xxxxx \
-e 1646
