#!/bin/ksh
#===========================================================================#
#===========================================================================#
# Script to query to the database the processed luminosity 
# (offline luminosity in nb-1 still not corrected by the factor 1.019)
# of the runs / runsection ranges in a given good run list (first argument)
# for a given dataset (second argument)
#===========================================================================#
#===========================================================================#
# The Script should be executable
# => chmod u+=rwx SAMLumi.ksh
#===========================================================================#
# Required inputs:
# - First argument  ($1): existing file corresponding to your good run list
#                         (can be copy for the good run list web page)
# - Second argument ($2): identifier of the dataset you are interested in 
#   => For instance bhmu then
#      - bhmu0d will be considered for runs < 190000
#      - bhmu0h will be considered for 190000 < runs < 203800
#      - bhmu0i will be considered for 203800 < runs < 228597
#      - bhmu0j will be considered for 228597 < runs < 252835
#      - bhmu0k will be considered for run > 252835
#      (see lines 169 to 172 of the script if you need to change that)
#===========================================================================#
# Ouput = a file called $1_OffLumi_$2.dat
# => Same than good run list input file ($1) but with a 4th column that 
# gives the processed luminosity (offline luminosity in nb-1, still not 
# corrected by the factor 1.019) of the run / runsection range specify on
# the 3 first columns in the chosen dataset ($2)
# => The last lines give the total processed luminosity
# (offline luminosity in nb-1 still not corrected by the factor 1.019)
#===========================================================================#
# Command line example:
# ./SAMLumi.ksh goodrun_em_mu_si.list bhmu
#===========================================================================#
#===========================================================================#
# DFCQuery provided by:   Randolph J. Herber    herber@fnal.gov
# Script written by:      Regis Lefevre         regis@fnal.gov
# Last modification:      July 4th 2006
#===========================================================================#
#===========================================================================#

goodrunlist=$1
dataset=$2

. ~cdfsoft/cdf2.shrc
setup cdfsoft2 6.1.2

if [[ -e $goodrunlist ]] ; then 
  echo "$goodrunlist file found"
else
  echo "$goodrunlist file not found in working directory" 
  echo "=> exit"
  exit
fi

wc -l $goodrunlist > mylinecountcheck
nline=` awk '{print $1}'  mylinecountcheck`
rm mylinecountcheck
echo "Number of lines to be treated $nline"
echo "=> The script usely take about 1 minute for 100 lines"
echo "=> You shall see the first crontrol print out (Treating line 10) in few seconds"
if [[ $nline -gt 5000 ]] ; then
  echo "Number of lines greater than 5000"
  echo "=> SAMLumi.ksh does not handle more than 5000 lines"
  echo "=> Please revisit the SAMLumi.ksh script (see lines 77 to 161)" 
  echo "=> exit"
  exit
fi

output=$1_OffLumi_$2.dat
if [[ -e $output ]] ; then
  rm $output
fi

runs=` awk '{print $1}' $goodrunlist `
mina=` awk '{print $2}' $goodrunlist `
maxa=` awk '{print $3}' $goodrunlist `

set -A minsa
set -A minsb
set -A minsc
set -A minsd
set -A minse
set -A maxsa
set -A maxsb
set -A maxsc
set -A maxsd
set -A maxse

i=0
for min in $mina ; do
 i=$(($i + 1))
 if [[ $i -le 1000 ]] ; then
  minsa[$i]=$min
 fi
 if [[ $i -gt 1000 ]] && [[ $i -le 2000 ]] ; then
  minsb[$(($i - 1000))]=$min
 fi
 if [[ $i -gt 2000 ]] && [[ $i -le 3000 ]] ; then
  minsc[$(($i - 2000))]=$min
 fi
 if [[ $i -gt 3000 ]] && [[ $i -le 4000 ]] ; then
  minsd[$(($i - 3000))]=$min
 fi
 if [[ $i -gt 4000 ]] && [[ $i -le 5000 ]] ; then
  minse[$(($i - 4000))]=$min
 fi
done

i=0
for max in $maxa ; do
 i=$(($i + 1))
 if [[ $i -le 1000 ]] ; then
  maxsa[$i]=$max
 fi
 if [[ $i -gt 1000 ]] && [[ $i -le 2000 ]] ; then
  maxsb[$(($i - 1000))]=$max
 fi
 if [[ $i -gt 2000 ]] && [[ $i -le 3000 ]] ; then
  maxsc[$(($i - 2000))]=$max
 fi
 if [[ $i -gt 3000 ]] && [[ $i -le 4000 ]] ; then
  maxsd[$(($i - 3000))]=$max
 fi
 if [[ $i -gt 4000 ]] && [[ $i -le 5000 ]] ; then
  maxse[$(($i - 4000))]=$max
 fi
done

i=0
t=0
td=0
th=0
ti=0
tj=0
tk=0

for run in $runs ; do

 i=$(($i + 1))

 j=$(print "$i%10"|bc)
 if [[ $j = 0 ]] ; then 
  echo "Treating line $i (over $nline lines)"
 fi

 if [[ $i -le 1000 ]] ; then
  min=${minsa[$i]}
  max=${maxsa[$i]}
 fi
 if [[ $i -gt 1000 ]] && [[ $i -le 2000 ]] ; then
  min=${minsb[$(($i - 1000))]}
  max=${maxsb[$(($i - 1000))]}
 fi
 if [[ $i -gt 2000 ]] && [[ $i -le 3000 ]] ; then
  min=${minsc[$(($i - 2000))]}
  max=${maxsc[$(($i - 2000))]}
 fi
 if [[ $i -gt 3000 ]] && [[ $i -le 4000 ]] ; then
  min=${minsd[$(($i - 3000))]}
  max=${maxsd[$(($i - 3000))]}
 fi
 if [[ $i -gt 4000 ]] && [[ $i -le 5000 ]] ; then
  min=${minse[$(($i - 4000))]}
  max=${maxse[$(($i - 4000))]}
 fi
 mini=$min
 maxi=$max

 if [[ $max = -1 ]] ; then 
   max=65535
 fi 

 datasethere="$dataset"0d
 if [[ $run -gt 190000 ]] ; then
   datasethere="$dataset"0h
 fi 
 if [[ $run -gt 203800 ]] ; then
   datasethere="$dataset"0i
 fi 
 if [[ $run -gt 228597 ]] ; then
   datasethere="$dataset"0j
 fi 
 if [[ $run -gt 252835 ]] ; then
   datasethere="$dataset"0k
 fi  


 DFCQuery "SELECT SUM(DR.LUM_INTEGRAL_OFFLINE) LOFF FROM DATA_FILES SMDF, DATA_FILES_LUMBLOCKS SMDFLMB, CDF2_RUNSECTIONS DR WHERE SMDF.FILE_ID IN (SELECT FILE_ID FROM DATA_FILES DF JOIN DATA_FILES_PARAM_VALUES DFPV USING(FILE_ID) JOIN PARAM_VALUES PV USING(PARAM_VALUE_ID) JOIN PARAM_TYPES PT USING(PARAM_TYPE_ID) JOIN PARAM_CATEGORIES PC USING(PARAM_CATEGORY_ID) WHERE PC.PARAM_CATEGORY='cdf' AND PT.PARAM_TYPE='dataset' AND PV.PARAM_VALUE='$datasethere') AND DR.RUN_NUMBER=$run AND (DR.SECTION_NUMBER BETWEEN $min AND $max) AND SMDF.FILE_ID IN (SELECT FILE_ID FROM DATA_FILES SMDF JOIN DATA_FILES_RUNS SMDFR USING(FILE_ID) JOIN RUNS SMR USING(RUN_ID) JOIN RUN_TYPES SMRT USING(RUN_TYPE_ID) WHERE SMR.RUN_NUMBER=$run) AND SMDF.FILE_ID=SMDFLMB.FILE_ID AND DR.ID BETWEEN SMDFLMB.LUM_MIN AND SMDFLMB.LUM_MAX AND (SMDF.FILE_CONTENT_STATUS_ID IS NULL OR SMDF.FILE_CONTENT_STATUS_ID=(SELECT FILE_CONTENT_STATUS_ID FROM FILE_CONTENT_STATUSES WHERE FILE_CONTENT_STATUS='good')) AND EXISTS (SELECT 1 FROM DATA_FILE_LOCATIONS SMDFL, DATA_STORAGE_LOCATIONS SMDSL WHERE SMDFL.FILE_ID=SMDF.FILE_ID AND SMDFL.LOCATION_ID=SMDSL.LOCATION_ID AND (SMDSL.LOCATION_TYPE!='tape' OR EXISTS (SELECT 1 FROM VOLUMES SMV WHERE SMV.VOLUME_ID=SMDFL.VOLUME_ID AND SMV.VOLUME_STATUS='online'))) ORDER BY ID" > myDFCQuerytmp

 lumi=` awk '{print $1}' myDFCQuerytmp`

 ls -ltr myDFCQuerytmp > myDFCQuerytmpCheck
 size=` awk '{print $5}' myDFCQuerytmpCheck`
 if [[ $size = 2 ]] ; then 
  lumi=0
 fi
 rm myDFCQuerytmp myDFCQuerytmpCheck

 if [[ $datasethere = "$dataset"0d ]] ; then
   td=$(print "$td+$lumi"|bc)
 fi
 if [[ $datasethere = "$dataset"0h ]] ; then
   th=$(print "$th+$lumi"|bc)
 fi
 if [[ $datasethere = "$dataset"0i ]] ; then
   ti=$(print "$ti+$lumi"|bc)
 fi
 if [[ $datasethere = "$dataset"0j ]] ; then
   tj=$(print "$tj+$lumi"|bc)
 fi
 if [[ $datasethere = "$dataset"0k ]] ; then
   tk=$(print "$tk+$lumi"|bc)
 fi 

 t=$(print "$t+$lumi"|bc)

 printf "%.3f" $lumi > myDFCQuerytmp
 lumi=` awk '{print $1}' myDFCQuerytmp`
 rm myDFCQuerytmp

# echo $run $mini $maxi $lumi 
 echo $run $mini $maxi $lumi >> $output

done

printf "%.0f" $t > myDFCQuerytmp
t=` awk '{print $1}' myDFCQuerytmp`
printf "%.0f" $td > myDFCQuerytmp
td=` awk '{print $1}' myDFCQuerytmp`
printf "%.0f" $th > myDFCQuerytmp
th=` awk '{print $1}' myDFCQuerytmp`
printf "%.0f" $ti > myDFCQuerytmp
ti=` awk '{print $1}' myDFCQuerytmp`
printf "%.0f" $tj > myDFCQuerytmp
tj=` awk '{print $1}' myDFCQuerytmp`
printf "%.0f" $tk > myDFCQuerytmp
tk=` awk '{print $1}' myDFCQuerytmp`
rm myDFCQuerytmp

echo "Total processed luminosity = $t nb-1"
echo "=> $td nb-1 in 0d dataset + $th nb-1 in 0h dataset + $ti nb-1 in 0i dataset + $tj nb-1 in 0j dataset + $tk nb-1 in 0k dataset" 
echo "(offline luminosity still not corrected by the factor 1.019)"

echo "Total processed luminosity = $t nb-1"                                                                                          >> $output
echo "=> $td nb-1 in 0d dataset + $th nb-1 in 0h dataset + $ti nb-1 in 0i dataset + $tj nb-1 in 0j dataset + $tk nb-1 in 0k dataset" >> $output
echo "(offline luminosity still not corrected by the factor 1.019)"                                                                  >> $output

exit
