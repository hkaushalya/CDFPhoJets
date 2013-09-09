#!/bin/sh

if [ $1 -eq 11 ]
then

grep "datasets have" job_1.out
echo -n "Total Processed: " 
grep "TRM:01" *.out | awk '{f+=$6}END{print f}'
echo -n "Total Saved: " 
grep "P2J:02" *.out | awk '{f+=$6}END{print f}'

fi

##########################################

if [ $1 -eq 2 ]
then
echo -n "tgz files = "
ls *.tgz |wc -l
echo -n "out files = "
ls *.out |wc -l
echo " extracting *.out files"
`~/samantha/RootTools/extractfiles \*.out`
grep "datasets have" job_1.out
echo -n "BeginJob: chained: "
grep "BeginJob: chained" *.out | awk '{f+=$5}END{print f}'
echo -n "Total Processed: "
grep "TRM:01" *.out | awk '{f+=$6}END{print f}'
echo -n "Total Saved: "
grep "P2J:02" *.out | awk '{f+=$6}END{print f}'
echo " extracting *.err files"
`~/samantha/RootTools/extractfiles \*.err`
echo " extracting ROOT files"
`~/samantha/RootTools/extractfiles \*.root`

fi

#########################################

if [ $1 -eq 3 ]
then

for file in `ls -1tr | grep out$`  #preserve time order
do
echo " ===================== $file "
grep "datasets have" $file
echo -n "BeginJob: chained: "
grep "BeginJob: chained" $file
echo -n "Total Processed: "
grep "TRM:01" $file
echo ""
#grep "TRM:01" $file | awk '{f+=$6}END{print f}'
done
fi
