#! /bin/tcsh

#wrote this to extract met preictions
#when 6 jet filter modules with MetSig cuts
# 2,3,4,5,6, and 7 is run in one job.
#Tested this using Number worksheet (elog#935) -12-12-2008.

set tmpfile17_2="tmp17_2.tmp"
set tmpfile17_3="tmp17_3.tmp"
set tmpfile17_4="tmp17_4.tmp"
set tmpfile17_5="tmp17_5.tmp"
set tmpfile17_6="tmp17_6.tmp"
set tmpfile17_7="tmp17_7.tmp"

set tmpfile18_2="tmp18_2.tmp"
set tmpfile18_3="tmp18_3.tmp"
set tmpfile18_4="tmp18_4.tmp"
set tmpfile18_5="tmp18_5.tmp"
set tmpfile18_6="tmp18_6.tmp"
set tmpfile18_7="tmp18_7.tmp"

set tmpfile19_2="tmp19_2.tmp"
set tmpfile19_3="tmp19_3.tmp"
set tmpfile19_4="tmp19_4.tmp"
set tmpfile19_5="tmp19_5.tmp"
set tmpfile19_6="tmp19_6.tmp"
set tmpfile19_7="tmp19_7.tmp"

cd $PWD
#/bin/rm -i tmp*.tmp
/bin/rm -rf tmp*.tmp

set filelist=`ls -1 *.out*`
#set filelist=`ls -1 $1`
foreach file ( $filelist )
	echo $file "::"
	grep -m1 "TMJ:17" $file >> $tmpfile17_2
	grep -m2 "TMJ:17" $file | tail -1 >> $tmpfile17_3 
	grep -m3 "TMJ:17" $file | tail -1 >> $tmpfile17_4 
	grep -m4 "TMJ:17" $file | tail -1 >> $tmpfile17_5 
	grep -m5 "TMJ:17" $file | tail -1 >> $tmpfile17_6 
	grep -m6 "TMJ:17" $file | tail -1 >> $tmpfile17_7 

	grep -m1 "TMJ:18" $file >> $tmpfile18_2
	grep -m2 "TMJ:18" $file | tail -1 >> $tmpfile18_3 
	grep -m3 "TMJ:18" $file | tail -1 >> $tmpfile18_4 
	grep -m4 "TMJ:18" $file | tail -1 >> $tmpfile18_5 
	grep -m5 "TMJ:18" $file | tail -1 >> $tmpfile18_6 
	grep -m6 "TMJ:18" $file | tail -1 >> $tmpfile18_7 

	grep -m1 "TMJ:19" $file >> $tmpfile19_2
	grep -m2 "TMJ:19" $file | tail -1 >> $tmpfile19_3 
	grep -m3 "TMJ:19" $file | tail -1 >> $tmpfile19_4 
	grep -m4 "TMJ:19" $file | tail -1 >> $tmpfile19_5 
	grep -m5 "TMJ:19" $file | tail -1 >> $tmpfile19_6 
	grep -m6 "TMJ:19" $file | tail -1 >> $tmpfile19_7 
end

echo " MET SIG 2 ========"
set dataobserved=`cat $tmpfile17_2 | awk '{f+=$6}END{print f}'`
#echo "Data observed            = $dataobserved" 
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_2 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_2 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_2 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_2 | awk '{f+=$8}END{print f}' 


echo " MET SIG 3 ========"
#echo -n "Data observed            = " 
#cat $tmpfile17_3 | awk '{f+=$6}END{print f}' 
set dataobserved=`cat $tmpfile17_3 | awk '{f+=$6}END{print f}'`
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_3 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_3 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_3 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_3 | awk '{f+=$8}END{print f}' 

echo " MET SIG 4 ========"
#echo -n "Data observed            = " 
#cat $tmpfile17_4 | awk '{f+=$6}END{print f}' 
set dataobserved=`cat $tmpfile17_4 | awk '{f+=$6}END{print f}'`
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_4 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_4 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_4 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_4 | awk '{f+=$8}END{print f}' 

echo " MET SIG 5 ========"
#echo -n "Data observed            = " 
#cat $tmpfile17_5 | awk '{f+=$6}END{print f}' 
set dataobserved=`cat $tmpfile17_5 | awk '{f+=$6}END{print f}'`
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_5 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_5 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_5 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_5 | awk '{f+=$8}END{print f}' 

echo " MET SIG 6 ========"
#echo -n "Data observed            = " 
#cat $tmpfile17_6 | awk '{f+=$6}END{print f}' 
set dataobserved=`cat $tmpfile17_6 | awk '{f+=$6}END{print f}'`
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_6 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_6 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_6 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_6 | awk '{f+=$8}END{print f}' 

echo " MET SIG 7 ========"
#echo -n "Data observed            = " 
#cat $tmpfile17_7 | awk '{f+=$6}END{print f}' 
set dataobserved=`cat $tmpfile17_7 | awk '{f+=$6}END{print f}'`
#echo -n "Prediction+/-syst+/-stat = " 
set prediction=`cat $tmpfile18_7 | awk '{f+=$6}END{print f}'`
set syst=`cat $tmpfile18_7 | awk '{f+=($8*$8)}END{print sqrt(f)}'`
set stat=`cat $tmpfile18_7 | awk '{f+=($10*$10)}END{print sqrt(f)}'`
#echo "$prediction +/- $syst +/- $stat"
echo "<td>$dataobserved</td><td>$prediction +/- $syst +/- $stat</td>"
#echo -n "Events Before cuts       = " 
#cat $tmpfile19_7 | awk '{f+=$8}END{print f}' 

