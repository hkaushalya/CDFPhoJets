#! /bin/tcsh

set tmpfile17_1="tmp17_1.tmp"
set tmpfile17_2="tmp17_2.tmp"
set tmpfile17_3="tmp17_3.tmp"
set tmpfile18_1="tmp18_1.tmp"
set tmpfile18_2="tmp18_2.tmp"
set tmpfile18_3="tmp18_3.tmp"
set tmpfile20_1="tmp20_1.tmp"
set tmpfile20_2="tmp20_2.tmp"
set tmpfile20_3="tmp20_3.tmp"

rm -i tmp*.tmp

#set filelist=`ls -1 *.out`
set filelist=`ls -1 $1`
foreach file ( $filelist )
	echo $file "::"
	grep -m1 "TMJ:17" $file >> $tmpfile17_1
	grep -m2 "TMJ:17" $file | tail -1 >> $tmpfile17_2 
	grep -m3 "TMJ:17" $file | tail -1 >> $tmpfile17_3 
	grep -m1 "TMJ:18" $file >> $tmpfile18_1
	grep -m2 "TMJ:18" $file | tail -1 >> $tmpfile18_2
	grep -m3 "TMJ:18" $file | tail -1 >> $tmpfile18_3
	grep -m1 "TMJ:20" $file >> $tmpfile20_1
	grep -m2 "TMJ:20" $file | tail -1 >> $tmpfile20_2
	grep -m3 "TMJ:20" $file | tail -1 >> $tmpfile20_3
end

echo " MET SIG 4 ========"
#echo -n "Events after cuts        = " 
#cat $tmpfile20_1 | awk '{f+=$8}END{print f}' 
echo -n "Data observed            = " 
cat $tmpfile17_1 | awk '{f+=$6}END{print f}' 
echo -n "Prediction+/-syst+/-stat = " 
cat $tmpfile18_1 | awk '{f+=$6}END{print f}' 
cat $tmpfile18_1 | awk '{f+=($8*$8)}END{print sqrt(f)}' 
cat $tmpfile18_1 | awk '{f+=($10*$10)}END{print sqrt(f)}' 

echo " MET SIG 5 ========"
echo -n "Events after cuts        = " 
cat $tmpfile20_2 | awk '{f+=$8}END{print f}' 
echo -n "Data observed            = " 
cat $tmpfile17_2 | awk '{f+=$6}END{print f}' 
echo -n "Prediction+/-syst+/-stat = " 
cat $tmpfile18_2 | awk '{f+=$6}END{print f}' 
cat $tmpfile18_2 | awk '{f+=($8*$8)}END{print sqrt(f)}' 
cat $tmpfile18_2 | awk '{f+=($10*$10)}END{print sqrt(f)}' 

echo " MET SIG 6 ========"
echo -n "Events after cuts        = " 
cat $tmpfile20_3 | awk '{f+=$8}END{print f}' 
echo -n "Data observed            = " 
cat $tmpfile17_3 | awk '{f+=$6}END{print f}' 
echo -n "Prediction+/-syst+/-stat = " 
cat $tmpfile18_3 | awk '{f+=$6}END{print f}' 
cat $tmpfile18_3 | awk '{f+=($8*$8)}END{print sqrt(f)}' 
cat $tmpfile18_3 | awk '{f+=($10*$10)}END{print sqrt(f)}' 
