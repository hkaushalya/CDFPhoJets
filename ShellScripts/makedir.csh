#! /bin/tcsh
echo [`date`] "Starting makedir script"
#cd /cdf/scratch/samantha/
#rm -rf $1
#mkdir $1
cd $1
#rm -rf *
#scp ~/.rootrc .
tar xzf ~/mycafjob.tgz
ls -l
echo sleeping
sleep 2
echo makedir [`date`] done
exit
