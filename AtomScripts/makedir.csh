#! /bin/tcsh
echo [`date`] "Starting makedir script"
cd /cdf/scratch/samantha/
rm -rf $1
pwd
mkdir $1
cd $1
pwd
rm -rf *
#scp ~/.rootrc .
#tar xzf ~/mycafjob.tgz
tar xzf ~/cpr.tgz
ls -l
echo sleeping
sleep 2
echo makedir [`date`] done
exit
