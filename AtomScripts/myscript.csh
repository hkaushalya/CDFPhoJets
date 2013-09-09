#! /bin/tcsh
#I got this from Ray to submit and run jobs on Atom cluster remotely.
echo [`date`] "Starting myscript"
#echo "jobno" $1
#echo "dataset" $2
#echo "subset" $3
#echo "newdir" $4
# pick up the parent process ticket
setenv KRB5CCNAME /tmp/krb5${USER}
echo "myscript klist"
klist

# start actual work
#source ~cdfsoft/cdf2.cshrc
#setup cdfsoft2 6.1.4

#cd $4
cd $2
echo "Got file=$1, dir=$2"
pwd
ls
#which root.exe
# root commands here 
#root -b -q 'cafFlatNtupleMaker.C+ ('$1','$2','$3')'
#root -b -q 'cafCosmicSkim.C+ ('$1','$2','$3')'
echo "myscript got file $1"
$2/$1
echo sleeping
sleep 5

echo [`date`] done

exit
