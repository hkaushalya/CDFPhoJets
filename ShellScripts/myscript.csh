#! /bin/tcsh
echo [`date`] "Starting myscript"
echo "jobno" $1
echo "dataset" $2
echo "subset" $3
echo "cafsections $4"
echo "newdir" $5
# pick up the parent process ticket
setenv KRB5CCNAME /tmp/krb5${USER}
echo "myscript klist"
klist

# start actual work
source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 6.1.4

cd $5
which root.exe
# root commands here 
root -b -q 'cafFlatNtupleMaker.C+('$1','$2','$3','$4')'
#root -b -q 'cafCosmicSkim.C+ ('$1','$2','$3')'

echo sleeping
sleep 20

echo [`date`] done

exit
