#! /bin/tcsh
echo [`date`] "Starting myscript_temp"
echo "Nevts   = $1"
echo "pdfset  = $2"
echo "cafsection   = $3"
echo "Ncafsections = $4"
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
./samantha/pdftools/run_pdf_syst_caf.csh $1 $2 $3 $4
echo sleeping
sleep 20

echo [`date`] done

exit
