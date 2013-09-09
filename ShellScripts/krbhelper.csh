#! /bin/tcsh
echo [`date`] Starting krbhelper
set KRB5CCFILE=`echo ${KRB5CCNAME} | /usr/bin/awk -F: '{if($1=="FILE") print $2}'`
echo "krbhelper copy "${KRB5CCFILE} /tmp/krb5${USER}
/bin/cp ${KRB5CCFILE} /tmp/krb5${USER}
echo "krbhelper klist"
klist
# assume the arguments are the command to be run in the background
echo
echo ">>>>>>> krb helper: running command: " $* \&
echo
$* &
