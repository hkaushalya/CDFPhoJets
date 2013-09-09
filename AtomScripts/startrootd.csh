#! /bin/tcsh
#Start a rootd server in ATOM node.
echo [`date`] "Starting start rootd ..."

# pick up the parent process ticket
setenv KRB5CCNAME /tmp/krb5${USER}
klist

# start actual work
source ~cdfsoft/cdf2.cshrc
setup cdfsoft2 6.1.4
rootd -p 5151
#sleep 2
ps aux | grep rootd
echo "rootd started."
exit
