#! /bin/tcsh
echo [`date`] "Starting continous monitor of ATOM jobs (for 100 times, once 5 minute)"

foreach job ( 1 2 3 4 5 6 7 8 9 10 )
	foreach time ( 1 2 3 4 5 6 7 8 9 10 )
		./probe.csh
		echo "sleeping ..."
		sleep 600
		date
	end
end
echo "Continuous probing is comlete."
