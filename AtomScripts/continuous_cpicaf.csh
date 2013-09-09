#! /bin/tcsh
###################################
# created this to copy/move files from my icaf
# account when jobs are un on caf and output
# is sent to icaf. Not complete yet.
# What I want is to time to time (every 5 min or so)
# check if any files are there, move them to 
# a new location so icaf will have disk space
# for the rest of the job. Espcially in stn making
# for cp2 calibration.
# For this I cannot just copy any file delete it.
# I got to check if the file is copied completely from
# caf to icaf. For this I was hoping to check the 
# size of the file. I'll loop and stat the file
# untill its size is contant and them move it.
# Or it is better if the time stamp is <1min from
# current time, skip it till next time.
#
#
###################################

echo [`date`] "Starting continous file copy/rm from fcdficaf account"
set a = 1
set node="fcdficaf3.fnal.gov"
set srcdir="/cdf/spool/samantha/"
set destdir="samantha@praseodymium.fnal.gov:/cdf/scratch/samantha/p27stn/"
set destdir2="samantha@promethium.fnal.gov:/cdf/scratch/samantha/p27stn/"
set destdir3="samantha@holmium.fnal.gov:/cdf/scratch/samantha/p27stn/"

set DESTDIR=$destdir
while ( $a > 0 )
		#set file = `ssh $node "ls -1 ${srcdir}"`
		set file = `ssh $node "find ${srcdir}* -amin +5"`
		#echo $file
		#echo "size = ${#file[*]}"
		foreach f ( $file )
			#echo "for loop begin"
			echo "copying ${f} to ${DESTDIR}"
			#ssh ${node} "ls -l ${srcdir}| grep ${file} | cut -d ' ' -f5"
			#set oldsize = `ssh ${node} "ls -l ${srcdir}| grep ${file} | cut -d ' ' -f5"`
			#set newsize = 0
			#set ready = 1
			#echo "${f} size old/new= ${oldsize} ${newsize}"
			#while ( $ready > 0 )
			#	echo "while loop : sleeping 2"
			#	sleep 2
			#	@ newsize = `ssh ${node} "ls -l ${srcdir} | grep ${file} | cut -d ' ' -f5"`
			#	@ ready = $newsize - $oldsize
			#	echo "${f} size old/new= ${oldsize} ${newsize}"
			#	echo "ready = $ready"
			#	@ oldsize = $newsize
			#end
			#echo "${node}:${f} ${destdir}"

			#try 3 locations for copying
			set stat1 = `ssh ${node} "scp ${node}:${f} ${destdir}"`
			if ( $stat1 > 0 ) then
				echo "destination 1 : $destdir failed!. trying detsination 2"
				set stat2 = `ssh ${node} "scp ${node}:${f} ${destdir2}"`
				@ stat1 = $stat2

				if ( $stat2 > 0 ) then
					echo "destination 2 : $destdir failed!. trying detsination 3"
					set stat3 = `ssh ${node} "scp ${node}:${f} ${destdir3}"`
					@ stat1 = $stat3
				end
			end

			if ( $stat1 == 0 ) then
				echo "Removing file ${f} with cp status = ${stat}"
				`ssh ${node} "rm -rf ${f}"` 
			end

		end
		#`ssh rsync -av samantha@${node}:${srcdir} praseodymium.fnal.gov:/cdf/scratch/samantha/p27stn/`
		#echo "sleeping ...600s"
		sleep 120
		date
end
echo "Continuous probing is comlete."
