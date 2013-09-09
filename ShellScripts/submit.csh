#! /bin/tcsh
set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium Thulium Ytterbium"
#set nodes="nbay01 nbay02 nbay03 nbay04 nbay05"
#set nodes="Bromine chlorine fluorine oxygen sulfur nitrogen Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"

#scp ~/stn_rel/.rootrc lanthanum:~/
#scp ~/oldhome/caf/mycafjob.tgz holmium:~/
#scp -v /mnt/autofs/misc/nbay02.b/samantha/caf/mycafjob.tgz holmium:~/
set dataset=10
set subset=5              #last jobs i started, change this for new one
set cafsections=200

set newdir="/cdf/scratch/samantha/stuple"
#foreach node ( $nodes )
#	echo "making dir $newdir in $node"
#   ssh $node "~/makedir.csh $newdir >> stuple.mkdirlog "
#end

# index of jobs
#set jobno=0
#  foreach node ( $nodes )
#	 echo $jobno
#    set id = `printf "%02i_%02i_%02i" $dataset $subset $jobno`
#    ssh $node "~/krbhelper.csh ~/myscript.csh $jobno $dataset $subset $newdir >&! ${newdir}/dataset_${id}.log "
#    @ jobno = $jobno + 1
#	 sleep 5
#  end
#
#
#
# index of jobs
#exit
set jobno=176
# 2 copies of the job per node
foreach j (1 2)
  set i=0
  foreach node ( $nodes )
    # jobno_nodeno_copy
   # set id = `printf "%02i_%02i_%02i" $jobno $i $j`
#    ssh $node "~/krbhelper.csh ~/myscript.csh $jobno $dataset $subset $newdir >&! ${newdir}/dataset${dataset}_${subset}_${id}.log "
    set id = `printf "%02i_%02i_%02i" $dataset $subset $jobno`
	 echo "${jobno}::${node}::${id}"
    ssh $node "~/krbhelper.csh ~/myscript.csh $jobno $dataset $subset $cafsections $newdir >&! ${newdir}/dataset_${id}.log "
    @ i = $i + 1
    @ jobno = $jobno + 1
	 sleep 2
  end
end

#now renice all the spawned jobs
foreach node ( $nodes )
	set list=`ssh $node "ps ax -u samantha | grep "root.exe" | grep -v "grep" | cut -c 1-5"`
	foreach ps ( $list ) 
		echo "renice job ${ps}"
		ssh $node "renice 19 ${ps}"
	end
end

