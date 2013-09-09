#! /bin/tcsh
#set nodes="Lanthanum Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium Actinium"
set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"
#set nodes="Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"
#set nodes="Samarium Europium"
#set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium"
#set nodes="Erbium"

#scp ~/oldhome/caf/mycafjob.tgz lanthanum:~/
#scp ~/cp2stnrel/cpr.tgz ytterbium:~/

#scp ~/samantha/AtomScripts/makedir.csh ytterbium:~/
#scp ~/samantha/AtomScripts/krbhelper.csh ytterbium:~/
#scp ~/samantha/AtomScripts/myscript.csh ytterbium:~/
set newdir="/cdf/scratch/samantha/p27stn"

#MAKE NEW DIRECTORY FOR THIS CALIBRATION PERIOD
#foreach node ( $nodes )
#	echo "makedir:::$node"
# ssh $node "~/makedir.csh $newdir >> & ! cp2p27.mkdirlog "
#end

set firstrun=284835
set lastrun=287261
set run1=$firstrun
set run2=0
set stream="a"
#set stream="b"
set step=250
#echo "run1=$run1 run2=$run2"
#index of jobs

#foreach job ( 1 2 )
foreach job ( 1 )

foreach node ( $nodes )
	echo ":::::: job=$job :: $node"
#	echo ":::::: $node"

	@ run2 = $run1 + $step - 1
	if ($run2 > $lastrun) then
		@ run2 = $lastrun
	endif

	echo "run1-run2=$run1 to $run2"
	set runFile = "runCaf_${stream}${run1}_${run2}.sh"
	set tclFile = "P27_${stream}${run1}_${run2}.tcl"
	set logFile = "${stream}${run1}_${run2}.log"
	echo "runFile=$runFile tclFile=${tclFile} logFile=${logFile}"
	
	scp ~/cp2stnrel/cpr2/$runFile ${node}:${newdir}
	scp ~/cp2stnrel/cpr2/cprtcl/$tclFile ${node}:${newdir}/cprtcl
	ssh $node "~/krbhelper.csh ~/myscript.csh $runFile $newdir >&! ${newdir}/${logFile}"

	
	#ssh $node "~/krbhelper.csh ~/myscript.csh $runFile $newdir >&! ${newdir}/a270419_270500.log"
	#ssh $node "~/krbhelper.csh ~/startrootd.csh"
	#ssh $node "~/krbhelper.csh ~/startrootd.csh >&! ~/${node}_rootd.log"
	
	@ run1 = $run1 + $step
	if ( $run1 > $lastrun ) then
		break
	endif

	#wait sometime for the first job to setup sam stuff
	if ( ($job == 1) && ($node == "Praseodymium") ) then
		echo "sleeping 500 $job $node"
		sleep 500
	endif
	echo "sleeping 5 (default)"
	sleep 5

	end

end
