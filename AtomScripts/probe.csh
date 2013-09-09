#! /bin/tcsh

#set nodes="Lanthanum Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium Actinium"
set nodes="Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"
#set nodes="Praseodymium Promethium" 
#set nodes="Praseodymium" 
#set nodes="Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium"

set run=266528
set step=60


foreach node ( $nodes )
	echo ">>>>>>>>>> ${node} <<<<<<<<<<<<<"
	#ssh $node "/cdf/atom/home/samantha/makedir.csh"
	#foreach j (0 1)
	#echo "${node}::$j"
	
	#echo $run
	#scp ${HOME}/cp2stnrel/cpr2/*${run}* $node.fnal.gov:~/scratch/p20stnhv/
	#scp ${HOME}/cp2stnrel/cpr2/cprtcl/*${run}* $node.fnal.gov:~/scratch/p20stnhv/cprtcl/
	#set max=0
	#@ max = $run + $step - 1
	#echo $max
	#set filea="runCaf_a${run}_${max}.sh"
	#echo $filea
	#ssh $node "pwd"
   #ssh $node "~/krbhelper.csh /cdf/scratch/samantha/p20stnhv/${filea} >&! /cdf/scratch/samantha/p20stnhv/a${run}_${max}.log"
	#@ run = $run + $step
	
	
	#if [ $1 -eq 1 ]
	#then 
    	#ssh $node "ps -al r -U $USER | grep root"
   	#ssh $node "ps -al -U $USER | grep root"

	#elif [ $1 -eq 2 ]
	#then
	#	set list=`ssh $node "ls /cdf/scratch/samantha/zee/dataset*"`
	# 	foreach file ( $list ) 
	# 		echo $file
   #		ssh $node "~/summary $file"
	# 	end

#	else 
		
		#set list=`ssh $node "ls /cdf/scratch/samantha/phodata/dataset*"`
		#set list=`ssh $node "ls /cdf/scratch/samantha/cosmicskim/dataset*"`
		#set list=`ssh $node "ls /cdf/scratch/samantha/zee/dataset*"`
		#set list=`ssh $node "ls -1 /cdf/scratch/samantha/cp2stn4/*.log"`
		#set list=`ssh $node "ls -1 /cdf/scratch/samantha/p20stn/*.log"`
		#set list=`ssh $node "ls -1 /cdf/atom/home/samantha/p21stn/*.log"`
		#set list=`ssh $node "ls -1 /cdf/scratch/samantha/p23stn/*.log"`
		set list=`ssh $node "ls -1 /cdf/scratch/samantha/p25stn/*.log"`
		ssh $node "ls -ltr /cdf/scratch/samantha/p25stn/*.log"
		echo ">>> List of stntuples <<<<"
		#ssh $node "ls -1 /cdf/scratch/samantha/p24stn/*.stn*"
		echo -n "Stnfiles = "
		ssh $node "ls -1 /cdf/scratch/samantha/p24stn/*.stn* | wc -l"
		#set list=`ssh $node "top -b -n 1"`
	 	foreach file ( $list ) 
			#echo $file
	 		#grep "stnmaker" $file
   		#ssh $node "~/processed $file"
   		#ssh $node "~/stnf_processed $file"
   		#ssh $node "cat $file | grep Error"
   		ssh $node "tail -10 $file"
   		#ssh $node "cat $file | grep  'Your job'"
	 	end
		#ssh $node "ls -1 /cdf/scratch/samantha/p20stn/*.stn*"
		

	#fi
	#end

	#sleep 1
end
