#! /bin/tcsh

set nodes="Lanthanum Praseodymium Promethium Samarium Europium Gadolinium Terbium Holmium Erbium  Thulium Ytterbium Actinium"
dir="/cdf/scratch/samantha/zee/"
file="dataset*"

foreach node ( $nodes )
	echo $node "::"
   ssh $node "~/checkandlog $dir $file"
	sleep 1
end
