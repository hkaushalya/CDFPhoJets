#!/home/cdfsoft/products/perl/v5_005/Linux+2/bin/perl
#
#=========================================================================
# New example script for gegtuple in 4.11.2
#
# NB now reconstruction modules are by default disabled.
# Also beamline access is updated so you should use CalibrationManager
# to specify a used set, not MyModule talk-to to specify a historycode.
#
#    1.  Just to build the ROOT library:
#            addpkg -h ElectronUser
#    2.  To run the executable:
#        Fix memory leak, update Phoenix for CAF, get right plug calibs
#            addpkg -h ElectronUser
#            addpkg -h PADSMods
#            addpkg -h PADSObjects
#            addpkg PhoenixMods
#            cvs update -r 1.14 PhoenixMods/src/PhoenixTracking.cc
#            cvs update -r 1.8 PhoenixMods/PhoenixMods/PhoenixTracking.hh
#            addpkg CalorObjects
#            cvs co -r 1.79 CalorObjects/src/Calib.cc
#========================================================================
#run this script gives me the tcl files (1 for each job) and
# one shell script runCaf_.sh;
# setup cdfsoft development 
# CafGui
#_


#=======================================================================
# *********************** CVS Commit Log *******************************
# $Log: loop_physr.pl,v $
# Revision 1.4  2009/04/05 20:27:18  samantha
# Enabled CALQ module and makeMet to get SumEt,MEt info to find a better criteria for the cp2 HV Off. Commented out the CalibrationManager which load L3 calibs as I do not need them. There are few stuff commented out that have not fully implemented yet.
#
#=======================================================================

#========================================================================
# Write a file called myfilesetlist.lst, with one fileset per line
# (or file, and modify $infile and DHInput talk-to accordingly
#========================================================================
#$listname = "aug03";
#$runnum=246125;
$listname = "P23";
$runnum_min = 272470;
$runnum_max = 274055;
$nodes = 10;		#specify the number of nodes(computers) you'll run on
$jobspernode = 2; #specify the number of jobs per node

#$step = ($runnum_max - $runnum_min)/($nodes * $jobspernode);
#print $step;
$step = 80;       #specify the runs per job

#$datalist = "$listname.lst";
# Or you can put your fileset lists in a directory
#$datalist = "dataStreams/$listname.lst";
#$datalist = "dataStreams/mytest_all.lst";
$datalist = "dataStreams/mytest.lst";
open(LIST,"<$datalist");
chomp(@filesetlist=<LIST>);
close(LIST);
print "going\n";
$loop = 1;
while ($loop ne 0) {
$min = $runnum_min;
$runnum_min += $step;  #need an offset of 1 to remove overlap, which done at the end of the loop
$max = $runnum_min - 1;
if ($max >= $runnum_max) {
	$max = $runnum_max;
	$loop = 0;
}
print "min,max = $min   $max  loop=$loop\n";
foreach $fileseti (@filesetlist){
$fileset="${fileseti}${min}_${max}";
$cafsh = "runCaf_$fileset.sh";
print "$fileset\n";

#========================================================================
# Make a directory for tcl, here called tcltemp
#========================================================================
$tclfile="cprtcl/$listname\_$fileset.tcl";
$infile = $fileset;
# or these can all be something like
#$infile = root://mycomputer//mylocation/$fileset*;
$outfile = "$listname\_$fileset.stn";
$logfile = "$listname\_$fileset.log";
$lumlog  = "$listname\_$fileset.lumlog";
$stream  = "${fileseti}physr";
open OUT, ">$tclfile" or die "Cannot open $out for write";
print OUT <<EOF;

talk DHInput 
 dropList set CosmicRayInfo CdfMuonView::CosmicMuons
 dropList add ZVertexColl::ZVertexColl
 requireCatalog set F
 
 include dataset ${stream} run>=$min run<=$max
 #selectEvents set run=$runnum
 report set 5000
 show include
exit
#
#loads L3 calibrations- see elog#1088
#talk CalibrationManager
#  ProcessName set L3_PHYSICS_CDF
#exit
module enable  CalqModule
#we do not use this. and to overcome the need for a new Calib.cc file for very recent runs, disabled this. see elog https://hep.baylor.edu/elog/samantha/532 - sam,05-07-2008
module enable  CalorimetryModule
#module disable  CalorimetryModule
module enable InitStntuple StntupleMaker FillStntuple
path disable AllPath
path create  StntuplePath    ManagerSequence

# see above comment on CalorimetryModule
#path append  StntuplePath    CalqModule
path append  StntuplePath    CalqModule        \\
                             CalorimetryModule
path append  StntuplePath    SimInitManager            \\
                             Prereq                    \\
                             InitStntuple              \\
                             StntupleFilter            \\
                             StntupleMaker             \\
                             FillStntuple
path enable StntuplePath
module talk StntupleMaker
 histfile $outfile
 processName   set PROD
 jetCollName set JetCluModule-cone0.4 JetCluModule-cone0.7 JetCluModule-cone1.0
 splitMode set -1
#                                   lines below can be uncommented if necessary
  makeL3Summary set 1
  makeCalData   set 1
# Ray has set this to 0 to fix the need for new Calib.cc file for very recent runs. see above comment. -sam -05-07-2008
  makeCprData   set 0 
  makeCp2Data   set 1
  makeCcrData   set 1
  makeCesData   set 1
  makePesData   set 0
  makePesCorrectedData set 0
  makeClcData   set 0
  makeCmuData   set 0
  makeCmpData   set 0
  makeCmxData   set 0
  makeElectrons set 0
  makeDcasData  set 0
  makeJets      set 0
  makeMet       set 1
  makeMuons     set 0
  makeCosmic    set 0
  makePhotons   set 0
  makeClusters  set 0
  makeTaus      set 0
  makeTracks    set 1
  makeTrigger   set 0
  makeVertices  set 1
  makeZVertices set 1
  makeGenp      set 0
  makeObsp      set 0
#
# these are off:
#
  makeCotData       set 0
  makeTopTElectrons set 0
  makeTopLElectrons set 0
  makeFwdDetData    set 0
  makeTopTJets      set 0
  makeTopLJets      set 0
  makeTopMet        set 0
  makeL3Muons       set 0
  makeTopTMuons     set 0
  makeTopLMuons     set 0
  makeSvxData       set 0
  makeSiStrips      set 0
  #makeTags          set 0
  makeTrackLinks    set 0
  makeTrigSim       set 0
  makeTopSummary    set 0
  makeXft           set 0
  storeXftHits      set 0
  storeXftPixels    set 0
  storeXftTrack     set 0
  makeSvt           set 0
  makeSiGeantIsect  set 0
  makeSiIsect       set 0
  makeSecVtxTag     set 0
  makeJetProb       set 0
  makePi0s          set 0

exit
module disable ConfigManager
action off ProcAction
#begin -nev 1000
begin
show
exit

EOF
    $cafoutfile="${listname}_${fileset}";

open OUT, ">$cafsh" or die "Cannot open $out for write";
print OUT <<EOF;
echo 'path = ' $PATH
export PATH=/usr/krb5/bin:/usr/bin:/bin:${PATH}

export inum=\$1
export caftclfile="cprtcl/${listname}\_${fileset}.tcl"

source ~cdfsoft/cdf2.shrc
./setups.sh
setup cdfsoft2 6.1.4

echo 'LD_LIBRARY_PATH = ' \$\{LD_LIBRARY_PATH\}
export LD_LIBRARY_PATH=.:./shlib:\${LD_LIBRARY_PATH}
export USE_SAM_METADATA=1
#
./stnmaker_prod.exe \$caftclfile

#
X=\$?

#for myfile in `ls *.stn*`
#	do
	#scp \$myfile samantha\@nbay02.fnal.gov:/mnt/autofs/misc/nbay03.a/samantha/cp2stn_jan30/
#done
#	fcp -c stnmaker_prod.exe cdfprd_cal\@fcdfdata401:/cdf/local/disk01/cdfprd_cal/cpr2/fcp/
#	/usr/krb5/bin/rcp -N stnmaker_prod.exe cdfprd_cal\@fcdfdata401:/cdf/local/disk01/cdfprd_cal/cpr2/rcp/
#
#/bin/rm -rf stnmaker_prod.exe
#/bin/rm -rf *.stn*
exit \$X

EOF

system("chmod a+x $cafsh");

#For LSF; don't try to use dcache like this any more (Nov 03)
#system("bsub -q long -o $logfile bin/IRIX6-KCC_4_0/gegtuple_test $tclfile");
#system("bsub -q short -o $logfile bin/IRIX6-KCC_4_0/gegtuple_test $tclfile");


}
}
