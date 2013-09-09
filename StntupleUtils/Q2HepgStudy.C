#include "TBenchmark.h"
#include "TSystem.h"
#include <iostream>
#include <sstream>
#include <string>
#include "TObjArray.h"
#include "TIterator.h"
#include "Stntuple/Stntuple/loop/TStnAna.hh"
#include "Stntuple/Stntuple/loop/TStnInputModule.hh"
#include "Stntuple/Stntuple/loop/TStnOutputModule.hh"
#include "Stntuple/Stntuple/loop/TStnCatalog.hh"
#include "Stntuple/base/base/TStnDataset.hh"
#include "samantha/samantha/Pho/PhoJetsHepgSyst.hh"

void Q2HepgStudy(int dataset, int subset, int debug = 0) 
{
	assert (dataset>=-1 && dataset<=1 && "Q2HepgStudy::Invalid dataset given! exiting.");
	std::cout << "I got dataset, subset =" << dataset << ", " <<  subset << std::endl;
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");


	TChain *ch = new TChain("STNTUPLE");
	std::string filename;

	if (dataset == 0)  //base
	{
		if (subset == 20)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_1");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_2");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_3");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_4");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_5");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_6");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_7");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_base.root_8");

			filename = "Q2BasePtHat20Results.root";
		} else if (subset == 60)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.245316_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/base/stn.246147_4.root");

			filename = "Q2BasePtHat60Results.root";

		} else if (subset == 100)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.246166_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/base/stn.245546_20.root");

			filename = "Q2BasePtHat100Results.root";

		} else if (subset ==160)
		{

			/*
			 * these files does not work either. location is on nbay03
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_10.root");
			
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.245559_8.root");
			*/
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246667_11.root");
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246670_13.root");
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246671_16.root");
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246672_18.root");
//			ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246672_19.root");
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246662_3.root");
			//ch->Add("root://nbay02.fnal.gov:5151//mnt/autofs/misc/nbay03.a/samantha/Take9_PtHat160/base/stn.246662_9.root");

 			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_1.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_10.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_12.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_14.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_15.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_17.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_2.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_20.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_4.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_5.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_6.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_7.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.245559_8.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246667_11.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246670_13.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246671_16.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246672_18.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246672_19.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246662_3.root");
			ch->Add("/data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base/stn.246663_9.root");

			
/*
 * for some reason root does not like these files!
 			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.245559_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246667_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246670_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246671_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246672_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246672_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246662_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/base2/stn.246662_9.root");
*/
			
			filename = "Q2BasePtHat160Results.root";
		} else if (subset == 210)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.245572_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.246690_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.246691_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.246692_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/base/stn.246692_19.root");

			filename = "Q2BasePtHat210Results.root";
		}

	} else if (dataset==1)  //more Q2 MC sample
	{
		std::cout << "Using Q2 Up sample." << std::endl;
		/*ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set1of3.root");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set1of3.root_1");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set1of3.root_2");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set1of3.root_3");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set1of3.root_4");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set2of3.root");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set3of3.root");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set3of3.root_1");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set3of3.root_2");
		  ch->Add("root://praseodymium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyMoreQ2Set3of3.root_3");
		  */

		if (subset == 20)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_1");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_10");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_11");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_12");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_13");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_14");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_2");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_3");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_4");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_5");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_6");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_7");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_8");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2up.root_9");
			filename = "Q2UpPtHat20Results.root";

		} else if (subset == 60)
		{
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_1.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_10.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_11.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_12.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_13.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_14.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_15.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_16.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_17.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_18.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_19.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_2.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_20.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_3.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_4.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_5.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_6.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_7.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_8.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2up/stn.245342_9.root");

			filename = "Q2UpPtHat60Results.root";


		} else if (subset ==100)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.245554_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.246626_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.246627_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2up/stn.246658_15.root");

			filename = "Q2UpPtHat100Results.root";

		} else if (subset == 160)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.245567_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246673_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246673_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246673_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246676_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246676_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246676_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246676_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246676_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2up/stn.246679_3.root");

			filename = "Q2UpPtHat160Results.root";
		} else if (subset == 210)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2up/stn.246055_9.root");

			filename = "Q2UpPtHat210Results.root";
		}

	} else if (dataset == -1) //less Q2
	{
		std::cout << "Using Q2 Down sample." << std::endl;
		/*ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set1of3.root");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set1of3.root_1");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set1of3.root_2");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set1of3.root_3");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set1of3.root_4");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set2of3.root");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set2of3.root_1");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set2of3.root_2");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set3of3.root");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set3of3.root_1");
		  ch->Add("root://thulium.fnal.gov:5151//cdf/scratch/samantha/stnrel/stn.1392PythiaHepgOnlyLessQ2Set3of3.root_2");
		  */
		if (subset == 20)
		{

			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_1");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_10");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_11");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_12");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_13");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_2");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_3");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_4");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_5");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_6");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_7");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_8");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b//samantha/MCSamples/Take5_PtHat20/Stntuples/stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_Q2down.root_9");
			filename = "Q2DownPtHat20Results.root";
		} else if (subset == 60)
		{
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_1.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_10.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_11.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_12.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_13.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_14.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_15.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_16.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_17.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_18.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_19.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_2.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_20.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_3.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_4.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_5.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_6.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_7.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_8.root");
			ch->Add("~/MCSamples/Take6_PtHat60/Q2down/stn.245344_9.root");


			filename = "Q2DownPtHat60Results.root";

		} else if (subset ==100)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.245556_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/Q2down/stn.246659_18.root");

			filename = "Q2DownPtHat100Results.root";

		} else if (subset == 160)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.245568_9.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246680_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246682_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246682_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246682_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246684_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246685_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246685_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246685_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246686_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/Q2down/stn.246686_19.root");

			filename = "Q2DownPtHat160Results.root";
		} else if (subset == 210)
		{
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_1.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_10.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_11.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_12.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_13.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_14.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_15.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_16.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_17.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_18.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_19.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_2.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_20.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_3.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_4.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_5.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_6.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_7.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_8.root");
			ch->Add("root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/Q2down/stn.246056_9.root");

			filename = "Q2DownPtHat210Results.root";
		}

	}

	TStnAna *ap = new TStnAna(ch);
	ap->GetInputModule()->SetPrintLevel(1);   // print file name as they are opened


	TH1::AddDirectory(kFALSE); //do not add these histos to memory 


	/******************************************************/
	// define all module here
  /******************************************************/
  PhoJetsHepgSyst *pho = new PhoJetsHepgSyst();
  	pho->SetJetCone(0); //0(0.4), 1(0.7), 2(1.0)
		
	ap->AddModule(pho,1);
	
	if (debug) ap->Run(debug);
	else ap->Run();
  	
	
	ap->SaveHist(filename.c_str(),2);
	delete pho;
	delete ch;
	std::cout << "/************* JOB SUMMARY ******************************/" <<std::endl;
	std::cout << "::::::::: Created ROOT file      => "<< filename << std::endl;
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}
void Q2HepgStudy() 
{
	Q2HepgStudy(-1, 60,0); 
	Q2HepgStudy(1, 60,0); 
	
	/*for (int i=-1; i<=1; ++i)
	{
		Q2HepgStudy(i, 20,0);
		Q2HepgStudy(i, 60,0); 
		Q2HepgStudy(i, 100,0);
		Q2HepgStudy(i, 160,0);
		Q2HepgStudy(i, 210,0);
	}
	*/
}
