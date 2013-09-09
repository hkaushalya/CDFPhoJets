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
#include "TAuthenticate.h"
/////////////////////////////////////////////////////////////////////////
// TO study the ISR/FSR systematics using only hepg info. 10-10-2009
/////////////////////////////////////////////////////////////////////////

void IsrFsrHepgStudy(int dataset, const float pt_hat, int debug = 0) 
{
	assert (dataset>=0 && dataset<=4 && "IsrFsrHepgStudy::Invalid dataset given! exiting.");
	std::cout << "I got dataset = " << dataset << std::endl;
	std::cout << "Pt-hat        = " << pt_hat << std::endl;
	
	assert ( (pt_hat == 20. || pt_hat ==60. || pt_hat == 100
				|| pt_hat == 160 || pt_hat == 210) && "no matching dataset for the choosen pt-hat");
	
	gROOT->ProcessLine(".!date");
	gBenchmark->Start("metphoana_time");


	TChain *ch = new TChain("STNTUPLE");
	std::stringstream filename;
	std::vector<std::string> filenames;
	std::string path;
	if (dataset==0)  //base MC sample
	{
		std::cout << "Using BASE sample ";
		std::cout << " with minPt-hat " << pt_hat << " GeV";

		if (pt_hat == 20)
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
		} else if (pt_hat == 60)
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


		} else if (pt_hat == 100)
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


		} else if (pt_hat ==160)
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
			
		} else if (pt_hat == 210)
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

		}

		std::cout << "" << std::endl;

		filename << "IsrFsrPythiaPtHat"<< pt_hat <<"_BaseResults.root";

	} else if (dataset == 1) //more ISR
	{
		std::cout << "Using more ISR sample " ;
		if (pt_hat ==20.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take5_PtHat20/Stntuples/";
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root");
			/*filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_1");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_2");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_3");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_4");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_5");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_6");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_7");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreISR.root_8");
*/
		} else if (pt_hat==60.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/moreISR/Job1/";
			filenames.push_back("stn.245318_1.root");
			filenames.push_back("stn.245318_10.root");
			filenames.push_back("stn.245318_11.root");
			filenames.push_back("stn.245318_12.root");
			filenames.push_back("stn.245318_13.root");
			filenames.push_back("stn.245318_14.root");
			filenames.push_back("stn.245318_15.root");
			filenames.push_back("stn.245318_16.root");
			filenames.push_back("stn.245318_17.root");
			filenames.push_back("stn.245318_18.root");
			filenames.push_back("stn.245318_19.root");
			filenames.push_back("stn.245318_2.root");
			filenames.push_back("stn.245318_20.root");
			filenames.push_back("stn.245318_3.root");
			filenames.push_back("stn.245318_4.root");
			filenames.push_back("stn.245318_5.root");
			filenames.push_back("stn.245318_6.root");
			filenames.push_back("stn.245318_7.root");
			filenames.push_back("stn.245318_8.root");
			filenames.push_back("stn.245318_9.root");

		} else if (pt_hat == 100.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/moreISR/";

			filenames.push_back("stn.245548_1.root");
			filenames.push_back("stn.245548_11.root");
			filenames.push_back("stn.245548_12.root");
			filenames.push_back("stn.245548_13.root");
			filenames.push_back("stn.245548_14.root");
			filenames.push_back("stn.245548_15.root");
			filenames.push_back("stn.245548_17.root");
			filenames.push_back("stn.245548_18.root");
			filenames.push_back("stn.245548_19.root");
			filenames.push_back("stn.245548_2.root");
			filenames.push_back("stn.245548_20.root");
			filenames.push_back("stn.245548_3.root");
			filenames.push_back("stn.245548_4.root");
			filenames.push_back("stn.245548_5.root");
			filenames.push_back("stn.245548_6.root");
			filenames.push_back("stn.245548_7.root");
			filenames.push_back("stn.245548_8.root");
			filenames.push_back("stn.245548_9.root");
			filenames.push_back("stn.246171_10.root");
			filenames.push_back("stn.246172_16.root");

		} else if (pt_hat == 160.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/moreISR/";
			filenames.push_back("stn.245561_11.root");
			filenames.push_back("stn.245561_13.root");
			filenames.push_back("stn.245561_14.root");
			filenames.push_back("stn.245561_16.root");
			filenames.push_back("stn.245561_17.root");
			filenames.push_back("stn.245561_19.root");
			filenames.push_back("stn.245561_2.root");
			filenames.push_back("stn.245561_20.root");
			filenames.push_back("stn.245561_8.root");
			filenames.push_back("stn.245561_9.root");


		} else if (pt_hat == 210.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/moreISR/";
			filenames.push_back("stn.245574_1.root");
			filenames.push_back("stn.245574_12.root");
			filenames.push_back("stn.245574_13.root");
			filenames.push_back("stn.245574_14.root");
			filenames.push_back("stn.245574_17.root");
			filenames.push_back("stn.245574_18.root");
			filenames.push_back("stn.245574_19.root");
			filenames.push_back("stn.245574_2.root");
			filenames.push_back("stn.245574_20.root");
			filenames.push_back("stn.245574_3.root");
			filenames.push_back("stn.245574_4.root");
			filenames.push_back("stn.245574_5.root");
			filenames.push_back("stn.245574_6.root");
			filenames.push_back("stn.245574_7.root");
			filenames.push_back("stn.245574_8.root");

		}
		std::cout << "" << std::endl;
		filename <<"IsrFsrPythiaPtHat"<< pt_hat <<"_MoreISRResults.root";

	} else if (dataset == 2) // less ISR
	{
		std::cout << "Using less ISR sample ";
		if (pt_hat ==20.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take5_PtHat20/Stntuples/";

			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_1");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_2");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_3");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_4");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_5");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_6");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_7");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessISR.root_8");

		} else if (pt_hat==60.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/lessISR/";
			filenames.push_back("stn.245319_1.root");
			filenames.push_back("stn.245319_10.root");
			filenames.push_back("stn.245319_12.root");
			filenames.push_back("stn.245319_13.root");
			filenames.push_back("stn.245319_14.root");
			filenames.push_back("stn.245319_15.root");
			filenames.push_back("stn.245319_16.root");
			filenames.push_back("stn.245319_18.root");
			filenames.push_back("stn.245319_19.root");
			filenames.push_back("stn.245319_2.root");
			filenames.push_back("stn.245319_20.root");
			filenames.push_back("stn.245319_3.root");
			filenames.push_back("stn.245319_4.root");
			filenames.push_back("stn.245319_5.root");
			filenames.push_back("stn.245319_6.root");
			filenames.push_back("stn.245319_7.root");
			filenames.push_back("stn.245319_8.root");
			filenames.push_back("stn.245319_9.root");
			filenames.push_back("stn.246164_11.root");
			filenames.push_back("stn.246165_17.root");


		} else if (pt_hat == 100.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/lessISR/";

			filenames.push_back("stn.245549_1.root");
			filenames.push_back("stn.245549_10.root");
			filenames.push_back("stn.245549_11.root");
			filenames.push_back("stn.245549_12.root");
			filenames.push_back("stn.245549_14.root");
			filenames.push_back("stn.245549_15.root");
			filenames.push_back("stn.245549_17.root");
			filenames.push_back("stn.245549_19.root");
			filenames.push_back("stn.245549_2.root");
			filenames.push_back("stn.245549_20.root");
			filenames.push_back("stn.245549_3.root");
			filenames.push_back("stn.245549_4.root");
			filenames.push_back("stn.245549_5.root");
			filenames.push_back("stn.245549_6.root");
			filenames.push_back("stn.245549_7.root");
			filenames.push_back("stn.245549_8.root");
			filenames.push_back("stn.245549_9.root");
			filenames.push_back("stn.246177_13.root");
			filenames.push_back("stn.246179_16.root");
			filenames.push_back("stn.246180_18.root");

		} else if (pt_hat == 160.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/lessISR/";
			filenames.push_back("stn.245562_1.root");
			filenames.push_back("stn.245562_11.root");
			filenames.push_back("stn.245562_12.root");
			filenames.push_back("stn.245562_17.root");
			filenames.push_back("stn.245562_19.root");
			filenames.push_back("stn.245562_2.root");
			filenames.push_back("stn.245562_20.root");
			filenames.push_back("stn.245562_3.root");
			filenames.push_back("stn.245562_4.root");
			filenames.push_back("stn.245562_8.root");
			filenames.push_back("stn.245562_9.root");

		} else if (pt_hat == 210.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/lessISR/";
			filenames.push_back("stn.245575_1.root");
			filenames.push_back("stn.245575_10.root");
			filenames.push_back("stn.245575_12.root");
			filenames.push_back("stn.245575_13.root");
			filenames.push_back("stn.245575_14.root");
			filenames.push_back("stn.245575_17.root");
			filenames.push_back("stn.245575_19.root");
			filenames.push_back("stn.245575_20.root");
			filenames.push_back("stn.245575_5.root");
			filenames.push_back("stn.245575_9.root");

		}

		std::cout << "" << std::endl;
		filename <<"IsrFsrPythiaPtHat"<< pt_hat <<"_LessISRResults.root";
	} else if (dataset == 3) // more FSR
	{
		std::cout << "Using more FSR sample";

		if (pt_hat ==20.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take5_PtHat20/Stntuples/";
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_1");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_2");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_3");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_4");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_5");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_6");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_7");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_moreFSR.root_8");

		} else if (pt_hat==60.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/moreFSR/";
			filenames.push_back("stn.245330_1.root");
			filenames.push_back("stn.245330_10.root");
			filenames.push_back("stn.245330_11.root");
			filenames.push_back("stn.245330_12.root");
			filenames.push_back("stn.245330_13.root");
			filenames.push_back("stn.245330_14.root");
			filenames.push_back("stn.245330_15.root");
			filenames.push_back("stn.245330_16.root");
			filenames.push_back("stn.245330_17.root");
			filenames.push_back("stn.245330_18.root");
			filenames.push_back("stn.245330_19.root");
			filenames.push_back("stn.245330_2.root");
			filenames.push_back("stn.245330_20.root");
			filenames.push_back("stn.245330_3.root");
			filenames.push_back("stn.245330_4.root");
			filenames.push_back("stn.245330_5.root");
			filenames.push_back("stn.245330_6.root");
			filenames.push_back("stn.245330_7.root");
			filenames.push_back("stn.245330_8.root");
			filenames.push_back("stn.245330_9.root");


		} else if (pt_hat == 100.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/moreFSR/";
			filenames.push_back("stn.245550_1.root");
			filenames.push_back("stn.245550_10.root");
			filenames.push_back("stn.245550_11.root");
			filenames.push_back("stn.245550_12.root");
			filenames.push_back("stn.245550_13.root");
			filenames.push_back("stn.245550_14.root");
			filenames.push_back("stn.245550_16.root");
			filenames.push_back("stn.245550_17.root");
			filenames.push_back("stn.245550_18.root");
			filenames.push_back("stn.245550_19.root");
			filenames.push_back("stn.245550_2.root");
			filenames.push_back("stn.245550_20.root");
			filenames.push_back("stn.245550_3.root");
			filenames.push_back("stn.245550_4.root");
			filenames.push_back("stn.245550_5.root");
			filenames.push_back("stn.245550_6.root");
			filenames.push_back("stn.245550_8.root");
			filenames.push_back("stn.245550_9.root");

		} else if (pt_hat == 160.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/moreFSR/";
			filenames.push_back("stn.245564_10.root");
			filenames.push_back("stn.245564_12.root");
			filenames.push_back("stn.245564_14.root");
			filenames.push_back("stn.245564_15.root");
			filenames.push_back("stn.245564_16.root");
			filenames.push_back("stn.245564_18.root");
			filenames.push_back("stn.245564_20.root");
			filenames.push_back("stn.245564_4.root");
			filenames.push_back("stn.245564_7.root");
			filenames.push_back("stn.245564_8.root");
			filenames.push_back("stn.245564_9.root");


		} else if (pt_hat == 210.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/moreFSR/";
			filenames.push_back("stn.245577_1.root");
			filenames.push_back("stn.245577_11.root");
			filenames.push_back("stn.245577_13.root");
			filenames.push_back("stn.245577_17.root");
			filenames.push_back("stn.245577_18.root");
			filenames.push_back("stn.245577_19.root");
			filenames.push_back("stn.245577_2.root");
			filenames.push_back("stn.245577_20.root");
			filenames.push_back("stn.245577_6.root");
			filenames.push_back("stn.245577_7.root");
			filenames.push_back("stn.245577_8.root");

		}

		std::cout << "" << std::endl;
		filename <<"IsrFsrPythiaPtHat"<< pt_hat <<"_MoreFSRResults.root";

	} else if (dataset == 4) // less FSR
	{
		std::cout << "Using less FSR sample";
		std::cout << " with minPt-hat " << pt_hat << " GeV";
		if (pt_hat ==20.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take5_PtHat20/Stntuples/";
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_1");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_2");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_3");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_4");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_5");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_6");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_7");
			filenames.push_back("stn.1429_PYTHIA_gammajet_ptHat20_HEPGonly_lessFSR.root_8");

		} else if (pt_hat==60.)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take6_PtHat60/lessFSR/";
			filenames.push_back("stn.245332_1.root");
			filenames.push_back("stn.245332_10.root");
			filenames.push_back("stn.245332_11.root");
			filenames.push_back("stn.245332_12.root");
			filenames.push_back("stn.245332_13.root");
			filenames.push_back("stn.245332_14.root");
			filenames.push_back("stn.245332_15.root");
			filenames.push_back("stn.245332_16.root");
			filenames.push_back("stn.245332_17.root");
			filenames.push_back("stn.245332_18.root");
			filenames.push_back("stn.245332_19.root");
			filenames.push_back("stn.245332_2.root");
			filenames.push_back("stn.245332_20.root");
			filenames.push_back("stn.245332_3.root");
			filenames.push_back("stn.245332_4.root");
			filenames.push_back("stn.245332_5.root");
			filenames.push_back("stn.245332_6.root");
			filenames.push_back("stn.245332_7.root");
			filenames.push_back("stn.245332_8.root");
			filenames.push_back("stn.245332_9.root");

		} else if (pt_hat == 100)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take8_PtHat100/lessFSR/";
			filenames.push_back("stn.245553_1.root");
			filenames.push_back("stn.245553_11.root");
			filenames.push_back("stn.245553_13.root");
			filenames.push_back("stn.245553_14.root");
			filenames.push_back("stn.245553_15.root");
			filenames.push_back("stn.245553_16.root");
			filenames.push_back("stn.245553_17.root");
			filenames.push_back("stn.245553_18.root");
			filenames.push_back("stn.245553_19.root");
			filenames.push_back("stn.245553_2.root");
			filenames.push_back("stn.245553_20.root");
			filenames.push_back("stn.245553_3.root");
			filenames.push_back("stn.245553_4.root");
			filenames.push_back("stn.245553_5.root");
			filenames.push_back("stn.245553_6.root");
			filenames.push_back("stn.245553_7.root");
			filenames.push_back("stn.245553_8.root");


		} else if (pt_hat == 160)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take9_PtHat160/lessFSR/";
			filenames.push_back("stn.245565_10.root");
			filenames.push_back("stn.245565_11.root");
			filenames.push_back("stn.245565_12.root");
			filenames.push_back("stn.245565_14.root");
			filenames.push_back("stn.245565_16.root");
			filenames.push_back("stn.245565_17.root");
			filenames.push_back("stn.245565_18.root");
			filenames.push_back("stn.245565_19.root");
			filenames.push_back("stn.245565_20.root");
			filenames.push_back("stn.245565_3.root");
			filenames.push_back("stn.245565_4.root");
			filenames.push_back("stn.245565_6.root");
			filenames.push_back("stn.245565_7.root");
			filenames.push_back("stn.245565_8.root");

		} else if (pt_hat == 210)
		{
			std::cout << " with minPt-hat " << pt_hat << " GeV";
			path = "root://nbay02.fnal.gov:5151//data/nbay02/b/samantha/MCSamples/Take10_PtHat210/lessFSR/";
			filenames.push_back("stn.245578_14.root");
			filenames.push_back("stn.245578_16.root");
			filenames.push_back("stn.245578_17.root");
			filenames.push_back("stn.245578_18.root");
			filenames.push_back("stn.245578_2.root");
			filenames.push_back("stn.245578_20.root");
			filenames.push_back("stn.245578_3.root");
			filenames.push_back("stn.245578_4.root");
			filenames.push_back("stn.245578_6.root");
			filenames.push_back("stn.245578_7.root");
			filenames.push_back("stn.245578_8.root");
			filenames.push_back("stn.245578_9.root");

		}
		std::cout << "" << std::endl;
		filename <<"IsrFsrPythiaPtHat"<< pt_hat <<"_LessFSRResults.root";
	}


	for (unsigned int i =0;i < filenames.size(); ++i)
	{
		std::stringstream file;
		file << path << filenames.at(i);
		ch->Add(file.str().c_str());
	}



	//if (debug)	ch->Print();
	
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
  	
	
	ap->SaveHist(filename.str().c_str(),2);
	delete pho;
	delete ch;
	std::cout << "/^^^^^^^^^^^^^ JOB SUMMARY ^^^^^^^^^^^^^^^^^^^^^^^^^^^^^^/" <<std::endl;
		std::cout << "::: Created ROOT file      => "<< filename.str() << std::endl;
	if (debug) 
	{
		std::cout << "::: DebugMode::Attaching " << filename.str() << " to a TBrowser." << std::endl;
		new TFile(filename.str().c_str()); 
	}
	gBenchmark->Show("metphoana_time");
	gROOT->ProcessLine(".!date");
	std::cout << "/********************************************************/" <<std::endl;
	
//	gROOT->ProcessLine(".q");
}

//run through all sample in one step
void IsrFsrHepgStudy(const float pt_hat) 
{
	for (int idataset =1; idataset<5; idataset++)
	{
		IsrFsrHepgStudy(idataset,pt_hat,0);
	}
}

void IsrFsrHepgStudy()
{
	IsrFsrHepgStudy(20); 
	IsrFsrHepgStudy(60); 
	IsrFsrHepgStudy(100); 
	IsrFsrHepgStudy(160); 
	IsrFsrHepgStudy(210); 
}


