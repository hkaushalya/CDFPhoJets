{
	TChain ch("STNTUPLE");
	//TList *filelist = TFileInfo::CreateListMatching("root://remotehost/path/to/files/*.root");
//	TList *filelist = TFileInfo::CreateListMatching("root://holmium.fnal.gov:5151//cdf/scratch/samantha/p27_calibexe_stn/cpr2stntuple_gb.*.root");
//	if (filelist) {
//	   ch.AddFileInfoList(filelist);
//		   delete filelist;
//	}


	dir = gSystem->OpenDirectory("root://holmium.fnal.gov:5151//cdf/scratch/samantha/p27_calibexe_stn/");
	char *ent;
	while ((ent = gSystem->GetDirEntry(dir))) {
		TString fn = Form("root://holmium.fnal.gov:5151//cdf/scratch/samantha/p27_calibexe_stn/%s", ent);
		std::cout << "file = " << fn << std::endl;
		if (fn.EndsWith(".root")) {
			FileStat_t st;
			if (!gSystem->GetPathInfo(fn, st) && R_ISREG(st.fMode))
				ch.Add(fn);
		}
	}
	gSystem->FreeDirectory(dir); 

	std::cout << "Total Entries = " << ch.GetEntries() << std::endl;
} 
