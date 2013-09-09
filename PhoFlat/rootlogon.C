{
#include <iostream>
	std::cout << "Loading " << gSystem->pwd() << "/rootlogon.C" << std::endl;
	std::string path = "PhoFlat";
	std::string include = gSystem->GetIncludePath ();

	// this is to detect and fix some stupid inconsistency between
	// different root versions.
	std::string prefix = "-I", postfix = "";
	if (include.find ("-I\"") != std::string::npos)
	{
		prefix = "-I\"";
		postfix = "\"";
	};

	include += " " + prefix + path + postfix;
	gSystem->SetIncludePath (include.c_str());
	gSystem->Load("$ROOTSYS/lib/libPhysics.so");
	//gSystem->Load("$ROOTSYS/lib/libEG.so");
	//gStyle->SetOptStat("neoui");

//my general canvas settings
  gStyle->SetOptStat(11111111);
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);
				 
	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);
	//gStyle->SetTitleOffset(0.1);
	gStyle->SetTitleColor(kBlue,"XYZ");
	gStyle->SetTitleColor(kBlue);
	gStyle->SetTitleSize(0.04,"XYZ");
//	gStyle->SetTitleSize(0.1);

	//gSystem->Load("../RootTools/CommonTools_cc.so");

};
