//------------------------------------------------------------------------------
//  rootlogon.C: a sample ROOT logon macro allowing use of ROOT script 
//               compiler in CDF RunII environment. The name of this macro file
//               is defined by the .rootrc file
//
//  USESHLIBS variable has to be set to build Stntuple libraries locally:  
//
//  setenv USESHLIBS 1
//
//  Feb 12 2001 P.Murat
//------------------------------------------------------------------------------
{
#include <iomanip.h>
#include <time.h>
#include <iostream>             


	std::cout << "Loading ~/stn_rel/rootlogon.C" << std::endl;

	// the line below tells ROOT script compiler
	// where to look for the include files

	gSystem->SetIncludePath(" -I./include -I$CDFSOFT2_DIR/include");

	// Setting up SAM to use with ROOT
	gSystem->Load("libsam_cpp_api.so");
	gSystem->Load("libdiskcache_samroot.so"); 

	// load in ROOT physics vectors and event
	// generator libraries


	gSystem->Load("$ROOTSYS/lib/libPhysics.so");
	gSystem->Load("$ROOTSYS/lib/libEG.so");
	gSystem->Load("$ROOTSYS/lib/libDCache.so");

	// load a script with the macros
	char command[200];

	sprintf(command,"%s/Stntuple/scripts/global_init.C",
			gSystem->Getenv("CDFSOFT2_DIR"));

	gInterpreter->LoadMacro(command);

	// STNTUPLE shared libraries are assumed to be
	// built in the private test release area with 
	// USESHLIBS environment variable set 
	// we always need libStntuple_loop, but the
	// other 2 libs should be loaded in only if
	// we're running bare root

	const char* exec_name = gApplication->Argv(0);

	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_base.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_obj.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_alg.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_loop.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_photon.so");

	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libZMutility.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libExceptions.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libJetUser_forRoot.so");

	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_photon.so"); // 3/28/07 sam: need this and PhoBC for CES/CPR method
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libPhoBC.so"); // 3/23/07 sam: for flat ntuples //for Highlevel EM objects
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libPho.so");  // for ray's package of GoodRun
	//gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/shlib/$BFARCH/libDiPho.so"); 

	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_geom.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_val.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libStntuple_ana.so");

	//std::cout << "loading sam's stuff in ~/stn_rel/rootlogon.C" << std::endl;
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libSamantha_utils.so");
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libSamantha_obj.so");

	// since Pho2Jets is tangled with TPDFReweight, need to load pdf tools stuff before libPho.
	// better use the driver script to do PDFsyst stuff (or set
	// the env vars manually and run). look under samantha/pdftools/
	std::cout << "Loading ./samantha/pdftools/libpdftools.so" << std::endl;
	gSystem->Load("./samantha/pdftools/libpdftools.so");
	std::cout << "Compiling/Loading ./samantha/pdftools/TPDFReweight.cc" << std::endl;
	gSystem->CompileMacro("./samantha/pdftools/TPDFReweight.cc","k");
	std::cout << "Loading libSamantha_Pho.so" << std::endl;
	gSystem->Load("$BASERELEASE/shlib/Linux2_SL-GCC_3_4/libSamantha_Pho.so");
	//gSystem->Load("./libSamantha_MetModel.so");
	//std::cout << "loading slave" << std::endl;
	//gSystem->Load("./cafFlatNtupleMaker_C.so");
	//gSystem->Load("./runPho2Jets_caf_C.so");
	//gSystem->Load("./runPhoJetCount_C.so");

	//TGaxis::SetMaxDigits(3); // sam: limits the number of digits shown in axis labels

	// print overflows/underflows in the stat box
	gStyle->SetOptStat(11111111);
	// print fit results in the stat box
	gStyle->SetOptFit(1110);
	// this line reports the process ID which simplifies
	// debugging

	//TAuthenticate::SetGlobalUser(gSystem->Getenv("USER"));
	 //TAuthenticate::SetGlobalPasswd("aa");
	 //std::cout << "USER= " << 	TAuthenticate::GetGlobalUser() << std::endl;
	 //std::cout << "PASS= " << 	TAuthenticate::GetPromptUser() << std::endl;
	// TAuthenticate::PromptPasswd("\n");
	 //TAuthenticate::PromptUser
	 //TAuthenticate::ProcessUser(gSystem->Getenv("USER"));
	//gInterpreter->ProcessLine(".! ps | grep root");

	//canvas settings
	gStyle->SetCanvasColor (10);
	gStyle->SetCanvasBorderSize (0);
	gStyle->SetCanvasBorderMode (0);

	gStyle->SetPadColor (10);
	gStyle->SetFillColor (10);
	gStyle->SetTitleFillColor (10);
	gStyle->SetTitleBorderSize (0);
	gStyle->SetStatColor (10);
	gStyle->SetStatBorderSize (1);

}
