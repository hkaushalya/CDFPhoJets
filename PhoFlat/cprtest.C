{	
	TChain *ch = 0;
	std::vector<std::string> files;

	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0scd_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sdd_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sed_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0see_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0seh_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze0sej_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze1s6d_NoJetVtxCuts.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze1sad_NoJetVtxCuts_1.root");
	files.push_back("~/STUPLES2/StupleV7_MCPythia_Zee_ze1sad_NoJetVtxCuts.root");

	for (int i=0; i<files.size() ; ++i)
	{
		std::cout << "i = " << i << std::endl;
		ch = new TChain("Stuple");
		ch->Add(files.at(i).c_str());
		new TCanvas();
		ch->Draw("pho_CprWgt");

	}
}
