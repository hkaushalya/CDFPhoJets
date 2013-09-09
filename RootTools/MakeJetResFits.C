// ROOT macro to gererate Jet Resolution
// hists

void MakeJetResFits()
{
	gSystem->CompileMacro("~/samantha/RootTools/JetResFitPerEPerEtaBin.C","f");

	for (int e = 0; e < 3; ++e)
	{
		for (int eta = 0; eta <= 2; eta++)
		{
			std::cout << "E,Eta = " << e << "\t"<< eta << std::endl;
			//std::cout << "E,Eta = " << e << "\t"<< 1 << std::endl;
			//std::cout << "E,Eta = " << 2 << "\t"<< 1 << std::endl;
			JetResFitPerEPerEtaBin(e,eta);
			//JetResFitPerEPerEtaBin(e,1);
		}
	}


}
