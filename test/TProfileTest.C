#include<iostream>
#include "TProfile.h"


//how do you setbincontent and error for TProfiles.
//if you just set them and draw it,nothing shows up in the hist!
//you have to setbinentries to some values to see things!

void TProfileTest()
{

	TProfile *hist = new TProfile("tprof","test" , 10,0,100);

	//hist->SetBinContent(1, 100);
	//hist->SetBinError(1,10);
	//hist->SetBinEntries(1,1);
	for (int i=1; i<=5; i++)
	{
		hist->Fill(i,i);
		float val = hist->GetBinContent(1);
		float err = hist->GetBinError(1);
		float entries = hist->GetBinEntries(1);
		std::cout << "i["<< i << "] val[" << val 
					<< "] err[" << err << "] entries[" << entries << "]" << std::endl;
//		hist->SetBinEntries(i,2);
//		std::cout << "\ti["<< i << "] val[" << val 
//					<< "] err[" << err << "] entries[" << entries << "]" << std::endl;
		//hist->Fill(i+0.1);
		//hist->Fill(i+0.2);
	}
	hist->Draw();
	hist->Print();
}

void RebinTest()
{
	TH1F* hist = new TH1F("hist","rebin test",10,0,10);
	hist->Sumw2();
	
	std::cout << "i \t val1/err1 \t val2/err2" << std::endl;
	for (int i=0;i<15;++i)
	{
		float fval = 0.1 * i;
		hist->Fill(fval);
		float val1 = hist->GetBinContent(1);
		float err1 = hist->GetBinError(1);
		float val2 = hist->GetBinContent(2);
		float err2 = hist->GetBinError(2);
		std::cout << i << "\t" << val1 << " / " << err1 
				<< " \t " << val2 << " / " << err2 << std::endl;
	}
	hist->Draw("E");

	hist->Rebin(2);
	
		float val1 = hist->GetBinContent(1);
		float err1 = hist->GetBinError(1);
		float val2 = hist->GetBinContent(2);
		float err2 = hist->GetBinError(2);
		std::cout << "\n\t" << val1 << " / " << err1 
				<< " \t " << val2 << " / " << err2 << std::endl;

}
