{

TH1F* h = new TH1F("hist","hist",5,0,5);

hist->Fill(-1);
hist->Fill(-1);
hist->Fill(1);
hist->Fill(1);
hist->Fill(3);
hist->Fill(31);
hist->Fill(31);
hist->Fill(31);
hist->Print("all");
std::cout << "int = " << hist->Integral() << std::endl;
std::cout << "int = " << hist->Integral("width") << std::endl;
double w =0;
for (int bin=1; bin <= hist->GetNbinsX(); bin++) 
{
	std::cout << "bin/val " << bin << "/" << hist->GetBinContent(bin) << std::endl;
	w+= hist->GetBinContent(bin) * hist->GetBinWidth(bin);
}
std::cout << " w = " << w << std::endl;
}
